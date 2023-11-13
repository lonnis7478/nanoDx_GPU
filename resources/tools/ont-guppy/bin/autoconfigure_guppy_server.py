#! /usr/bin/env python3

import argparse
import logging
import os
from pathlib import Path
import platform
import queue
import signal
import subprocess
import sys
import threading
import time

logger = logging.getLogger("autoconfigure_guppy_server")
# Set up a default NullHandler in case we don't end up using another one
# Taken from http://docs.python-guide.org/en/latest/writing/logging/
logger.addHandler(logging.NullHandler())

# Timeouts are in seconds
AUTO_PORT_ID_TIMEOUT = 10
SERVER_SHUTDOWN_TIMEOUT = 3


class MKFormatter(logging.Formatter):
    """ A formatter that outputs log records matching MinKNOW's standard, which
    is this::

    <date_time>    <level>: <message_tag> (<logger_name>)
        <extra_data_key>: <extra_data_value>

    Where the date_time field is *almost* like Python's default asctime, but
    uses decimal milliseconds instead of a comma, and <level> is right
    justified eight spaces. We'll just use a normal message instead of
    message_tag.

    """
    format_string = '%(asctime)s.%(msecs)06d %(levelname)8s: %(message)s (%(name)s)'

    def __init__(self):
        super().__init__(fmt=self.format_string, datefmt='%Y-%m-%d %H:%M:%S')

    def format(self, record):
        """ Take a standard log record, check for an "extra_log_data"
        dictionary provided as part of the record, and, if the field is present,
        output those key-value pairs as extra lines in the log.

        Example output should be something like this::

            2021-03-23 14:36:31.348083    INFO: mk_manager_starting (mk_manager)
                guppy_version: 5.0.0+e4819c3ab
                config_version: Unknown

        """
        message = super().format(record)
        if (hasattr(record, 'extra_log_data') and
                type(record.extra_log_data) is dict):
            for key, value in record.extra_log_data.items():
                extra_info = "\n    {}: {}".format(key, value)
                message += extra_info
        return message


def _set_up_logging() -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(MKFormatter())
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.info("logging_initialised")


def _check_server_executable(path_to_server: Path) -> None:
    if not path_to_server.exists():
        log_details = {"expected_server_path": path_to_server}
        logger.error("cannot_find_server_exe", extra={"extra_log_data":
                                                      log_details})
        sys.exit(1)
    else:
        # We can't resolve a non-existant path, hence why we have two
        # log_details sections.
        log_details = {"expected_server_path": path_to_server.resolve()}
        logger.info("found_server_exe", extra={"extra_log_data":
                                               log_details})


def _get_available_threads() -> int:
    total_threads = os.cpu_count()
    # MinKNOW needs two threads for itself, if at all possible
    return max(total_threads - 2, 1)


def _construct_server_args(server_exe: Path, log_dir: Path,
                           available_threads: int, port: int,
                           auto_port: bool, use_tcp: bool) -> list:
    args = [str(server_exe),
            "--log_path", str(log_dir),
            "--config", "dna_r9.4.1_450bps_fast.cfg",
            "--num_callers", "1",
            "--cpu_threads_per_caller", str(available_threads),
            "--ipc_threads", "3"]
    if auto_port:
        logger.info("using_auto_port")
        args += ["--port", "auto"]
    else:
        args += ["--port", str(port)]
    if use_tcp:
        args += ["--use_tcp"]
    return args


def _extract_port_from(process):
    # Pipe process' stdout through a queue on a separate thread so we can make
    # it non-blocking. If we don't do this then stdout.readline() will block
    # when there's no stdout left, and that could end up causing the script to
    # hang.
    # https://stackoverflow.com/a/4896288
    def enqueue_output(out, queue, stop):
        for line in iter(out.readline, ''):
            if stop.is_set():
                break
            queue.put(line)

    stdout_queue = queue.Queue()
    stop_event = threading.Event()
    stdout_thread = threading.Thread(target=enqueue_output,
                                     args=(process.stdout, stdout_queue,
                                           stop_event))
    # This will cause the thread to be shut down with the script itself
    stdout_thread.daemon = True
    stdout_thread.start()

    port = None
    now = time.time()
    while(port is None and
          process.poll() is None and
          time.time() - now < AUTO_PORT_ID_TIMEOUT):
        try:
            line = stdout_queue.get_nowait().strip()
        except queue.Empty:
            line = ''
        # We're looking for:
        #   Starting server on port: <port>
        if line.startswith("Starting server on port:"):
            port = line.split()[-1]
        time.sleep(0.1)
    stop_event.set()
    return port


def _launch_server(server_exe: Path, log_dir: Path,
                   available_threads: int, port: int,
                   auto_port: bool, use_tcp: bool) -> int:

    # We'll use this to intercept ctrl-c and try and shut everything down
    # cleanly.
    # https://stackoverflow.com/a/43787607
    class SIGINT_handler():
        def __init__(self):
            self.SIGINT = False

        def signal_handler(self, _0, _1):
            self.SIGINT = True

    args = _construct_server_args(server_exe, log_dir, available_threads, port,
                                  auto_port, use_tcp)
    logger.info("launching_server",
                extra={"extra_log_data": {"args": args}})

    handler = SIGINT_handler()
    signal.signal(signal.SIGINT, handler.signal_handler)

    process = subprocess.Popen(args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)

    try:
        if auto_port:
            port = _extract_port_from(process)
            if port is None:
                logger.error("could_not_identify_auto_port")
                handler.SIGINT = True
            else:
                logger.info("server_running_on_port",
                            extra={"extra_log_data": {"port": port}})
        while process.poll() is None and not handler.SIGINT:
            time.sleep(1)
    finally:
        logger.info("shutting_down_server")
        if platform.system() == "Windows":
            os.kill(process.pid, signal.CTRL_C_EVENT)
        else:
            process.send_signal(signal.SIGINT)
        # Give the server a few seconds to shut down cleanly, then get more
        # aggressive.
        now = time.time()
        while(process.poll() is None and
              time.time() - now < SERVER_SHUTDOWN_TIMEOUT):
            time.sleep(1)
        if process.poll() is None:
            process.terminate()

    result_details = {"returncode": process.returncode,
                      "stdout": [x.strip() for x in process.stdout.readlines()],
                      "stderr": [x.strip() for x in process.stderr.readlines()]}
    if process.returncode == 0:
        logger.info("server_exited_gracefully",
                    extra={"extra_log_data": result_details})
    else:
        logger.error("server_terminated_with_error",
                     extra={"extra_log_data": result_details})
    return process.returncode


def main(server_exe: Path, log_dir: Path, disable_logging: bool,
         port=5555, max_threads=None, auto_port: bool=False, use_tcp: bool=False) -> None:
    if not disable_logging:
        _set_up_logging()
    _check_server_executable(server_exe)
    available_threads = _get_available_threads()
    if max_threads is not None:
        available_threads = min(available_threads, max(max_threads, 1))
    returncode = _launch_server(server_exe, log_dir, available_threads, port,
                                auto_port, use_tcp)
    sys.exit(returncode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatically configure a CPU-only version of "
        "guppy_basecall_server, based on the resources available.")
    parser.add_argument("--server_exe",
                        help="Path to the guppy_basecall_server executable. If "
                        "unspecified then this directory will be searched.",
                        default=None,
                        type=Path)
    parser.add_argument("--log_path",
                        help="Directory to output logs to. If unspecified then "
                        "a directory will be created in this one.",
                        default=None,
                        type=Path)
    parser.add_argument("--disable_logging",
                        help="Turn off logging.",
                        action="store_true")
    parser.add_argument("--port",
                        help="Port for the server to use. Defaults to 5555.",
                        default=5555,
                        type=int)
    parser.add_argument("--max_threads",
                        help="Maximum threads for the server to use. Will "
                        "always be at least one, and will be capped at "
                        "nproc - 2",
                        default=None,
                        type=int)
    # Instruct guppy to choose its own port
    parser.add_argument("--auto_port",
                        help=argparse.SUPPRESS,
                        action="store_true")
    parser.add_argument("--use_tcp",
                        help="Prefer tcp connection over Unix socket files. "
                        "Does nothing on windows.",
                        action="store_true")
    args = parser.parse_args()
    if args.server_exe is None:
        args.server_exe = Path.joinpath(Path(__file__).parent,
                                        "guppy_basecall_server")
        system_name = platform.system()
        if "Windows" in system_name or "WIN" in system_name:
            args.server_exe = args.server_exe.with_suffix('.exe')
    if args.log_path is None:
        args.log_path = Path.joinpath(Path(__file__).parent,
                                     "guppy_logs")
    main(args.server_exe, args.log_path, args.disable_logging, args.port,
         args.max_threads, args.auto_port, args.use_tcp)
