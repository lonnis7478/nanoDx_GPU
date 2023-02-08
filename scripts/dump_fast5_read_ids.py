import sys
from ont_fast5_api.fast5_interface import get_fast5_file

f = sys.argv[1]

with get_fast5_file(f, mode="r") as f5:
  for read in f5.get_reads():
    print(read.get_read_id())

