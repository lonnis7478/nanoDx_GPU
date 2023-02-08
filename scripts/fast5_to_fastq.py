import sys
from ont_fast5_api.fast5_interface import get_fast5_file
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools

f = sys.argv[1]

with get_fast5_file(f, mode="r") as f5:
  for read in f5.get_reads():
    with Basecall1DTools(read, group_name=read.get_latest_analysis('Basecall_1D')) as basecall:
      fq = basecall.get_called_sequence('template', fastq = True)
      print(fq)
