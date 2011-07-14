#!/panvol1/simon/bin/python

from mlst_modules import rm_files
import subprocess
import os

rm_files(['run_mlst_velveth.*', 'run_mlst_velvetg.*', 'run_mlst_interleave.*', '*.interleaved', 'pbsjob.tmp*', 'run_mlst_velvetaccept.*', 'run_mlst_trim.*'])

if os.path.exists('trimmed'):
   subprocess.call('rm -r trimmed/', shell=True)
