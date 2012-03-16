#!/tools/opt/python/python2.7.2/bin/python2.7

from mlst_modules import rm_files
import subprocess
import os

rm_files(['run_mlst_velvet.*', 'run_mlst_interleave.*', '*.interleaved', 'pbsjob.tmp*', 'run_mlst_postprocess.*', 'run_mlst_trim.*', 'run_velvet*', '*.sh'])

if os.path.exists('trimmed'):
   subprocess.call('rm -r trimmed/', shell=True)
