#!/tools/opt/python/python2.7.2/bin/python2.7

# NB: Need to copy convertFastq2Fasta.py to /home/panfs/cbs/projects/cge/servers/MLST/assembly/
# also mlst_home is defined as /panvol1/simon/bin/mlst/ so it is actually not running the scripts at the server dir

import argparse
import logging
import os
from mlst_classes import Semaphore_file
import mlst_modules

def required_nargs_abspath(min,max):
   '''Enforces input to nargs to be between min and max long and sets abspath for input 1:'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not min<=len(value)<=max:
            msg='argument "{f}" requires between {min} and {max} arguments'.format(
               f=self.dest,min=min,max=max)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
         import os
         if type(value) == str:
            f_abs = os.path.abspath(value)
            setattr(args, self.dest, f_abs)
         elif type(value) == list:
            new_list = [value[0]]
            for f in value[1:]:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, value)
   return RequiredInterval

def required_nargs(min,max):
   '''Enforces input to nargs to be between min and max long'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not min<=len(value)<=max:
            msg='argument "{f}" requires between {min} and {max} arguments'.format(
               f=self.dest,min=min,max=max)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
   return RequiredInterval


def set_abspath():
   '''Returns absolute path of file to argparse'''
   class SetAbspath(argparse.Action):
      def __call__(self, parser, args, filenames, option_string=None):
         import os
         if type(filenames) == str:
            f_abs = os.path.abspath(filenames)
            setattr(args, self.dest, f_abs)
         elif type(filenames) == list:
            new_list = []
            for f in filenames:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, filenames)
   return SetAbspath

def set_abspath_velvet():
   '''Returns absolute path of file to argparse'''
   class SetAbspath(argparse.Action):
      def __call__(self, parser, args, filenames, option_string=None):
         import os
         if type(filenames) == str:
            f_abs = os.path.abspath(filenames)
            setattr(args, self.dest, f_abs)
         elif type(filenames) == list:
            new_list = [filenames[0]]
            for f in filenames[1:]:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, filenames)
   return SetAbspath


def create_jobs(prog, args, sfile, logger):
   '''Create an msub command from prog and args'''
   
   import subprocess
   import os
   import mlst_modules
   
   paths = mlst_modules.setSystem()
   
   # create commands
   msub = '%sxmsub -d %s -l nodes=1:ppn=1,mem=256mb,walltime=172800,partition=%s -q %s -r y -N run_%s -O run_%s.out -E run_%s.err -t' % (paths['mlst_home'], os.getcwd(), args.partition, args.q, args.assembler, args.assembler, args.assembler)
   cmd = [msub, prog]
   
   # create parameters for assembler
   for key, value in vars(args).items():
      # special cases
      if key == 'assembler':
         continue
      if value == None:
         continue
      if key == 'wait':
         cmd.append('--sfile %s' % sfile)
         continue
      if type(value) == bool:
         if value == True: cmd.append('--%s' % key)
         continue
      if key == 'sample':
         cmd.append('--sample None')
         continue
      if key == 'add_solid':
         cmd.append('--add_solid "%s"' % value)
         continue
      if key == 'add_velveth':
         cmd.append('--add_velveth "%s"' % value)
         continue
      if key == 'add_velvetg':
         cmd.append('--add_velvetg "%s"' % value)
         continue
      
      # key-value paramters
      cmd.append('--%s' %key)
      if type(value) == list:
         cmd.append(' '.join(value))
      elif type(value) == str or type(value) == int or type(value) == float:
         cmd.append('%s' %value)
      else:
         raise ValueError('%s, %s is a %s, should be either list, string or int' % (key, value, type(value)))
   
   # submit job
   cmd = ' '.join(cmd)
   logger.info(cmd)
   job = subprocess.check_output(cmd, shell=True)
   job = job.strip('\n')
   print 'Jobs are spawned by: %s' % job
   


parser = argparse.ArgumentParser(prog='mlst_wrapper_denovo_2.0.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=25),
                                 description='''Wrapper to start denovo assembly for MLST''', 
                                 usage='%(prog)s module [options]')

# general (parent parser)
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--sample', help='name of run and output directory [mlst_run]', default='mlst_run')
parent_parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parent_parser.add_argument('--m', help='memory needed for assembly [7gb]', default='7gb')
parent_parser.add_argument('--partition', help='partition to run on (cge-cluster, uv) [cge-cluster]', default='cge-cluster')
parent_parser.add_argument('--q', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parent_parser.add_argument('--log', help='log level [INFO]', default='info')
parent_parser.add_argument('--wait', help='wait for assembly, ONLY www [False]', default=False, action='store_true')

# create subparsers
subparsers = parser.add_subparsers(dest='assembler')

# velvet
parser_velvet = subparsers.add_parser('velvet', help='Run velvet denovo assembly', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=110), usage='mlst_wrapper_denovo_2.0.py velvet [options]', description=
   '''Run Velvet denovo assembly, add_velvetg/add_velveth has to be added with quotes, eg: add_velvetg "-very_clean yes". NB: Trimming does not work on already interleaved files, This version of the wrapper only takes fasta and fastq as input\nDefault ksizes (auto) will run a range awer45 (read_length/3..read_length/3*2) with step 4 between \n''')
parser_velvet.add_argument('--short', help='short reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--shortPaired', help='short paired reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--short2', help='short2 reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--shortPaired2', help='short paired2 reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--long', help='long reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--longPaired', help='long paired reads', nargs='+', action=set_abspath_velvet())
parser_velvet.add_argument('--ksizes', help='kmers to run assemblies for (single (m) or m M s (min, max, step)) [auto]', nargs='+', default=['auto'])
parser_velvet.add_argument('--outpath', help='name of run, also output dir [assembly]', default='assembly')
parser_velvet.add_argument('--trim', help='should input files be trimmed (illumina only) [False]', default=False, action='store_true')
parser_velvet.add_argument('--min_contig_lgth', help='mininum length to report contig [100]', default=100, type=int)
parser_velvet.add_argument('--cov_cutoff', help='coverage cutoff for removal of low coverage (float) [None]', default=None, type=float)
parser_velvet.add_argument('--exp_cov', help='Expected mean coverage (None, float, auto) [auto]', default='auto')
parser_velvet.add_argument('--ins_length', help='insert size (reads included) [None]', default=None, type=int)
parser_velvet.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
parser_velvet.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)

# newbler
parser_newbler = subparsers.add_parser('newbler', help='Run newbler denovo assembly', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=110), usage='mlst_wrapper_denovo_2.0.py newbler [options]', description=
   '''Run Newbler (2.6) denovo assembly. Format can be: fasta, fastq, sff. If both fasta and qual files should be used for assembly only add the fasta file and make sure the qual file is in the same path with same name except .qual ending''')
parser_newbler.add_argument('--se', help='input sff/fasta/fastq files', nargs='+', action=set_abspath())
parser_newbler.add_argument('--pe', help='input paired end sff/fasta/fastq files', nargs='+', default=None, action=set_abspath())
parser_newbler.add_argument('--outpath', help='assembly output dir [denovo]', default='denovo')

# solid
parser_solid = subparsers.add_parser('solid', help='Run solid denovo assembly (uses velvet)', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=110), usage='mlst_wrapper_denovo_2.0.py solid [options]', description=
   ''' ''')
parser_solid.add_argument('--se', help='input single end csfasta file (F3, F3q)', nargs='+', action=set_abspath())
parser_solid.add_argument('--pe', help='input paired end csfasta files (F3, F3q, F5, F5q) ', nargs='+', default=None, action=set_abspath())
parser_solid.add_argument('--mp', help='input mate pair csfasta files (F3, F3q, R3, R3q)', nargs='+', default=None, action=set_abspath())
parser_solid.add_argument('--rf', help='input expected length of genome in bp', type=int, required=True)
parser_solid.add_argument('--ins_length', help='estimate of mate/paired end insert length eg. (1200/170)', type=int, required=True)
parser_solid.add_argument('--ins_length_sd', help='estimate of mate/paired end insert length eg. (300/30)', type=int, required=True)
parser_solid.add_argument('--add_solid', help='additional parameters to solid assembler', default=None)

args = parser.parse_args()
#args = parser.parse_args('velvet --shortPaired Kleb-10-213361_2_1_sequence.txt Kleb-10-213361_2_2_sequence.txt --ksizes 41 55 4 --trim'.split())
#args = parser.parse_args('velvet --shortPaired Kleb-10-213361_2.interleaved.fastq --trim --sample Kleb_auto'.split())
#args = parser.parse_args('velvet --short 110601_I238_FCB067HABXX_L3_ESCqslRAADIAAPEI-2_1.fq --ksizes 45 75 --sample E_coli_TY2482_illumina --trim'.split())
#args = parser.parse_args('velvet --shortPaired test_kleb_1.fq test_kleb_2.fq --ksizes 41 55 4 --sample test_kleb --cov_cutoff 8'.split())
#args = parser.parse_args('newbler --se life_unimuenster_sff/*.sff --sample test_newbler --wait'.split())
#args = parser.parse_args('solid --mp ecoli_600x_F3.csfasta ecoli_600x_F3.qual ecoli_600x_R3.csfasta ecoli_600x_R3.qual --rf 5000000 --ins_length 1300 --ins_length_sd 300 --m 7gb --sample solid_test --wait'.split())
#args = parser.parse_args('velvet --shortPaired test_kleb_1.fq test_kleb_2.fq --ksizes 41 55 4'.split())

# set pythonpath
os.environ['PYTHONPATH'] = '/panvol1/simon/lib/python/:/home/panfs/cbs/projects/cge/servers/MLST/assembly/'
paths = mlst_modules.setSystem()

# If working dir is given, create and move to working directory else run where program is invoked
if args.sample:
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   #os.chmod(args.sample, 0777)
   os.chdir(args.sample)
else:
   pass

# create log dir
if not os.path.exists('log'):
   os.makedirs('log')

# set logging
logger = logging.getLogger('mlst_wrapper_denovo.py')
hdlr = logging.FileHandler('mlst_wrapper_denovo.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# toggle waiting
if args.wait:
   import random
   import string
   rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
   sfile = 'semaphores/denovo%s' % rand
else:
   sfile = None

# create jobs
if args.assembler == 'velvet':
   prog = '%smlst_denovo_velvet.py' % paths['mlst_home']
   create_jobs(prog, args, sfile, logger)
elif args.assembler == 'newbler':
   prog = '%smlst_denovo_newbler.py' % paths['mlst_home']
   create_jobs(prog, args, sfile, logger)
elif args.assembler == 'solid':
   prog = '%smlst_denovo_solid.py' % paths['mlst_home']
   create_jobs(prog, args, sfile, logger)

# if waiting
if args.wait:
   s = Semaphore_file(sfile, 60, 2*86400)
   s.wait()
