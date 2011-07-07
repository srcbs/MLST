#!/panvol1/simon/bin/python2.7

import argparse
import logging
import os
import time
import subprocess
import pipelinemod
import re
import sys

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

parser = argparse.ArgumentParser(prog='mlst_wrapper_denovo_2.0.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=25),
                                 description='''Wrapper to start denovo assembly for MLST''', 
                                 usage='%(prog)s module [options]')

# general (parent parser)
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--sample', help='name of run and output directory', default=None)
parent_parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parent_parser.add_argument('--m', help='memory needed for assembly [2gb]', default='2gb')
parent_parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cge]', default='cge')
parent_parser.add_argument('--log', help='log level [INFO]', default='info')

# create subparsers
subparsers = parser.add_subparsers(dest='assembler')

parser_velvet = subparsers.add_parser('velvet', help='Run velvet denovo assembly', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=110), usage='mlst_wrapper_denovo_2.0.py velvet [options]')
parser_velvet.add_argument('--short', help='input read format and short reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--shortPaired', help='input read format and short paired reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--short2', help='input read format and short2 reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--shortPaired2', help='input read format and short paired2 reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--long', help='input read format and long reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--longPaired', help='input read format and long paired reads', nargs='+', action=required_nargs(0,3))
parser_velvet.add_argument('--ksizes', help='kmers to run assemblies for (single no or range) [33]', nargs='+', default=[33])
parser_velvet.add_argument('--outpath', help='name of run, also output dir [velvet_assembly]', default='velvet_assembly')
parser_velvet.add_argument('--min_contig_lgth', help='mininum length to report contig [100]', default=100, type=int)
parser_velvet.add_argument('--cov_cut', help='coverage cutoff for removal of low coverage (float) [None]', default=None, type=float)
parser_velvet.add_argument('--exp_cov', help='Expected mean coverage (None, float, auto) [auto]', default='auto')
parser_velvet.add_argument('--ins_length', help='insert size (reads included) [None]', default=None, type=int)
parser_velvet.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
parser_velvet.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)

parser_newbler = subparsers.add_parser('newbler', help='Run newbler denovo assembly', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=100), usage='mlst_wrapper_denovo_2.0.py newbler [options]')
parser_newbler.add_argument('--i', help='single end input files (454:sff/fa+qual/fastq/fa)', nargs='+', action='append', default=None)
parser_newbler.add_argument('--pe', help='if using paired reads, paired reads file [None]', nargs='+', action='append', default=None)

args = parser.parse_args()

# If working dir is given, create and move to working directory else run where program is invoked
if args.sample:
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   os.chmod(args.sample, 0777)
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


if args.assembler == 'velvet':
   from mlst_denovo_velvet import *
   start_assembly(args)
elif args.assembler == 'newbler':
   pass

   
