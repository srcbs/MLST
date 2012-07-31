#!/usr/bin/env python2.7

import argparse
import os

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


def newbler(args):
   '''Creating newbler calls'''

   def convert_fastq(args, paths):
      '''If input is fastq convert to fasta+qual'''
      
      cmds = []
      se = []
      pe = []
      # identify file inputs
      if args.se:
         se_ftypes = map(mlst_modules.set_filetype, args.se)
         for i,f in enumerate(args.se):
            if se_ftypes[i] == 'fastq':
               cmds.append('%smlst_fastq2fastaqual.py --i %s --p %s' % (paths['mlst_home'], f, os.path.split(f)[1]))
               se.append(os.path.split(f)[1]+'.fasta')
            elif se_ftypes[i] == 'fasta':
               if f.endswith('.gz'):
                  fnew = os.path.splitext(os.path.split(f)[1])[0]
                  cmds.append('''%spigz -dc -p %s %s > %s''' % (paths['pigz_home'], args.n, f, fnew))
                  se.append(fnew)
                  
                  # look for qual file (dont add to path because newbler will pick it up)
                  possible_qual = os.path.splitext(os.path.splitext(f)[0])[0] + '.qual.gz'
                  if os.path.exists(possible_qual):
                     qnew = os.path.split(os.path.splitext(os.path.splitext(f)[0])[0] + '.qual')[1]
                     cmds.append('''%spigz -dc -p %s %s > %s''' % (paths['pigz_home'], args.n, possible_qual, qnew))
               else:
                  se.append(f)
            else:
               se.append(f)
      
      if args.pe:   
         pe_ftypes = map(mlst_modules.set_filetype, args.pe)
         for i,f in enumerate(args.pe):
            if pe_ftypes[i] == 'fastq':
               cmds.append('%smlst_fastq2fastaqual.py --i %s --p %s' % (paths['mlst_home'], f, os.path.split(f)[1]))
               pe.append(os.path.split(f)[1]+'.fasta')
            elif pe_ftypes[i] == 'fasta':
               if f.endswith('.gz'):
                  fnew = os.path.splitext(os.path.split(f)[1])[0]
                  cmds.append('''%spigz -dc -p %s %s > %s''' % (paths['pigz_home'], args.n, f, fnew))
                  fnew = os.path.splitext(f)[0]
                  pe.append(fnew)
                  
                  # look for qual file (dont add to path because newbler will pick it up)
                  possible_qual = os.path.splitext(os.path.splitext(f)[0])[0] + '.qual.gz'
                  if os.path.exists(possible_qual):
                     qnew = os.path.split(os.path.splitext(os.path.splitext(f)[0])[0] + '.qual')[1]
                     cmds.append('''%spigz -dc -p %s %s > %s''' % (paths['pigz_home'], args.n, possible_qual, qnew))
               else:
                  se.append(f)

            else:
               pe.append(f)
      
      return cmds, se, pe
   
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   cmds = []
   cf = convert_fastq(args, paths)
   cmds.extend(cf[0])
   args.se = cf[1]
   args.pe = cf[2]
   
   cmds.append('newAssembly %s' % args.outpath)
   if args.se:
      for f in args.se:
         cmds.append('addRun -lib shotgun %s %s' % (args.outpath, f))
   if args.pe:
      for i,f in enumerate(args.pe):
         cmds.append('addRun -p -lib PE%i %s %s' % (i, args.outpath, f))
   cmds.append('runProject -cpu %s %s' % (args.n, args.outpath))
   
   # write in bash script
   fh = open('newbler.sh', 'w')
   fh.write('#!/bin/sh\n\n')
   for cmd in cmds:
      fh.write(cmd+'\n')
   fh.close()
   
   # return command (NB. add env. variable to run)
   return ['sh newbler.sh']

def newbler_stats(args):
   '''Create newbler stats calls'''
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   
   cmds = []
   assembly_path = '%s/assembly/' % args.outpath
   cmds.append('''perl -ne 'if ($_ =~ m/^>.+length=(\d+)/) { print $1, "\\n"}' %s > 454AllContigs.lengths ''' % (assembly_path + '454AllContigs.fna'))
   cmds.append('''perl -ne 'if ($_ =~ m/^>.+length=(\d+)/) { print $1, "\\n"}' %s > 454LargeContigs.lengths ''' % (assembly_path + '454LargeContigs.fna'))
   if args.pe:
      cmds.append('''perl -ne 'if ($_ =~ m/^>.+length=(\d+)/) { print $1, "\\n"}' %s > 454Scaffolds.lengths ''' % (assembly_path + '454Scaffolds.fna'))
      cmds.append('%sR --vanilla 454AllContigs.lengths 454LargeContigs.lengths 454Scaffolds.lengths assembly.stats.txt < %smlst_denovo_newbler_stats.R ' % (paths['R_home'], paths['mlst_home']))
   else:
      cmds.append('%sR --vanilla 454AllContigs.lengths 454LargeContigs.lengths NA assembly.stats.txt < %smlst_denovo_newbler_stats.R ' % (paths['R_home'], paths['mlst_home']))
   
   ## write semaphore
   if args.sfile and args.sfile != 'None': cmds.append('echo "done" > %s' % args.sfile)   
   
   # write in bash script
   fh = open('newbler_stats.sh', 'w')
   fh.write('#!/bin/sh\n\n')
   for cmd in cmds:
      fh.write(cmd+'\n')
   fh.close()
   
   # return command (NB. add env. variable to run)
   return ['sh newbler_stats.sh']

def start_assembly(args, logger):
   '''start newbler assembly'''
   
   import mlst_modules
   from mlst_classes import Moab
   from mlst_classes import Semaphore   
   import os
   
   # set queueing
   paths = mlst_modules.setSystem()
   home = os.getcwd()
   if args.partition == 'uv':
      cpuV = 'ncpus=%i,mem=%s,walltime=172800' % (args.n, args.m)
      cpuA = 'ncpus=1,mem=512mb,walltime=172800'
      cpuC = 'ncpus=1,mem=2gb,walltime=172800'
      cpuE = 'ncpus=1,mem=5gb,walltime=172800'
      cpuF = 'ncpus=2,mem=2gb,walltime=172800'
      cpuB = 'ncpus=16,mem=10gb,walltime=172800'      
   else:
      cpuV = 'nodes=1:ppn=%i,mem=%s,walltime=172800' % (args.n, args.m)
      cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
      cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
      cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
      cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
      cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   newbler_calls = newbler(args)
   newblerstats_calls = newbler_stats(args)
   
   # set environment variable (add newbler binaries to bin):
   env_var = 'PATH=/panvol1/simon/bin/454/bin/'
   
   # submit and release jobs
   print "Submitting jobs"
   newbler_moab = Moab(newbler_calls, logfile=logger, runname='run_mlst_newbler', queue=args.q, cpu=cpuV, env=env_var, partition=args.partition, host='cge-s2.cbs.dtu.dk')
   newblerstats_moab = Moab(newblerstats_calls, logfile=logger, runname='run_mlst_newblerstats', queue=args.q, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=newbler_moab.ids, partition=args.partition, host='cge-s2.cbs.dtu.dk')
   
   # release jobs
   newbler_moab.release('cge-s2.cbs.dtu.dk')
   newblerstats_moab.release('cge-s2.cbs.dtu.dk')
      

if __name__ == '__main__':
   
   import os
   import logging
   
   # create the parser
   parser = argparse.ArgumentParser(prog='mlst_denovo_newbler.py', formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog,max_help_position=50, width=110), usage='%(prog)s [options]', description='''
   Run Newbler (2.6) denovo assembly
   Format can be: fasta, fastq, sff
   If both fasta and qual files should be used for assembly only add the fasta file 
   and make sure the qual file is in the same path with same name except .qual ending''')
   
   parser.add_argument('--se', help='input sff/fasta/fastq files', nargs='+', action=set_abspath())
   parser.add_argument('--pe', help='input paired end sff/fasta/fastq files', nargs='+', default=None, action=set_abspath())
   parser.add_argument('--sample', help='name of run and output directory [newbler_assembly]', default='newbler_assembly')
   parser.add_argument('--outpath', help='assembly output dir [denovo]', default='denovo')
   parser.add_argument('--n', help='number of cpus for parallel run [4]', default=4, type=int)
   parser.add_argument('--m', help='memory needed for assembly [2gb]', default='2gb')
   parser.add_argument('--q', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cge]', default='cge')
   parser.add_argument('--partition', help='partition to run on (cge-cluster, uv) [cge-cluster]', default='cge-cluster')
   parser.add_argument('--log', help='log level [info]', default='info')
   parser.add_argument('--sfile', help='semaphore file for waiting [None]', default=None)
   
   args = parser.parse_args()
   #args = parser.parse_args(' --q cge --log info --m 7gb --n 4 --sample None --outpath denovo --se /panvol1/simon/projects/cge/test/GOS1.sff --sfile semaphores/denovoIICMLDK1EI'.split())
   
   # change to sample dir if set
   if args.sample and args.sample != 'None':
      if not os.path.exists(args.sample):
         os.makedirs(args.sample)
      os.chmod(args.sample, 0777)
      os.chdir(args.sample)
   else:
      pass
   
   # start logging
   logger = logging.getLogger('mlst_denovo_newbler.py')
   hdlr = logging.FileHandler('mlst_denovo_newbler.log')
   formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
   hdlr.setFormatter(formatter)
   logger.addHandler(hdlr) 
   if args.log == 'info':
      logger.setLevel(logging.INFO)
   
   # start assembly
   start_assembly(args, logger)
