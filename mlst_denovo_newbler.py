#!/panvol1/simon/bin/python2.7

import argparse

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
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   cmds = []
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
      cmds.append('%sR-2.12 --vanilla 454AllContigs.lengths 454LargeContigs.lengths 454Scaffolds.lengths assembly.stats.txt < %smlst_denovo_newbler_stats.R ' % (paths['R_home'], paths['mlst_home']))
   else:
      cmds.append('%sR-2.12 --vanilla 454AllContigs.lengths 454LargeContigs.lengths NA assembly.stats.txt < %smlst_denovo_newbler_stats.R ' % (paths['R_home'], paths['mlst_home']))
   
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
   newbler_moab = Moab(newbler_calls, logfile=logger, runname='run_mlst_newbler', queue=args.q, cpu=cpuV, env=env_var)
   newblerstats_moab = Moab(newblerstats_calls, logfile=logger, runname='run_mlst_newblerstats', queue=args.q, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=newbler_moab.ids)
   
   # release jobs
   newbler_moab.release()
   newblerstats_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(newblerstats_moab.ids, home, 'newbler', args.q, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   

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
   parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cge]', default='cge')
   parser.add_argument('--log', help='log level [info]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--se 454Reads.fna --pe E_coli_PE.sff --sample Ecoli'.split())
   
   # change to sample dir if set
   if args.sample:
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
