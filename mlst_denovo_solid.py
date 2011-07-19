#!/panvol1/simon/bin/python

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

def solid(args):
   '''Create solid assembly calls'''
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   
   cmd = '%sassemble.pl' % paths['solid_home']
   
   if args.se:
      arg = ' %s %i -numcores %i' % (' '.join(args.se), args.rf, args.n)
   elif args.pe:
      arg = ' %s %s %i -f5 %s -f5qv %s -ins_length %i -ins_length_sd %i -numcores %i ' % (args.pe[0], args.pe[1], args.rf, args.pe[2], args.pe[3], args.ins_length, args.ins_length_sd, args.n)
   elif args.mp:
      arg = ' %s %s %i -r3 %s -r3qv %s -ins_length %i -ins_length_sd %i -numcores %i ' % (args.mp[0], args.mp[1], args.rf, args.mp[2], args.mp[3], args.ins_length, args.ins_length_sd, args.n)
   else:
      raise ValueError('Input must be given by --se, --pe or --mp')
   
   # add extra commands
   if args.add_solid: arg = arg + ' ' + args.add_solid
   
   return [cmd+arg]
   
   
def start_assembly(args, logger):
   '''Start assembly of solid reads'''
   
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
   
   solid_calls = solid(args)
   
   # set environment variable (add newbler binaries to bin):
   env_var = 'denovo2=%s' % paths['solid_home']
   
   # submit and release jobs
   print "Submitting jobs"
   solid_moab = Moab(solid_calls, logfile=logger, runname='run_mlst_solid', queue=args.queue, cpu=cpuV, env=env_var)
   
   # release jobs
   solid_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(solid_moab.ids, home, 'solid', args.queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   

if __name__ == '__main__':
   
   import os
   import logging
   
   # create the parser
   parser = argparse.ArgumentParser(prog='mlst_denovo_solid.py', formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog,max_help_position=50, width=110), usage='%(prog)s [options]', description='''
   Run solid denovo assembly
   Input is given as --se, --pe or --mp and only should be used and the order of files must be the ones listed in help.
   If corresponding qual files are not uploaded they should be given by "none" in the input.
   F3=csfasta, F3q=qual, F5=pair.csfasta, F5q=pair.qual, R3=mate.csfasta, R3q=mate.qual
   ''')
   
   parser.add_argument('--se', help='input single end csfasta file (F3, F3q)', nargs='+', action=set_abspath())
   parser.add_argument('--pe', help='input paired end csfasta files (F3, F3q, F5, F5q) ', nargs='+', default=None, action=set_abspath())
   parser.add_argument('--mp', help='input mate pair csfasta files (F3, F3q, R3, R3q)', nargs='+', default=None, action=set_abspath())
   parser.add_argument('--rf', help='input expected length of genome in bp', type=int, required=True)
   parser.add_argument('--ins_length', help='estimate of mate/paired end insert length (1200/170)', type=int)
   parser.add_argument('--ins_length_sd', help='estimate of mate/paired end insert length (300/30)', type=int)
   parser.add_argument('--add_solid', help='additional parameters to solid assembler', default=None)
   parser.add_argument('--sample', help='name of run and output directory [solid_assembly]', default='solid_assembly')
   parser.add_argument('--n', help='number of cpus for parallel run [4]', default=4, type=int)
   parser.add_argument('--m', help='memory needed for assembly [2gb]', default='2gb')
   parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cge]', default='cge')
   parser.add_argument('--log', help='log level [info]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--se  --sample solid_test'.split())
   
   # change to sample dir if set
   if args.sample:
      if not os.path.exists(args.sample):
         os.makedirs(args.sample)
      os.chmod(args.sample, 0777)
      os.chdir(args.sample)
   else:
      pass
   
   # start logging
   logger = logging.getLogger('mlst_denovo_solid.py')
   hdlr = logging.FileHandler('mlst_denovo_solid.log')
   formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
   hdlr.setFormatter(formatter)
   logger.addHandler(hdlr) 
   if args.log == 'info':
      logger.setLevel(logging.INFO)
   
   # start assembly
   start_assembly(args, logger)
   