#!/tools/opt/python/python2.7.2/bin/python2.7

import argparse
import gzip

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

# should be run on trimmed data(!)
def set_kmersizes(args, no_reads=10000, step=4):
   '''Automatically sets kmersizes from avg read lengths'''
   
   from Bio.SeqIO.QualityIO import FastqGeneralIterator
   from Bio import SeqIO
   import sys
   
   def average(values):
      '''Computes the arithmetic mean of a list of numbers.'''
      return sum(values, 0.0) / len(values)
   
   def get_readlengths(f, no_reads):
      '''Get avg. readlengths of no_reads'''
            
      def get_filetype(f):
         '''Return filetype (fasta/fastq)'''
         
         if f.endswith('.gz'):
            inhandle = gzip.open(f, 'rb')
         else:
            inhandle = open(f, "r")
         
         line = inhandle.readline()
         if line.startswith(">"):
            out = 'fasta'
         elif line.startswith("@"):
            out = 'fastq'
         else:
            raise ValueError('Input must be fasta or fastq')
         return out
      
      # main
      if f.endswith('.gz'):
         fh = gzip.open(f, 'rb')
      else:
         fh = open(f, 'r')
      
      count = 0
      L = []
      if get_filetype(f) == 'fastq':
         for (title, sequence, quality) in FastqGeneralIterator(fh):
            if count < no_reads:
               L.append(len(sequence))
               count = count + 1
            else:
               return average(L)
         # if not returned yet (less than no_reads in the file return now)
         return average(L)
      elif get_filetype(f) == 'fasta':
         for rec in SeqIO.parse(fh, 'fasta'):
            if count < no_reads:
               L.append(len(rec.seq))
               count = count + 1
            else:
               return average(L)
         # if not returned yet (less than no_reads in the file return now)
         return average(L)
      else:
         raise ValueError('Input must be fasta or fastq')
   
   def floor_to_odd(number):
      number = int(number)
      if number%2==0:
         return number - 1
      else:
         return number
   
   # main
   files = []
   if args.short: files.extend(args.short)
   if args.short2: files.extend(args.short2)
   if args.shortPaired: files.extend(args.shortPaired)
   if args.shortPaired2: files.extend(args.shortPaired2)
   
   # remove any occurences of fasta, fastq, fastq.gz, fasta.gz
   files = filter (lambda a: a != 'fasta', files)
   files = filter (lambda a: a != 'fastq', files)
   files = filter (lambda a: a != 'fasta.gz', files)
   files = filter (lambda a: a != 'fastq.gz', files)
   
   # get average lengths
   avg_lengths = []
   for f in files:
      avg_f = get_readlengths(f, no_reads)
      avg_lengths.append(avg_f)
   
   # get average of files and set ksizes from this
   avg = average(avg_lengths)
   if avg <= 75: ksizes = [str(floor_to_odd(avg/3*1)), str(floor_to_odd(avg/3*2.5)), '4']
   elif avg <= 100: ksizes = [str(floor_to_odd(avg/3*1.25)), str(floor_to_odd(avg/3*2.25)), '4']
   elif avg <= 150: ksizes = [str(floor_to_odd(avg/3*1)), str(floor_to_odd(avg/3*2.25)), '4']
   
   # change min ksizes if it is smaller than 15. Memory requirements becomes very large at k < 15
   # change max ksizes if larger than 99, only compiled to 99 as max
   if int(ksizes[0]) < 15: ksizes[0] = '15'
   if int(ksizes[1]) > 99: ksizes[1] == '99'
   sys.stderr.write('Ksizes set to (min, max, step) %s\n' % ' '.join(ksizes))
   return ksizes
      

def illumina_trim(args, min_length, min_baseq, min_avgq, min_adaptor_match, keep_n):
   '''Create single end trim calls'''
   
   import os
   import mlst_modules
   paths = mlst_modules.setSystem()
   
   cmd = '%smlst_illumina_trim_h.py' % (paths['mlst_home'])
   calls = []
   if args.short:
      if args.short[0] == 'fastq' or args.short[0] == 'fastq.gz':
         outfiles_short = []   
         for i,f in enumerate(args.short):
            if i == 0: continue
            outfile_short = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
            outfiles_short.append(outfile_short)
            arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s ' % (f, min_length,
               min_baseq, min_avgq, min_adaptor_match, outfile_short)
            if keep_n: arg = arg + ' --keep_n'
            calls.append(cmd+arg)
         args.short[1:] = outfiles_short
   
   if args.short2:
      if args.short2[0] == 'fastq' or args.short2[0] == 'fastq.gz':
         outfiles_short2 = []   
         for i,f in enumerate(args.short2):
            if i == 0: continue
            outfile_short2 = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
            outfiles_short2.append(outfile_short2)
            arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i  --min_adaptor_match %i --o %s ' % (f, min_length,
               min_baseq, min_avgq, min_adaptor_match, outfile_short2)
            if keep_n: arg = arg + ' --keep_n'
            calls.append(cmd+arg)
         args.short2[1:] = outfiles_short2
   
   if args.shortPaired and args.shortPaired[0].find('fastq') > -1:
      outfiles_shortPaired = []
      if len(args.shortPaired) == 3:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired[2])[1] + '.trim.fq'
         outfiles_shortPaired.append(outfile_pe1)
         outfiles_shortPaired.append(outfile_pe2)
         arg = ' --i %s %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired[1], args.shortPaired[2], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      elif len(args.shortPaired) == 2:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '_1.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '_2.trim.fq'
         outfiles_shortPaired.append(outfile_pe1)
         outfiles_shortPaired.append(outfile_pe2)
         arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired[1], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      else:
         raise ValueError('Length of input to shortPaired is not correct')
      
      if keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
      args.shortPaired[1:] = outfiles_shortPaired
   
   if args.shortPaired2 and args.shortPaired[0].find('fastq') > -1:
      outfiles_shortPaired2 = []
      if len(args.shortPaired2) == 3:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired2[2])[1] + '.trim.fq'
         outfiles_shortPaired2.append(outfile_pe1)
         outfiles_shortPaired2.append(outfile_pe2)
         arg = ' --i %s %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired2[1], args.shortPaired2[2], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      elif len(args.shortPaired2) == 2:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '_1.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '_2.trim.fq'
         outfiles_shortPaired2.append(outfile_pe1)
         outfiles_shortPaired2.append(outfile_pe2)
         arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired2[1], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      else:
         raise ValueError('Length of input to shortPaired2 is not correct')
      if keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
      args.shortPaired2[1:] = outfiles_shortPaired2
   
   if len(calls) > 0:
      if not os.path.exists('trimmed'):
         os.makedirs('trimmed')
   return calls
        

def create_velvet_calls(args):
   '''Create velvet calls'''
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   
   # VELVETH CALLS
   # create calls, outpath, ksizes, format, readtypes, reads
   cmd = '%svelveth' % paths['velvet_home']
   velveth_calls = []
   if len(args.ksizes) == 1:
      arg = ' %s %s -create_binary ' % (args.outpath, args.ksizes[0])
      if args.short: arg = arg + ' -short -%s %s' % (args.short[0], ' '.join(args.short[1:]))
      if args.short2: arg = arg + ' -short2 -%s %s' % (args.short2[0], ' '.join(args.short2[1:]))
      if args.shortPaired:
         if len(args.shortPaired) == 2:
            arg = arg + ' -shortPaired -%s %s' % (args.shortPaired[0], args.shortPaired[1])
         elif len(args.shortPaired) == 3:
            arg = arg + ' -shortPaired -separate -%s %s %s' % (args.shortPaired[0], args.shortPaired[1], args.shortPaired[2])
      if args.shortPaired2:
         if len(args.shortPaired2) == 2:
            arg = arg + ' -shortPaired2 -%s %s' % (args.shortPaired2[0], args.shortPaired2[1])
         elif len(args.shortPaired) == 3:
            arg = arg + ' -shortPaired2 -separate -%s %s %s' % (args.shortPaired2[0], args.shortPaired2[1], args.shortPaired2[2])
      if args.long: arg = arg + ' -long -%s %s' % (args.long[0], ' '.join(args.long[1:]))
      if args.longPaired:
         if len(args.longPaired) == 2:
            arg = arg + ' -longPaired -%s %s' % (args.longPaired[0], args.longPaired[1])
         elif len(args.longPaired) == 3:
            arg = arg + ' -longPaired -separate -%s %s %s' % (args.longPaired[0], args.longPaired[1], args.longPaired[2])
      if args.add_velveth: arg = arg + ' %s' % args.add_velveth
      call = cmd + arg
      velveth_calls.append(call)
   
   elif len(args.ksizes) >= 2 and len(args.ksizes) <= 3:
      if len(args.ksizes) == 2:
         step = 2
      elif len(args.ksizes) == 3:
         step = args.ksizes[2]
      
      # create calls, outpath, ksizes, format, readtypes, reads
      for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
         arg = ' %s_%s %s -create_binary ' % (args.outpath, k, k)
         if args.short: arg = arg + ' -short -%s %s' % (args.short[0], ' '.join(args.short[1:]))
         if args.short2: arg = arg + ' -short2 -%s %s' % (args.short2[0], ' '.join(args.short2[1:]))
         if args.shortPaired:
            if len(args.shortPaired) == 2:
               arg = arg + ' -shortPaired -%s %s' % (args.shortPaired[0], args.shortPaired[1])
            elif len(args.shortPaired) == 3:
               arg = arg + ' -shortPaired -separate -%s %s %s' % (args.shortPaired[0], args.shortPaired[1], args.shortPaired[2])
         if args.shortPaired2:
            if len(args.shortPaired2) == 2:
               arg = arg + ' -shortPaired2 -%s %s' % (args.shortPaired2[0], args.shortPaired2[1])
            elif len(args.shortPaired) == 3:
               arg = arg + ' -shortPaired2 -separate -%s %s %s' % (args.shortPaired2[0], args.shortPaired2[1], args.shortPaired2[2])
         if args.long: arg = arg + ' -long -%s %s' % (args.long[0], ' '.join(args.long[1:]))
         if args.longPaired:
            if len(args.longPaired) == 2:
               arg = arg + ' -longPaired -%s %s' % (args.longPaired[0], args.longPaired[1])
            elif len(args.longPaired) == 3:
               arg = arg + ' -longPaired -separate -%s %s %s' % (args.longPaired[0], args.longPaired[1], args.longPaired[2])
         if args.add_velveth: arg = arg + ' %s' % args.add_velveth
         call = cmd + arg
         velveth_calls.append(call)
   else:
      raise ValueError('ksizes must be one value giving ksize, two values giving lower and upper limit (step will be 2) or three values giving lower limit, upper limit and step')  
   
   # VELVETG CALLS
   # create cmd
   cmd = '%svelvetg' % paths['velvet_home']
   cmds = []
   if len(args.ksizes) == 1:
      cmd = '%svelvetg %s' % (paths['velvet_home'], args.outpath)
      cmds.append(cmd)
   elif len(args.ksizes) >= 2 and len(args.ksizes) <= 3:
      if len(args.ksizes) == 2:
         step = 2
      elif len(args.ksizes) == 3:
         step = args.ksizes[2]
      
      for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
         cmd = '%svelvetg %s_%s' % (paths['velvet_home'], args.outpath, k)
         cmds.append(cmd)
   
   # create arg: cov_cutoff, exp_cov, ins_length, add_velvetg
   velvetg_calls = []
   # add other parameters
   for i in range(len(cmds)):
      arg = ' -min_contig_lgth %i' % args.min_contig_lgth
      if args.cov_cutoff: arg = arg + ' -cov_cutoff %f' % args.cov_cutoff
      if args.exp_cov != "None": arg = arg + ' -exp_cov %s' % args.exp_cov
      if args.ins_length: arg = arg + ' -ins_length %i' % args.ins_length
      if args.add_velvetg: arg = arg + ' %s' % args.add_velvetg
      velvetg_calls.append(cmds[i]+arg)
   
   # COMBINE IN SH-FILES #
   sh_calls = []
   for i in range(len(velveth_calls)):
      fh = open('velvet%i.sh' % i, 'w')
      fh.write('#!/bin/sh\n\n')
      fh.write(velveth_calls[i]+'\n')
      fh.write(velvetg_calls[i]+'\n')
      fh.close()
      sh_calls.append('sh velvet%i.sh' %i)
   return sh_calls
   

def postprocess(args):
   '''Determine best assembly, remove other assemblies, clean up and write semaphore file (if required)'''
   
   import mlst_modules
   paths = mlst_modules.setSystem()
   
   calls = []
   if len(args.ksizes) > 1:
      ## parse_assemblies
      cmd = '%sR --vanilla ' % paths['R_home']
      
      # set argument
      if len(args.ksizes) == 1:
         arg = ' %s %s' % (args.outpath, args.ksizes[0])
      elif len(args.ksizes) >= 2:
         if len(args.ksizes) == 2:
            step = 2
         elif len(args.ksizes) == 3:
            step = args.ksizes[2]
         
         arg_list = []
         for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
            out = '%s_%s/stats.txt %s' % (args.outpath, k, k)
            arg_list.append(out)
         arg = ' '.join(arg_list)
      
      call = cmd + arg + ' < %smlst_denovo_velvet_parse.R' % (paths['mlst_home'])
      calls.append(call)
      
      ## accept assembly
      call = '%smlst_denovo_velvet_accept.py %s' % (paths['mlst_home'], args.outpath)
      calls.append(call)
      
   ## clean
   call = '%smlst_denovo_velvet_clean.py' % (paths['mlst_home'])
   calls.append(call)
   
   ## write semaphore
   if args.sfile and args.sfile != 'None': calls.append('echo "done" > %s' % args.sfile)
      
   ## write in bash script
   fh = open('postprocess.sh', 'w')
   fh.write('#!/bin/sh\n\n')
   for call in calls:
      fh.write(call+'\n')
   fh.close()
   
   return ['sh postprocess.sh']
   

def start_assembly(args, logger):
   '''Start assembly'''
   
   import mlst_modules
   from mlst_classes import Moab
   from mlst_classes import Semaphore   
   import os
   
   # set queueing
   paths = mlst_modules.setSystem()
   home = os.getcwd()
   if args.partition == 'uv':
      cpuV = 'procs=%i,mem=%s,walltime=172800,flags=sharedmem' % (args.n, args.m)
      cpuA = 'procs=1,mem=512mb,walltime=172800,flags=sharedmem'
      cpuC = 'procs=1,mem=2gb,walltime=172800,flags=sharedmem'
      cpuE = 'procs=1,mem=5gb,walltime=172800,flags=sharedmem'
      cpuF = 'procs=2,mem=%s,walltime=172800,flags=sharedmem' % args.m
      cpuB = 'procs=16,mem=10gb,walltime=172800,flags=sharedmem'      
   else:
      cpuV = 'nodes=1:ppn=%i,mem=%s,walltime=172800' % (args.n, args.m)
      cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
      cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
      cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
      cpuF = 'nodes=1:ppn=2,mem=%s,walltime=172800' % args.m
      cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
      
   # set kmersizes (if auto)
   if args.ksizes == ['auto']:
      args.ksizes = set_kmersizes(args)
   
   # trimming calls
   if args.trim:
      illuminatrim_calls = illumina_trim(args, int(args.ksizes[0]), 15, 20, 15, False)
      if not os.path.exists('trimmed'):
         os.makedirs('trimmed')
      
   # velvet calls
   velvet_calls = create_velvet_calls(args)
      
   # velvet parse calls
   postprocess_calls = postprocess(args)
   
   # set environment variable:
   env_var = 'OMP_NUM_THREADS=%i' % int(args.n - 1)
   
   # submit and release jobs
   # NB: mlst_denovo_velvet is run from a compute node, it will then ssh to "host" and submit the jobs from there (cge-s2)
   print "Submitting jobs"
   
   # if trimming is needed
   if args.trim:
      illuminatrim_moab = Moab(illuminatrim_calls, logfile=logger, runname='run_mlst_trim', queue=args.q, cpu=cpuF, partition=args.partition, host='cge-s2.cbs.dtu.dk')
      velvet_moab = Moab(velvet_calls, logfile=logger, runname='run_mlst_velvet', queue=args.q, cpu=cpuV, depend=True, depend_type='all', depend_val=[1], depend_ids=illuminatrim_moab.ids, env=env_var, partition=args.partition, host='cge-s2.cbs.dtu.dk')
   # if no trimming
   else:
      velvet_moab = Moab(velvet_calls, logfile=logger, runname='run_mlst_velvet', queue=args.q, cpu=cpuV, env=env_var, partition=args.partition, host='cge-s2.cbs.dtu.dk')
   
   # submit job for postprocessing
   postprocess_moab = Moab(postprocess_calls, logfile=logger, runname='run_mlst_postprocess', queue=args.q, cpu=cpuA, depend=True, depend_type='conc', depend_val=[len(velvet_calls)], depend_ids=velvet_moab.ids, partition=args.partition, host='cge-s2.cbs.dtu.dk')
   
   # release jobs
   print "Releasing jobs"
   if args.trim and len(illuminatrim_calls) > 0: illuminatrim_moab.release(host='cge-s2.cbs.dtu.dk')
   velvet_moab.release('cge-s2.cbs.dtu.dk')
   postprocess_moab.release(host='cge-s2.cbs.dtu.dk')
   

if __name__ == '__main__':
   import os
   import logging
   
   # create the parser
   parser = argparse.ArgumentParser(prog='mlst_denovo_velvet.py', formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog,max_help_position=50, width=110), usage='%(prog)s [options]', description='''
   Run Velvet denovo assembly. Read input is given by eg. --short <format> <reads>
   Format can be: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam
   add_velvetg/add_velveth has to be added with quotes, eg: add_velvetg "-very_clean yes"
   if mate_pairs are used (in eg. shortPaired2) one should add: add_velvetg "-shortMatePaired2 yes -ins_length2 <MP size>"
   NB: Trimming does not work on already interleaved files\n''')
   
   # add the arguments
   parser.add_argument('--short', help='input read format and short reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--shortPaired', help='input read format and short paired reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--short2', help='input read format and short2 reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--shortPaired2', help='input read format and short paired2 reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--long', help='input read format and long reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--longPaired', help='input read format and long paired reads', nargs='+', action=required_nargs(0,3))
   parser.add_argument('--ksizes', help='kmers to run assemblies for (single (m) or m M s (min, max, step)) [auto]', nargs='+', default=['auto'])
   parser.add_argument('--sample', help='name of run and output directory [velvet_assembly]', default='velvet_assembly')
   parser.add_argument('--outpath', help='assembly output dir [assembly]', default='assembly')
   parser.add_argument('--trim', help='should input files be trimmed (illumina only) [False]', default=False, action='store_true')
   parser.add_argument('--min_contig_lgth', help='mininum length to report contig [100]', default=100, type=int)
   parser.add_argument('--cov_cutoff', help='coverage cutoff for removal of low coverage (float) [None]', default=None, type=float)
   parser.add_argument('--exp_cov', help='Expected mean coverage (None, float, auto) [auto]', default='auto')
   parser.add_argument('--ins_length', help='insert size (reads included) [None]', default=None, type=int)
   parser.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
   parser.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)
   parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
   parser.add_argument('--m', help='memory needed for assembly [2gb]', default='2gb')
   parser.add_argument('--partition', help='partition to run on (cge-cluster, uv) [cge-cluster]', default='cge-cluster')
   parser.add_argument('--q', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cge]', default='cge')
   parser.add_argument('--log', help='log level [info]', default='info')
   parser.add_argument('--sfile', help='semaphore file for waiting [None]', default=None)

   args = parser.parse_args()
   #args = parser.parse_args('--short fastq test_1.fq test_2.fq --ksizes 33 49 4 --outpath test'.split())
   #args = parser.parse_args('--short fastq.gz Kleb-10-213361_2.interleaved.fastq.test.gz --ksizes 41 55 4 --sample Kleb'.split())
   #args = parser.parse_args('--shortPaired fastq /panfs1/cge/data/cge_private/s.aureus/lane7_sample28_TG8130/s_7_1_sequence.txt /panfs1/cge/data/cge_private/s.aureus/lane7_sample28_TG8130/s_7_2_sequence.txt --ksizes 21 45 4 --sample TG8130'.split())
   #args = parser.parse_args('--shortPaired Kleb-10-213361_2_1_sequence.txt Kleb-10-213361_2_2_sequence.txt --ksizes 33 75 4 --sample kleb_wtrim --trim'.split())
   #args = parser.parse_args('--short /panvol1/simon/projects/cge/EHEC/illumina/bgi/110601_I238_FCB067HABXX_L3_ESCqslRAADIAAPEI-2_1.fq --outpath assembly'.split())
   #args = parser.parse_args('--sample None --outpath assembly --min_contig_lgth 100 --exp_cov auto --log info --sfile None --short /net/cxfs1-10g/home/projects3/squid/rapsearch/reads/DFPA.reads.fq --ksizes auto --partition cge-cluster --m 1gb --n 1 --q cge'.split())
   #args = parser.parse_args('--shortPaired fastq.gz test_1.fastq.gz test_2.fastq.gz --sample testgz'.split())
   
   
   # add_velveth and add_velvetg works from commandline, eg:
   # mlst_denovo_velvet.py --short fastq.gz interleaved.fastq.gz --ksizes 33 --sample Kleb --add_velvetg "-very_clean yes"
   
   # change to sample dir if set
   if args.sample and args.sample != 'None':
      if not os.path.exists(args.sample):
         os.makedirs(args.sample)
      os.chmod(args.sample, 0777)
      os.chdir(args.sample)
   else:
      pass
   
   if not os.path.exists('log'):
      os.makedirs('log')
   
   # start logging
   logger = logging.getLogger('mlst_denovo_velvet.py')
   hdlr = logging.FileHandler('mlst_denovo_velvet.log')
   formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
   hdlr.setFormatter(formatter)
   logger.addHandler(hdlr) 
   if args.log == 'info':
      logger.setLevel(logging.INFO)
   
   # start assembly
   start_assembly(args, logger)
