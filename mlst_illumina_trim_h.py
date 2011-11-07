#!/panvol1/simon/bin/python2.7

from __future__ import division
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import mlst_modules
import sys
import random
import string
import cdb
import multiprocessing

class FastqTrim:
   '''Trim/Filter paired end fastq:
      trim 3' trailing bases under (q)
      trim sequencing adaptors (Illumina) (a)
      filter reads with less than avg. Q quals (m)
      filter reads with less than L length (l)
      filter reads with Ns (keep_n)
   '''
   
   def __init__(self, f, o, l=25, q=20, m=20, keep_n=False, M=20, 
                a=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 
                'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 
                'ACACTCTTTCCCTACACGACGCTCTTCCGATCT']):
      self.f = f
      self.o = o
      self.l = l
      self.q = q
      self.m = m
      self.keep_n = keep_n
      self.M = M
      self.a = a
      self.paired = False
      self.interleaved = False
      
      # set adaptors, readtype, fastq format
      self.adaptors = self.set_adaptors()
      self.set_readtype()
      self.format = mlst_modules.set_filetype(self.f[0])
      if self.format == 'fastq':
         self.fqtype = mlst_modules.set_fqtype(self.f[0])
         sys.stderr.write('Header/Quality encoding: %s/%s\n' % (self.fqtype[0], self.fqtype[1]))
   
   def __repr__(self):
      msg = 'FastqTrim(%s, %s, %i, %i, %i, %s, %i, %s, %s)' % ("["+", ".join(self.f)+"]", "["+", ".join(self.o)+"]", self.l, self.q, self.m, self.keep_n, self.M, self.paired, self.interleaved)
      return msg
   
   def set_readtype(self):
      '''Set paired and/or interleaved flags'''
            
      def check_interleaved(self):
         '''Check if a file is interleaved paired or single'''
         
         fh = open(self.f[0], 'r')
         count = 0
         ends = []
         if self.fqtype[0] == 'Illumina1.4':
            for (title, sequence, quality) in FastqGeneralIterator(fh):
               ends.append(title[-2:])
               count += 1
               if count > 9: break
            
            # check that they are interleaved (eg. /1, /2, /1, /2)
            interleaved = True
            for i,e in enumerate(ends[::2]):
               if e != '/1':
                  interleaved = False
                  break
            
            for i,e in enumerate(ends[1::2]):
               if e != '/2':
                  interleaved = False
                  break
         elif self.fqtype[0] == 'Illumina1.8':
            for (title, sequence, quality) in FastqGeneralIterator(fh):
               ends.append(title)
               count += 1
               if count > 9: break
            
            interleaved = True
            for i,e in enumerate(ends[::2]):
               fields = e.split(' ')
               if fields[1].startswith('2'):
                  interleaved = False
                  break
            for i,e in enumerate(ends[1::2]):
               fields = e.split(' ')
               if fields[1].startswith('1'):
                  interleaved = False
                  break
         
         if interleaved:
            self.paired = True
            self.interleaved = True
         
         return self
      
      
      if len(self.f) == 2:
         self.paired = True
         self.interleaved = False
      elif len(self.f) == 1:
         self = check_interleaved(self)
   
   def set_adaptors(self):
      '''Set adaptors according to adaptor length requirements'''
      
      adaptors = []
      for a in self.a:
         if self.M == 0:
            adaptor = a
         elif (self.M) > len(a):
            adaptor = a
         else:
            adaptor = a[:self.M]
         adaptors.append(adaptor)
      return adaptors
   
   def filter_adaptor(self, title, sequence, quality):
      '''Filters adaptor and Ns in read'''
      
      has_adaptor = False
      has_N = False
      is_short = False
      
      # search for adaptor in reads
      for adaptor in self.adaptors:
         hit = sequence.find(adaptor)
         if hit > -1:
            sequence = sequence[:hit]
            quality = quality[:hit]
            has_adaptor = True
      
      # check if read contains Ns
      if sequence.find('N') > -1:
         has_N = True
      
      # check if too small then do not print out
      if len(sequence) < self.l:
         is_short = True
      
      # return read
      if has_N or is_short:
         return (None, None, None)
      else:
         return (title, sequence, quality)
      
   # the actual read trimming function
   def trim_qual(self, title, sequence, quality):
      '''Trims read based on 3' qualities and average quality'''
      
      # functions for read trimming
      def average(values):
         '''Computes the arithmetic mean of a list of numbers.'''
         return sum(values, 0.0) / len(values)
      
      def illumina2qual(qual_string):
         '''Convert illumina qualities (offset 64) to values.'''
         qual = []
         for q in qual_string:
            qual.append(ord(q)-64)
         return qual
      
      def sanger2qual(qual_string):
         '''Convert sanger qualities (offset 33) to values.'''
         qual = []
         for q in qual_string:
            qual.append(ord(q)-33)
         return qual
      
      def trim_index(qual, m):
         '''Determine position to trim from'''
         qs = qual[::-1]
         index = None
         for i,q in enumerate(qs):
            if q > m:
               index = len(qual)-i
               break
         return index
      
      # check if read is not valid
      if title == None:
         return (None, None, None)
      
      if self.fqtype[1] == 'Illumina':
         qual = illumina2qual(quality)
      elif self.fqtype[1] == 'Sanger':
         qual = sanger2qual(quality)
      else:
         raise ValueError('Fastq must be Illumina or Sanger')
      
      # check trailing qualities, should it be trimmed?
      t = trim_index(qual, self.q)
      if t:
         sequence = sequence[:t]
         quality = quality[:t]
         # do not print out if average quality is less than m or length is less than l
         a = average(qual[:t])
         if a < self.m or len(sequence) < self.l:
            return (None, None, None)
         else:
            return (title, sequence, quality)
      else:
         # do not print out if average quality is less than m or length is less than l
         a = average(qual)
         if a < self.m or len(sequence) < self.l:
            return (None, None, None)
         else:
            return (title, sequence, quality)
   
   def write_pairs(self, f1, f2):
      '''Parse through two paired files and only write if both pairs are present'''
      
      def intersect(a, b):
         '''Intesection between lists'''
         return list(set(a) & set(b))
      
      def rm_files(patterns):
         '''Remove files using glob given as list of patterns'''
         
         import glob
         import os
         
         for p in patterns:
            files = glob.glob(p)
            if len(files) == 0:
               pass
            else:
               map(os.remove, files)
      
      def write_out(db_common, f, o):
         '''Write out reads'''
         
         fh = open(f, 'r')
         out = open(o, 'w')
         written_count = 0
         total_count = 0
         if self.fqtype[0] == 'Illumina1.4':
            for (title, sequence, quality) in FastqGeneralIterator(fh):
               total_count += 1
               if db_common.has_key(title[:-2]):
                  out.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
                  written_count += 1
         elif self.fqtype[0] == 'Illumina1.8':
            for (title, sequence, quality) in FastqGeneralIterator(fh):
               total_count += 1
               if db_common.has_key(title.split(' ')[0]):
                  out.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
                  written_count += 1
         
         sys.stderr.write('%s: Total %i, Written %i (%.1f%%)\n' % (f, total_count, written_count, written_count/total_count*100))
         fh.close()
         out.close()
      
      def create_db(f, db_fname):
         '''Write out db of headers'''
         
         fh = open(f, 'r')
         if self.fqtype[0] == 'Illumina1.4':
            fh_headers = (x.strip()[1:-2] for i, x in enumerate(fh) if not (i % 4))
         elif self.fqtype[0] == 'Illumina1.8':
            fh_headers = (x.split(' ')[0][1:] for i, x in enumerate(fh) if not (i % 4))
         else:
            sys.stderr.write('Header encoding not determined: %s\n' % self.fqtype[0])
         
         db = cdb.cdbmake(db_fname, db_fname + '.tmp')
         for h in fh_headers:
            db.add(h, 'T')
         db.finish()
         del(db)
      
      ## get headers from both trimmed files ##
      # strip the /2 or /1 and grab only the headers
      # write in dbm to minimze memory usage
      
      # create db's (parallel)
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(36))
      db1_fname = 'db1_%s' % rand
      db2_fname = 'db2_%s' % rand
      
      jobs = []
      p = multiprocessing.Process(target=create_db, args=(f1, db1_fname, ))
      p.start()
      jobs.append(p)
      
      p = multiprocessing.Process(target=create_db, args=(f2, db2_fname, ))
      p.start()
      jobs.append(p)
      
      # wait for jobs to finish
      for job in jobs:
         job.join()
      
      ## get headers that are in both trimmed files ##
      db1 = cdb.init(db1_fname)
      db2 = cdb.init(db2_fname)
      common = intersect(db1.keys(), db2.keys())
      
      dbcommon_fname = 'dbcommon_%s' % rand
      db_common = cdb.cdbmake(dbcommon_fname, dbcommon_fname + '.tmp')
      for h in common:
         db_common.add(h, 'T')
      db_common.finish()
      del(db_common)
      
      # open common db
      db_common = cdb.init(dbcommon_fname)
      jobs = []
      p = multiprocessing.Process(target=write_out, args=(db_common, f1, self.o[0]))
      p.start()
      jobs.append(p)
      
      p = multiprocessing.Process(target=write_out, args=(db_common, f2, self.o[1]))
      p.start()
      jobs.append(p)
      
      # wait for jobs to finish
      for job in jobs:
         job.join()
      
      rm_files([db1_fname, db2_fname, dbcommon_fname, f1, f2])
   
   def trim(self, paired=False, interleave=False):
      '''Start trimming of reads'''
      
      def trim_file(self, f, f_out):
         written = 0
         total = 0      
         fh_in = open(f, 'r')
         fh_out = open(f_out, 'w')
         for (title, sequence, quality) in FastqGeneralIterator(fh_in):
            total += 1
            (title, sequence, quality) = self.filter_adaptor(title, sequence, quality)
            (title, sequence, quality) = self.trim_qual(title, sequence, quality)
            if title != None:
               if len(sequence) != len(quality):
                  raise ValueError('sequence and quality not of the same length\n%s\n%s\n' % (sequence, quality))
               fh_out.write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
               written += 1
         fh_out.close()
         fh_in.close()
         sys.stderr.write('%s: written %s of %s (%.1f%%)\n' % (f, written, total, written/total))
      
      def trim_interleaved_file(self, f, f_out):
         written = 0
         total = 0      
         fh_in = open(f, 'r')
         fh_out = [open(f_out[0], 'w'), open(f_out[1], 'w')]
         for (title, sequence, quality) in FastqGeneralIterator(fh_in):
            (title, sequence, quality) = self.filter_adaptor(title, sequence, quality)
            (title, sequence, quality) = self.trim_qual(title, sequence, quality)
            if title != None:
               if len(sequence) != len(quality):
                  raise ValueError('sequence and quality not of the same length\n%s\n%s\n' % (sequence, quality))
               fh_out[total%2].write('@%s\n%s\n+\n%s\n' % (title, sequence, quality))
               written += 1
            total += 1
         fh_out[0].close()
         fh_out[1].close()
         fh_in.close()
         sys.stderr.write('%s: written %s of %s (%.1f%%)\n' % (f, written, total, written/total))
            
      if self.paired:
         # set tmp files
         rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(36))
         f_out = ['tmp0.'+rand, 'tmp1.'+rand]
         
         # check if file is interleaved
         
         if self.interleaved:
            sys.stderr.write('Trimming paired end interleaved\n')
            trim_interleaved_file(self, self.f[0], f_out)
         else:
            # create and start multiprocess
            sys.stderr.write('Trimming paired end\n')
            jobs = []
            for i in range(2):
               p = multiprocessing.Process(target=trim_file, args=(self, self.f[i], f_out[i], ))
               p.start()
               jobs.append(p)
            
            # wait for jobs to finish
            for job in jobs:
               job.join()
            
         # write files
         self.write_pairs(f_out[0], f_out[1])
      else:
         sys.stderr.write('Trimming single end\n')
         trim_file(self, self.f[0], self.o[0])

     
if __name__ == '__main__':
      
   parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=110), description='''
      Trim SE or PE files for low q bases, adaptor sequences, Ns
      ''')
   
   # add the arguments
   parser.add_argument('--i', help='input single or paired end files', nargs='+', required=True)
   parser.add_argument('--min_length', help='minimum length of a read to keep pairs [25]', type=int, default=25)
   parser.add_argument('--min_baseq', help='chomp bases with quality less than [20]', default=20, type=int)
   parser.add_argument('--min_avgq', help='minimum average quality of read [20]', default=20, type=int)
   parser.add_argument('--adaptors', help='adaptor sequence to clip [Illumina adaptors]', nargs='+', default=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'])
   parser.add_argument('--keep_n', help='do not remove sequences containing N', default=False, action='store_true')
   parser.add_argument('--min_adaptor_match', help='minimum length of match to adaptor (0=all of adaptor) [20]', default=20, type=int)
   parser.add_argument('--o', help='output files', nargs='+', required=True)
   parser.add_argument('--log', help='log level [INFO]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args(''.split())
   #args = parser.parse_args('--i kleb_test_2.fq --l 25 --q 20 --o kleb_test_2.trim.fq'.split())
   #args = parser.parse_args('--i Kleb-10-213361_2_1_sequence.txt Kleb-10-213361_2_2_sequence.txt --M 15 --o Kleb-10-213361_2_1_sequence.trim.fq Kleb-10-213361_2_2_sequence.trim.fq '.split())
   
   # create instance
   fqtrim = FastqTrim(args.i, args.o, args.min_length, args.min_baseq, args.min_avgq, args.keep_n, args.min_adaptor_match, args.adaptors)
   
   # start trimming
   if len(fqtrim.f) == 2: 
      fqtrim.trim(paired=True)
   elif len(fqtrim.f) == 1: 
      fqtrim.trim()
   else:
      raise ValueError('Only 1 (single end) or 2 (paired end) files must be given')
