import platform

###   SET SYSTEM PATH   ###

def setSystem():
   '''Identifies architecture and sets system specific
      variables such as executable paths'''
      
   plat = platform.uname()
   
   sys = {}
         
   if plat[4] == 'x86_64':
      sys['bin_home'] = '/panvol1/simon/bin/'
      sys['python2.7_home'] = '/panvol1/simon/bin/'
      sys['java_home'] = '/panvol1/simon/bin/jre1.6.0_21/bin/'
      sys['blastall_home'] = '/panfs/saqqaq/bin/'
      sys['blastdb'] = '/panvol1/simon/databases/blastdb/'
      sys['blat_home'] = '/panfs/saqqaq/bin/'
      sys['fastxtoolkit_home'] = '/panvol1/simon/bin/'
      sys['bowtie_home'] = '/panvol1/simon/bin/bowtie-0.12.5/'
      sys['bwa_home'] = '/panvol1/simon/bin/bwa-0.5.9/'
      sys['samtools_svn_home'] = '/panvol1/simon/bin/samtools_svn/'
      sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.16/'
      sys['picard_home'] = '/panvol1/simon/bin/picard-tools-1.26/'
      sys['bedtools_home'] = '/panvol1/simon/bin/BEDTools-2.8.3/'
      sys['velvet_home'] = '/panvol1/simon/bin/velvet_1.1.02/'
      sys['newbler'] = '/panvol1/simon/bin/454/bin/'
      sys['R_home'] = '/tools/bin/'
      sys['mlst_home'] = '/panvol1/simon/bin/mlst/'
   else:
      raise ValueError('Platform not identified')
   
   return(sys)

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

def set_filetype(f):
   '''Detects filetype from fa, fq and sff input file'''
   
   inhandle = open(f, "r")
   line = inhandle.readline()
   if line.startswith(">"):
      out = 'fasta'
   elif line.startswith("@"):
      out = 'fastq'
   else:
      inhandle = open(f, "rb")
      line = inhandle.readline()
      if line.startswith(".sff"):
         out = 'sff'
      else:
         raise ValueError('Input must be fasta, fastq or sff')
   
   return out

def set_fqtype(f):
   '''Detects sanger or illumina format from fastq'''
   
   # Open fastq, convert ASCII to number, check if number is above or below certain thresholds to determine format
   
   from Bio.SeqIO.QualityIO import FastqGeneralIterator
   
   type = 'not determined'
   inhandle = open(f, 'r')
   for (title, sequence, quality) in FastqGeneralIterator(inhandle):
      qs = map(ord, quality)
      for q in qs:
         if q > 73:
            type = 'Illumina'
            break
         elif q < 59:
            type = 'Sanger'
            break
      if type != 'not determined':
         break
   
   if type == 'not determined':
      raise ValueError('Fastq format not identified, are you sure it is sanger/illumina?')
   else:
      return type

