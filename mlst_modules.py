import platform
import os
import re
import subprocess
import time

###   SET SYSTEM PATH   ###

def setSystem():
   '''Identifies architecture and sets system specific
      variables such as executable paths'''
      
   plat = platform.uname()
   
   sys = {}
         
   if plat[4] == 'x86_64':
      sys['host'] = 'protein-s0'
      sys['bin_home'] = '/panvol1/simon/bin/'
      sys['python2.7_home'] = '/panvol1/simon/bin/'
      sys['java_home'] = '/panvol1/simon/bin/jre1.6.0_21/bin/'
      sys['phymm_home'] = '/panvol1/simon/bin/phymm/'
      sys['pyscripts_home'] = '/panvol1/simon/bin/pipeline/'
      sys['genobox_home'] = '/panvol1/simon/bin/genobox/'
      sys['blastall_home'] = '/panfs/saqqaq/bin/'
      sys['blastdb'] = '/panvol1/simon/databases/blastdb/'
      sys['blat_home'] = '/panfs/saqqaq/bin/'
      sys['fastxtoolkit_home'] = '/panvol1/simon/bin/'
      sys['bowtie_home'] = '/panvol1/simon/bin/bowtie-0.12.5/'
      sys['bwtindexes_home'] = '/panvol1/simon/bin/bowtie-0.12.5/indexes/'
      sys['hs_ref_GRCh37'] = '/panvol1/simon/databases/hs_ref37/'
      sys['hs_ref_ncbi36'] = '/panvol1/simon/databases/hs_ref36/'
      sys['bwa_home'] = '/panvol1/simon/bin/bwa-0.5.9/'
      sys['bwaindexes_home'] = '/panvol1/simon/bin/bwa-0.5.8a/indexes/'
      #sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.12.a/'
      sys['samtools_svn_home'] = '/panvol1/simon/bin/samtools_svn/'
      sys['samtools_home'] = '/panvol1/simon/bin/samtools-0.1.16/'
      sys['picard_home'] = '/panvol1/simon/bin/picard-tools-1.26/'
      sys['bedtools_home'] = '/panvol1/simon/bin/BEDTools-2.8.3/'
      sys['greengenes_home'] = '/panvol1/simon/databases/greengenes/'
      sys['mothur_home'] = '/panvol1/simon/bin/Mothur-1.11.0/'
      sys['rna_hmm3_home'] = '/panvol1/simon/bin/rna_hmm3/'
      sys['hmm-3.0_home'] = '/panvol1/simon/bin/hmmer-3.0/'
      sys['phymm_home'] = '/panvol1/simon/bin/phymm/'
      sys['prodigal_home'] = '/panvol1/simon/bin/prodigal/'
      sys['rnammer_home'] = '/panvol1/simon/bin/rnammer-1.2/'
      sys['cge_genomes'] = '/panvol1/simon/databases/cge_genomes/'
      sys['velvet_home'] = '/panvol1/simon/bin/velvet_1.1.02/'
      #sys['newbler'] = '/panfs/saqqaq/bin/newbler/bin/'
      sys['newbler'] = '/panvol1/simon/bin/454/bin/'
      sys['stampy'] = '/panvol1/simon/bin/stampy-1.0.6/'
      sys['chain_index'] = '/panvol1/simon/databases/chainmap/'
      sys['taxonomy_ncbi'] = '/panvol1/simon/databases/taxonomy/'
      sys['R_home'] = '/tools/bin/'
   else:
      raise ValueError('Platform not identified')
   
   return(sys)
