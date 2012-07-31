#!/panvol1/simon/bin/python2.7

import argparse
import gzip
from Bio import SeqIO

def fastq2fasta(f, p):
   '''Convert fastq to fasta'''
   
   # open handle
   if f.endswith('.gz'):
      fh = gzip.open(f, 'rb')
   else:
      fh = open(f, 'r')
   
   # convert
   SeqIO.convert(fh, "fastq", p+'.fasta', "fasta")
   fh.close()

def fastq2qual(f, p):
   '''Convert fastq to qual'''
   
   # open handle
   if f.endswith('.gz'):
      fh = gzip.open(f, 'rb')
   else:
      fh = open(f, 'r')
   
   SeqIO.convert(fh, "fastq", p+'.qual', "qual")
   fh.close()


if __name__ == '__main__':
   parser = argparse.ArgumentParser(description=
      '''Convert fastq to fasta + qual files''')
   
   parser.add_argument('--i', help='input fastq file', required=True)
   parser.add_argument('--p', help='output prefix', required=True)
   
   args = parser.parse_args()
   #args = parser.parse_args('--i test_1.fastq.gz --p tmp'.split())
   
   fastq2fasta(args.i, args.p)
   fastq2qual(args.i, args.p)