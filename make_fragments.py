#!/usr/bin/env python
from argparse import ArgumentParser
import os
import glob
from multiprocessing import Pool
parser  = ArgumentParser()
import sys
from shutil import copy

fragment_tools = '/work/robetta/src/rosetta/fragment_tools/'
psipred = '/work/robetta/src/psipred3.21/'

if not 'SPARKSDIR' in os.environ:
  os.environ[ 'SPARKSXDIR' ] = '/work/robetta/src/rosetta_server/bin/sparks-x'

def chckmkdir(path,out = sys.stdout):
  if not os.path.exists(path):
    os.mkdir(path)
    return True
  else:
    return False

def generate_fragments(pdb, args):
  currentdir = os.getcwd()
  basename, extension = os.path.splitext(pdb)
  pdbname = os.path.split(basename)[-1]
  exec_dir = currentdir +"/" +  pdbname + "_fragments"
  exec_dirs = [] 
  dir_maked  = chckmkdir(exec_dir)
  if dir_maked:
    copy(pdb, exec_dir + "/00001.pdb")
    os.chdir(exec_dir)
    os.system("/work/wangyr/scripts/pdb_util/pdb2fasta.py  00001.pdb >  00001.fasta")
    if args.single_psipred:
      if args.csbuild:
        os.system(psipred + "/runpsipred_csbuild_single 00001.fasta")
        os.system("/work/robetta/src/csbuild/csbuild -i  00001.fasta -I fas -D /work/robetta/src/csbuild/csblast-2.2.3/data/K4000.crf -o  00001.check -O chk")
        os.system(fragment_tools + "/make_fragments_tjTest.pl  -nosam -id 00001 00001.fasta -psipredfile 00001.csb.ss2 -nopsiblast_profile")
      else:
        os.system(psipred + "/runpsipred_single 00001.fasta")
        os.system(fragment_tools + "/make_fragments_tjTest.pl  -nosam -id 00001 00001.fasta -psipredfile 00001.ss2 -nopsiblast_profile")
    else:
      os.system(fragment_tools + "/make_fragments.pl  -nosam -nohoms -id 00001 00001.fasta")
    if(args.eval_method == 'nobu'):
      os.system("/work/nobuyasu/mini2/bin/r_frag_quality.default.linuxgccrelease -database /work/nobuyasu/minirosetta_database -in:file:native 00001.pdb -f  00001.200.9mers") 
      os.system("/work/nobuyasu/scripts/count_goodfrag.pl --s=frag_qual.dat > good_fragments_count.txt")
    if(args.eval_method == 'tj'):
      os.system("/work/brunette/src/Rosetta/main/source/bin/fragment_rmsd.linuxgccrelease -database /work/brunette/src/Rosetta/main/database -in:file:s 00001.pdb -in:file:frag3 00001.200.9mers > rmsd_9mers.txt")
  os.chdir(currentdir)

def generate_fragments_single_arg(pdb_args):
    generate_fragments(pdb_args[0],pdb_args[1])
#Note: the multiprocessor support for pool takes only a single argument in python 2.7. If we were running python 3.3 I believe starmap could be used making this function unncessary


if __name__ == '__main__':
  parser.add_argument('-n_proc', default=10, type=int, help='number of cores to use, by default 10')
  parser.add_argument('-single_psipred', default = True,  help='use psipred single sequence prediction')
  parser.add_argument('-no-single_psipred',dest='single_psipred',action='store_false')
  parser.add_argument('-csbuild',default= True ,help='use csbuild to generate the checkpoint file only valid for single sequence')
  parser.add_argument('-no-csbuild',dest='csbuild',action='store_false')
  parser.add_argument('-pdbs', nargs='+', help='pdbs to proccess')
  parser.add_argument('-pdbs_folder', type=str, help=' directory where the pdbs to process exist. (Must be labeled *.pdb) ')
  parser.add_argument('-eval_method',choices=['nobu','tj','none'], default='nobu', help='evaluation method <nobu>,<tj>,<none> only valid options')

  args = parser.parse_args()

  pdbs = glob.glob(str(args.pdbs_folder)+"/*.pdb")
  if (len(pdbs) == 0):
    if(args.pdbs):
      pdbs = args.pdbs
    else:
      print("NO PDBS ENTERED")
      sys.exit()

  cmd_list = [( [pdb,args] ) for pdb in pdbs ]
  myPool    = Pool( processes = args.n_proc )
  myResults = myPool.map( generate_fragments_single_arg,cmd_list )
