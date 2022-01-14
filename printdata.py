#!/usr/bin/env python3
# -*- coding: utf-8 -*-
prog_ver = 'printdata v1.1 Copyright (c) 2019-2022 Matjaz Rihtar'
import sys, os, glob
import ntpath, argparse
import traceback
from pprint import pprint

import pickle, json

# -----------------------------------------------------------------------------
def ntdirname(path):
  try:
    head, tail = ntpath.split(path)
    dirname = head or ntpath.dirname(head)
  except: dirname = ''
  if len(dirname) == 0: dirname = '.'
  if dirname.endswith(os.sep):
    return dirname
  else:
    return dirname + os.sep
# ntdirname

def ntbasename(path):
  try:
    head, tail = ntpath.split(path)
    basename = tail or ntpath.basename(head)
  except: basename = ''
  return basename
# ntbasename

# -----------------------------------------------------------------------------
def print_matrix(data, fd=sys.stdout, rev=False, hdr=False):
  sdata = [[str(val) for val in row] for row in data]
  if hdr:
    n = len(data[0])
    cols = []
    for ii in range(n):
      cols.append('{}'.format(ii+1))
    sdata.append(cols)
    lens = [max(map(len, col)) for col in zip(*sdata)]
    sdata = sdata[:-1]
  else:
    lens = [max(map(len, col)) for col in zip(*sdata)]
  fmt = '\t'.join('{{:>{}}}'.format(x) for x in lens)
  table = [fmt.format(*row) for row in sdata]
  if rev: table = reversed(table)
  if hdr:
    header = fmt.format(*cols)
    print(header, file=fd)
  print('\n'.join(table), file=fd)
# print_matrix

# =============================================================================
if __name__ == '__main__':
  global prog

# argv = ['printdata.py', 'data\s17_w068_1arc_v3.pickle']
  where = ntdirname(sys.argv[0])
  prog = ntbasename(sys.argv[0]).replace('.py', '').replace('.PY', '')

  parser = argparse.ArgumentParser(description='Prints data from a pickle as a matrix',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inp_files', metavar='<input.pickle>', nargs='+',
                       help='Input file(s) in pickle format')
  parser.add_argument('-c', action='store_true', default=False, dest='header',
                       help='Print column header')
  parser.add_argument('-r', action='store_true', default=False, dest='reverse',
                       help='Print reversed matrix')

  args = parser.parse_args()
  header = args.header
  reverse = args.reverse

  try:
    nfiles = 0
    for arg in args.inp_files:
      files = glob.glob(arg)
      for file in files:
        nfiles += 1

        fd = open(file, 'rb')
        sys.stderr.write('Reading {}\n'.format(file))
        data = pickle.load(fd)
        fd.close()

        print_matrix(data, rev=reverse, hdr=header)

    if nfiles == 0:
      raise FileNotFoundError('No files found')
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    errmsg = '{}({}): {}'.format(prog, exc_tb.tb_lineno, exc[-1].strip())
    sys.exit(errmsg)
