#!/usr/bin/env python3
# -*- coding: utf-8 -*-
prog_ver = 'pickle2json v1.1 Copyright (c) 2019-2020 Matjaz Rihtar'
import sys, os
import ntpath
import traceback

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

# =============================================================================
if __name__ == '__main__':
  prog = ntbasename(sys.argv[0]).replace('.py', '').replace('.PY', '')
  if len(sys.argv) < 2:
    sys.exit('Usage: {} <file.pickle>'.format(prog))

  try:
    ppath = sys.argv[1]
    fd = open(ppath, 'rb')
    sys.stderr.write('Reading {}\n'.format(ppath))
    data = pickle.load(fd)
    fd.close()

    jpath = ppath.replace('.pickle', '.json')
    if not jpath.endswith('.json'):
      jpath += '.json'

    if os.path.isfile(jpath):
      raise FileExistsError('File {} already exists'.format(jpath))

    fd = open(jpath, 'w', encoding='utf-8')
    sys.stderr.write('Writing {}\n'.format(jpath))
    json.dump(data, fd)
    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    errmsg = '{}({}): {}'.format(prog, exc_tb.tb_lineno, exc[-1].strip())
    sys.exit(errmsg)
