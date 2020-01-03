#!/usr/bin/python3
# -*- coding: utf-8 -*-
prog_ver = 'json2pickle v1.0 Copyright (c) 2019-2020 Matjaz Rihtar'
import sys, os
import ntpath
import traceback

import pickle, json

# -----------------------------------------------------------------------------
def ntdirname(path):
  try:
    head, tail = ntpath.split(path)
    dirname = head or ntpath.dirname(head)
  except: dirname = '.'
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
    sys.exit('Usage: {} <file.json>'.format(prog))

  try:
    jpath = sys.argv[1]
    fd = open(jpath, 'r', encoding='utf-8')
    sys.stderr.write('Reading {}\n'.format(jpath))
    data = json.load(fd)
    fd.close()

    ppath = jpath.replace('.json', '.pickle')
    if not ppath.endswith('.pickle'):
      ppath += '.pickle'

    if os.path.isfile(ppath):
      raise FileExistsError('File {} already exists'.format(ppath))

    fd = open(ppath, 'wb')
    sys.stderr.write('Writing {}\n'.format(ppath))
    pickle.dump(data, fd)
    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    errmsg = '{}({}): {}'.format(prog, exc_tb.tb_lineno, exc[-1].strip())
    sys.exit(errmsg)
