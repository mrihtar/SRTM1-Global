#!/usr/bin/env python3
# -*- coding: utf-8 -*-
prog_ver = 'prepare v1.12 Copyright (c) 2019-2020 Matjaz Rihtar'
# py_ver = sys.version_info.major
import sys, os, glob, re
import ntpath, argparse
import traceback
from pprint import pprint

import math
from struct import unpack
import pickle, json

import numpy as np
import scipy as sp
from scipy.interpolate import griddata, interp2d
from scipy.spatial import cKDTree

tiff_tags = {
    256: 'image width',
    257: 'image length',
    258: 'bits per sample',
    259: 'compression',
    262: 'photometric interpretation',
    273: 'strip offsets',
    277: 'samples per pixel',
    278: 'rows per strip',
    279: 'strip byte counts',
    282: 'X resolution',
    283: 'Y resolution',
    284: 'planar configuration',
    296: 'resolution unit',
    320: 'color map',
    339: 'sample format',
  33550: 'model pixel scale tag',
  33922: 'model tiepoint tag',
  34735: 'geokey directory tag',
  34736: 'geo double params tag',
  34737: 'geo ascii params tag',
  42112: 'GDAL metadata',
  42113: 'GDAL nodata'
}

tiff_types = {
   1: [1, 'byte'],
   2: [1, 'ascii'],
   3: [2, 'short'],
   4: [4, 'long'],
   5: [8, 'rational'],
   6: [1, 'sbyte'],
   7: [1, 'undefined'],
   8: [2, 'sshort'],
   9: [4, 'slong'],
  10: [8, 'srational'],
  11: [4, 'float'],
  12: [8, 'double']
}

gk_ids = {
  1024: 'GTModelTypeGeoKey',
  1025: 'GTRasterTypeGeoKey',
  2048: 'GeodeticCRSGeoKey',
  2049: 'GeodeticCitationGeoKey',
  2054: 'GeogAngularUnitsGeoKey',
  2057: 'EllipsoidSemiMajorAxisGeoKey',
  2059: 'EllipsoidInvFlatteningGeoKey'
}

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
def decode_smr16(x):
  if x == 0xffff:
    return None # no data
  elif x >= 0x8000:
    return 0x8000-x # signed magnitude
  return x
# decode_smr16

# -----------------------------------------------------------------------------
def load_dted(fpath):
  data = None; raster_type = None
  try:
    fd = open(fpath, 'rb')
    sys.stderr.write('Reading {}\n'.format(fpath))

    #
    # Parse UHL header
    #
    # check magic and version
    if fd.read(3) != b'UHL':
      raise SyntaxError('wrong magic')
    if fd.read(1) != b'1':
      raise SyntaxError('wrong version')

    longitude = fd.read(8).decode('utf-8')
    latitude = fd.read(8).decode('utf-8')

    lon_ival = int(fd.read(4))/10
    lat_ival = int(fd.read(4))/10

    tmp = fd.read(4)
    if tmp == b'NA  ': abs_vacc = None
    else: abs_vacc = int(tmp)

    usc = fd.read(3)
    uref = fd.read(12)

    num_lon = int(fd.read(4))
    num_lat = int(fd.read(4))
    data = np.empty([num_lon, num_lat])

    mult_acc = int(fd.read(1))

    reserved = fd.read(24)

    #
    # Skip other headers
    #
    fd.seek(648, 1)  # DSI
    fd.seek(2700, 1) # ACC

    #
    # Read data blocks
    #
    for ii in range(num_lon):
      # check magic of record
      magic_b = fd.read(1)
      magic = unpack('>B', magic_b)[0]
      if magic != 0xAA: # 170
        raise SyntaxError('wrong magic in data block')

      seq_b = fd.read(3)
      seq = unpack('>I', b'\x00' + seq_b)[0]

      lon_cnt_b = fd.read(2)
      lon_cnt = unpack('>H', lon_cnt_b)[0]
      if lon_cnt != ii:
        raise ValueError('unexpected longitude number: {}'.format(lon_cnt))

      lat_cnt_b = fd.read(2)
      lat_cnt = unpack('>H', lat_cnt_b)[0]
      if lat_cnt != 0:
        raise ValueError('latitude count not zero')

      # read elevations
      rowdata_b = fd.read(2*num_lat)

      # check values with checksum
      # (checksum is calculated as sum of unsigned bytes of whole row)
      checksum = unpack('>I', fd.read(4))[0]
      whole_row_b = magic_b + seq_b + lon_cnt_b + lat_cnt_b + rowdata_b
      rowsum = sum(unpack('>' + (8+2*num_lat)*'B', whole_row_b))
      if rowsum != checksum:
        raise ValueError('checksum failed on longitude {}, should be {:08X} but is {:08X}'.format(ii, checksum, rowsum))

      # add row to matrix, this time values are
      # interpreted as big endian 16-bit (signed magnitude)
      row = unpack('>' + num_lat*'H', rowdata_b)
      for lon in range(0, num_lon):
        data[lon][ii] = decode_smr16(row[lon])

    raster_type = 2 # PixelIsPoint (0,0)

    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
    data = None; raster_type = None
  return data, raster_type
# load_dted

# -----------------------------------------------------------------------------
def decode_sig16(x):
  if x == -32767 or x == -9999: # ALOS has -9999 for void
    return None # no data
  return x
# decode_sig16

# -----------------------------------------------------------------------------
def load_bil(fpath):
  data = None; raster_type = None
  try:
    fd = open(fpath, 'rb')
    sys.stderr.write('Reading {}\n'.format(fpath))

    # these should be read from .hdr
    num_lon = 3601
    num_lat = 3601
    data = np.empty([num_lon, num_lat])

    for ii in range(num_lon):
      # read elevations
      rowdata_b = fd.read(2*num_lat)

      row = unpack('<' + num_lat*'h', rowdata_b)
      data[num_lon-ii-1] = [ decode_sig16(x) for x in row ]

    raster_type = 2 # PixelIsPoint (0,0)

    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
    data = None; raster_type = None
  return data, raster_type
# load_bil

# -----------------------------------------------------------------------------
def load_tif(fpath):
  data = None; raster_type = None
  try:
    fd = open(fpath, 'rb')
    sys.stderr.write('Reading {}\n'.format(fpath))

    # check byte order and version
    byteord_b = fd.read(2)
    if byteord_b == b'II': le = True # little endian
    elif byteord_b == b'MM': le = False # big endian
    else: raise SyntaxError('wrong byte order')
    #print('Read byte order: {}'.format(str(byteord_b)))

    version_b = fd.read(2)
    if le: version = unpack('<H', version_b)[0]
    else: version = unpack('>H', version_b)[0]
    if version != 0x2A: # 42
      raise SyntaxError('wrong version')

    while True:
      # offset to Image File Directory (IFD)
      ifd_b = fd.read(4)
      if le: ifd = unpack('<I', ifd_b)[0]
      else: ifd = unpack('>I', ifd_b)[0]
      #print('Read IFD offset: {:08X}'.format(ifd))
      if ifd == 0:
        break

      fd.seek(ifd)

      # entry count
      entcnt_b = fd.read(2)
      if le: entcnt = unpack('<H', entcnt_b)[0]
      else: entcnt = unpack('>H', entcnt_b)[0]
      #print('Read IFD entry count: {}'.format(entcnt))

      # for each entry
      for ii in range(0, entcnt):
        # tag
        tag_b = fd.read(2)
        if le: tag = unpack('<H', tag_b)[0]
        else: tag = unpack('>H', tag_b)[0]
        #if tag in tiff_tags:
        #  print('\n({}) tag: {} {}'.format(ii, tag, tiff_tags[tag]))
        #else: 
        #  print('\n({}) tag: {} unknown'.format(ii, tag))

        # type code
        tipe_b = fd.read(2)
        if le: tipe = unpack('<H', tipe_b)[0]
        else: tipe = unpack('>H', tipe_b)[0]
        if tipe in tiff_types:
          tipe_len = tiff_types[tipe][0]
          #print('({}) type: {} {}'.format(ii, tipe, tiff_types[tipe][1]))
        else:
          tipe_len = None
          #print('({}) type: {} unknown'.format(ii, tipe))

        # count field
        count_b = fd.read(4)
        if le: count = unpack('<I', count_b)[0]
        else: count = unpack('>I', count_b)[0]
        #print('({}) count: {}'.format(ii, count))

        # data pointer/field
        offset_b = fd.read(4)
        if le: offset = unpack('<I', offset_b)[0]
        else: offset = unpack('>I', offset_b)[0]

        entlen = count * tipe_len
        #if entlen <= 4:
        #  print('({}) value: {}'.format(ii, offset))
        #else:
        #  print('({}) offset: {:08X}'.format(ii, offset))

        if tag == 256: num_lon = offset   # image width
        elif tag == 257: num_lat = offset # image length
        elif tag == 278: rows = offset    # rows per strip
        elif tag == 273: # strip offsets
          if entlen <= 4:
            strips = [ offset ]
          else:
            oldpos = fd.tell()
            fd.seek(offset)
            strips_b = fd.read(tipe_len * count)
            if le: strips = unpack('<' + count*'I', strips_b)
            else: strips = unpack('>' + count*'I', strips_b)
            fd.seek(oldpos)
        elif tag == 33922: # model tiepoint tag
          oldpos = fd.tell()
          fd.seek(offset)
          mttags_b = fd.read(tipe_len * count)
          if le: mttags = unpack('<' + count*'d', mttags_b)
          else: mttags = unpack('>' + count*'d', mttags_b)
          #print(mttags)
          fd.seek(oldpos)
        elif tag == 33550: # model pixel scale tag
          oldpos = fd.tell()
          fd.seek(offset)
          mpstags_b = fd.read(tipe_len * count)
          if le: mpstags = unpack('<' + count*'d', mpstags_b)
          else: mpstags = unpack('>' + count*'d', mpstags_b)
          #print(mpstags)
          fd.seek(oldpos)
        elif tag == 34735: # geokey directory tag
          oldpos = fd.tell()
          fd.seek(offset)
          gk_header_b = fd.read(tipe_len * 4)
          if le: gk_header = unpack('<' + 4*'H', gk_header_b)
          else: gk_header = unpack('>' + 4*'H', gk_header_b)

          gk_num = gk_header[3]
          #print('Read number of GKs: {}'.format(gk_num))
          for ii in range(gk_num):
            gk_key_b = fd.read(tipe_len * 4)
            if le: gk_key = unpack('<' + 4*'H', gk_key_b)
            else: gk_key = unpack('>' + 4*'H', gk_key_b)

            gk_id = gk_key[0]
            #if gk_id in gk_ids:
            #  print('\n(GK:{}) id: {} {}'.format(ii, gk_id, gk_ids[gk_id]))
            #else:
            #  print('\n(GK:{}) id: {} unknown'.format(ii, gk_id))
            gk_tagloc = gk_key[1]
            #print('(GK:{}) tagloc: {}'.format(ii, gk_tagloc))
            gk_count = gk_key[2]
            #print('(GK:{}) count: {}'.format(ii, gk_count))
            gk_offset = gk_key[3]
            #if gk_tagloc == 0:
            #  print('(GK:{}) value: {}'.format(ii, gk_offset))
            #else:
            #  print('(GK:{}) offset: {:08X}'.format(ii, gk_offset))

            if gk_id == 1025: # GTRasterTypeGeoKey
              # 1 = PixelIsArea (0.5,0.5)
              # 2 = PixelIsPoint (0,0)
              raster_type = gk_offset
              if raster_type == 1:
                sys.stderr.write('Raster type is PixelIsArea (0.5,0.5)\n')

          fd.seek(oldpos)

    data = np.empty([num_lon, num_lat])

    nr = 0
    for ii in range(len(strips)):
      offset = strips[ii]
      fd.seek(offset)

      for jj in range(rows):
        # read elevations
        rowdata_b = fd.read(2*num_lat)
        if le: row = unpack('<' + num_lat*'h', rowdata_b)
        else: row = unpack('>' + num_lat*'h', rowdata_b)
        data[num_lon-nr-1] = [ decode_sig16(x) for x in row ]
        nr += 1

    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
    data = None; raster_type = None
  return data, raster_type
# load_tif

# -----------------------------------------------------------------------------
def griddata2d(X, Y, Z, fill=None):
  Zn = Z
  try:
    # filter out nan points
    Zm = np.ma.masked_invalid(Z)
    XX, YY = np.meshgrid(X, Y)
    X1 = XX[~Zm.mask]
    Y1 = YY[~Zm.mask]
    Z1 = Zm[~Zm.mask]

    # interpolate nan points
    sys.stderr.write('Interpolating with griddata\n')
    if fill is None:
      Zn = griddata((X1, Y1), Z1.ravel(), (XX, YY), method='linear') # nearest, linear, cubic
    else:
      Zn = griddata((X1, Y1), Z1.ravel(), (XX, YY), method='linear', fill_value=fill)
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return Zn
# griddata2d

# -----------------------------------------------------------------------------
def gridknn2d(X, Y, Z):
  Zn = Z
  try:
    # list of all Z indices
    XYitems = np.indices(Z.shape).T
    XYitems = XYitems.reshape((-1, XYitems.shape[-1]))

    # list of Z nan indices
    nanitems = np.argwhere(np.isnan(Z))

    # list of Z non-nan indices
    dims = np.maximum(XYitems.max(0), nanitems.max(0))+1
    a = np.ravel_multi_index(XYitems.T, dims)
    b = np.ravel_multi_index(nanitems.T, dims)
    XY = XYitems[~np.in1d(a, b)]

    sys.stderr.write('Creating K-D Tree\n')
    kdtree = cKDTree(XY, leafsize=32)

    Zn = np.copy(Z)

    # calculate wighted average knn value for each nan value in Z
    for ind in nanitems:
      x = ind[0]; y = ind[1]
      # query 12 neighbours, exact values, euclidean distance
      dist, idx = kdtree.query(ind, k=12, eps=0, p=2, distance_upper_bound=np.inf)
      dist[dist == 0] = 1.0 # nearest distance can be 0!

      Zxy = XY[idx]
      # get values from Z via query-returned indices
      values = np.take(Z, np.ravel_multi_index(Zxy.T, Z.shape))
      # calculate weighted average knn value
      value = np.average(values, weights=1.0/dist)

      # replace nan value with weighted average knn value
      Zn[x,y] = value
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return Zn
# gridknn2d

# -----------------------------------------------------------------------------
def fill_missing(lat0, lon0, Z, raster_type):
  global no_interp

  Zn = None; changed = False
  try:
    size = len(Z[0])
    sys.stderr.write('Starting point: {} {}, size: {} x {}\n'.format(lat0, lon0, size, size))

    if raster_type == 1: # PixelIsArea (0.5,0.5)
      endp = False
    else: # PixelIsPoint (0,0)
      endp = True # srtm last line is repeated as first line in next tile

    if lon0 < 0:
      alon0 = abs(lon0)
      X = np.linspace(alon0, alon0-1, num=size, endpoint=endp, dtype=np.float_)
    else:
      X = np.linspace(lon0, lon0+1, num=size, endpoint=endp, dtype=np.float_)
    #print('X = {}'.format(X))
    if lat0 < 0:
      alat0 = abs(lat0)
      Y = np.linspace(alat0, alat0-1, num=size, endpoint=endp, dtype=np.float_)
    else:
      Y = np.linspace(lat0, lat0+1, num=size, endpoint=endp, dtype=np.float_)
    #print('Y = {}'.format(Y))

    Z = np.array(Z).astype(np.float_)
    #Z = np.array(np.random.uniform(0.0, 3000.0, (size, size))).astype(np.float_)
    #Z[0][0] = np.nan
    #Z[0][1] = np.nan
    #Z[2][2] = np.nan

    Zn = Z

    nanitems = np.isnan(Zn).sum()
    if nanitems > 0 and not no_interp:
      sys.stderr.write('{} NaN values, interpolating\n'.format(nanitems))

      # interpolate points
      Zn = griddata2d(X, Y, Zn)
      changed = True
    
      nanitems = np.isnan(Zn).sum()

    if nanitems > 0:
      sys.stderr.write('{} NaN values, extrapolating\n'.format(nanitems))

      # extrapolate points
      Zn = gridknn2d(X, Y, Zn)
      changed = True
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return Zn, changed
# fill_missing

# -----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

def plot_tile(data, x0, xinc, y0, yinc, title=None, aspect=None):
  try:
    plt.rc('figure', figsize=(9, 9))
    D = np.array(data)
    fig, ax = plt.subplots()
    cmap = cm.terrain
    cmap.set_bad(color='black')
    im = ax.imshow(D, origin='lower', cmap=cmap)
    if aspect is not None:
      ax.set_aspect(aspect)
    #else:
    #  mercator_aspect = 1 / math.cos(math.radians(central_lat))
    #  ax.set_aspect(mercator_aspect)
    plt.colorbar(im, fraction=0.0457, pad=0.04, label='Elevation [m]')

    xmin = 0; xmax = len(D[0])
    xv = []; xvs = []
    xstep = int((xmax - xmin) / 8)
    if x0 < 0:
      tmp = xmin; xmin = xmax + 1; xmax = tmp
      xstep = -xstep
    for x in range(xmin, xmax, xstep):
      xv.append(x)
      xvs.append('{:.1f}'.format(x0 + x * xinc))
    plt.xticks(xv, xvs)
    plt.xlabel('Longitude')

    ymin = 0; ymax = len(D)
    yv = []; yvs = []
    ystep = int((ymax - ymin) / 5)
    if y0 < 0:
      tmp = ymin; ymin = ymax + 1; ymax = tmp
      ystep = -ystep
    for y in range(ymin, ymax, ystep):
      yv.append(y)
      yvs.append('{:.1f}'.format(y0 + y * yinc))
    plt.yticks(yv, yvs)
    plt.ylabel('Latitude')

    if title is not None: plt.title(title)
    plt.show()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
# plot_tile

# -----------------------------------------------------------------------------
def procfile(fpath):
  rc = 0
  try:
    fname = ntbasename(fpath)

    data = None
    if fname.endswith('.dt2'):
      data, raster_type = load_dted(fpath)
      ppath = fpath.replace('.dt2', '.pickle')
    elif fname.endswith('.bil'):
      data, raster_type = load_bil(fpath)
      ppath = fpath.replace('.bil', '.pickle')
    elif fname.endswith('.tif'):
      data, raster_type = load_tif(fpath)
      ppath = fpath.replace('.tif', '.pickle')
    else:
      raise NotImplementedError('unknown file type')

    if data is None:
      raise ImportError('error loading {}'.format(fpath))

    lat0 = 0; lon0 = 0

    fn = fname.split('_')
    if fn[0].startswith('n') or fn[0].startswith('s'): # SRTM1
      lat0 = int(fn[0][1:])
      if fn[0].startswith('s'): lat0 = -lat0
    if fn[1].startswith('e') or fn[1].startswith('w'): # SRTM1
      lon0 = int(fn[1][1:])
      if fn[1].startswith('w'): lon0 = -lon0
    if fn[0].startswith('N') or fn[0].startswith('S'): # ALOS
      lat0 = int(fn[0][1:4])
      if fn[0][0].startswith('S'): lat0 = -lat0
      lon0 = int(fn[0][5:8])
      if fn[0][4].startswith('W'): lon0 = -lon0

    if plot:
      plot_tile(data, lon0, 1/3600, lat0, 1/3600, '{} (raw)'.format(fname))

    data, changed = fill_missing(lat0, lon0, data, raster_type)

    if plot and changed:
      plot_tile(data, lon0, 1/3600, lat0, 1/3600, '{} (interpolated)'.format(fname))

    sys.stderr.write('Writing {}\n'.format(ppath))
    fd = open(ppath, 'wb')
    pickle.dump(data.tolist(), fd)
    fd.close()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
    rc = 1
  return rc
# procfile

# =============================================================================
def main(argv):
  global where, prog, no_interp, plot

# argv = ['prepare.py', 'data\s17_w068_1arc_v3.dt2']
  where = ntdirname(argv[0])
  prog = ntbasename(argv[0]).replace('.py', '').replace('.PY', '')

  parser = argparse.ArgumentParser(description='Prepares SRTM/ALOS data files for processing',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inp_files', metavar='<data_file>', nargs='+',
                       help='Input file(s) in DTED/BIL/TIFF format')
  parser.add_argument('-noi', action='store_true', default=False, dest='no_interp',
                       help='Don\'t do griddata interpolation')
  parser.add_argument('-p', action='store_true', default=False, dest='plot',
                       help='Plot read SRTM/ALOS data')

  args = parser.parse_args()
  no_interp = args.no_interp
  plot = args.plot

  rc = 0

  nfiles = 0
  for arg in args.inp_files:
    files = glob.glob(arg)
    for file in files:
      nfiles += 1
      rc += procfile(file)

  if nfiles == 0:
    sys.exit('No files found')

  return rc
# main

# -----------------------------------------------------------------------------
if __name__ == '__main__':
  rc = main(sys.argv)
  sys.exit(rc)
