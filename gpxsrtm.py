#!/usr/bin/python3
# -*- coding: utf-8 -*-
prog_ver = 'gpxsrtm v1.6 Copyright (c) 2019-2020 Matjaz Rihtar'
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

import gpxpy
from gpxpy.gpx import GPXBounds, GPXWaypoint

# latitude/longitude in GPX files is always in WGS84 datum
# WGS84 defines the Earth semi-major axis as 6378.137 km
EARTH_RADIUS = 6378.137 * 1000

gdata = None

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

# -----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D # needed for 3d projection

def plot_tile(data, title=None):
  try:
    plt.rc('figure', figsize=(9, 9))
    D = np.array(data)
    cmap = cm.terrain
    cmap.set_bad(color='black')
    plt.imshow(D, origin='lower', cmap=cmap)
    if title is not None:
      plt.title(title)
    plt.show()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
# plot_tile

def plot_elev(X, Y, Z, title=None):
  try:
    Zm = np.ma.masked_invalid(Z)
    XX, YY = np.meshgrid(X, Y)
    X1 = XX[~Zm.mask]
    Y1 = YY[~Zm.mask]
    Z1 = Zm[~Zm.mask]
    
    # interpolate missing points (nearest, linear, cubic)
    GD = griddata((X1, Y1), Z1.ravel(), (XX, YY), method='cubic')

    # make plot smoother (linear, cubic, quintic)
    #spline = interp2d(X, Y, GD, kind='linear')
    #Xi = np.linspace(X.min(), X.max(), 100)
    #Yi = np.linspace(Y.min(), Y.max(), 100)
    #Zi = spline(Xi, Yi)

    plt.rc('figure', figsize=(10,10))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ticks = np.linspace(GD.min(), GD.max(), 1000, endpoint=True)

    #ax.contour(Xi, Yi, Zi, ticks[1:-1], cmap='coolwarm')
    ax.contour(X, Y, GD, ticks[1:-1], cmap='coolwarm')
    
    ax.scatter3D(X1, Y1, Z1, color='black')

    ax.ticklabel_format(useOffset=False)
    if title is not None:
      plt.title(title) # fontsize=16

    plt.show()
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
# plot_elev

# -----------------------------------------------------------------------------
# Haversine distance in meters
# See http://www.movable-type.co.uk/scripts/latlong.html

def rad(x):
  return x * math.pi / 180

def distance(lat1, lon1, lat2, lon2):
  delta_lat = rad(lat1 - lat2)
  delta_lon = rad(lon1 - lon2)
  lat1 = rad(lat1)
  lat2 = rad(lat2)

  a = pow(math.sin(delta_lat/2), 2) + \
      math.cos(lat1) * math.cos(lat2) * pow(math.sin(delta_lon/2), 2)
  c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
  d = EARTH_RADIUS * c
  return d
# distance

# -----------------------------------------------------------------------------
def get_elev_idw(point):
  elev = None
  try:
    lat = point.latitude
    lon = point.longitude

    inc = 1/3600

    lat_deg = math.floor(lat)
    lon_deg = math.floor(lon)
    #print('lat = {:.9f}'.format(lat))
    #print('lon = {:.9f}'.format(lon))
    lat_sec = lat % 1
    lon_sec = lon % 1
    #print('lat_sec = {:.9f}'.format(lat_sec))
    #print('lon_sec = {:.9f}'.format(lon_sec))
    lat_idx = int(lat_sec * 3600)
    lon_idx = int(lon_sec * 3600)
    #print('lat_idx = {}'.format(lat_idx))
    #print('lon_idx = {}'.format(lon_idx))

    data = gdata[lat_deg][lon_deg]

    lat1 = lat_deg + lat_idx * inc
    lon1 = lon_deg + lon_idx * inc
    elev1 = data[lat_idx][lon_idx]
    #print('Base point at ({:.9f},{:.9f}), elev: {:.3f}'.format(lat1, lon1, elev1))

    lat2 = lat_deg + lat_idx * inc
    lon2 = lon_deg + (lon_idx + 1) * inc
    elev2 = data[lat_idx][lon_idx+1]

    lat3 = lat_deg + (lat_idx + 1) * inc
    lon3 = lon_deg + lon_idx * inc
    elev3 = data[lat_idx+1][lon_idx]

    lat4 = lat_deg + (lat_idx + 1) * inc
    lon4 = lon_deg + (lon_idx + 1) * inc
    elev4 = data[lat_idx+1][lon_idx+1]

    # Inverse Distance Weighting (IDW) interpolation
    w = 0
    elev = 0

    dist1 = distance(lat, lon, lat1, lon1)
    if dist1 == 0.0: dist1 = 0.01
    #print('dist1: {:.3f}, elev: {:.3f}'.format(dist1, elev1))
    w += 1/dist1
    elev += elev1/dist1

    dist2 = distance(lat, lon, lat2, lon2)
    if dist2 == 0.0: dist2 = 0.01
    #print('dist2: {:.3f}, elev: {:.3f}'.format(dist2, elev2))
    w += 1/dist2
    elev += elev2/dist2

    dist3 = distance(lat, lon, lat3, lon3)
    if dist3 == 0.0: dist3 = 0.01
    #print('dist3: {:.3f}, elev: {:.3f}'.format(dist3, elev3))
    w += 1/dist3
    elev += elev3/dist3

    dist4 = distance(lat, lon, lat4, lon4)
    if dist4 == 0.0: dist4 = 0.01
    #print('dist4: {:.3f}, elev: {:.3f}'.format(dist4, elev4))
    w += 1/dist4
    elev += elev4/dist4

    if w > 0:
      elev = elev/w

    s = '{:.2f}'.format(elev) # don't use round()
    elev = float(s)

    #X = np.array([lon1, lon, lon2])
    #Y = np.array([lat1, lat, lat3])
    #Z = np.array([[elev1, np.nan, elev2], [np.nan, elev, np.nan], [elev3, np.nan, elev4]])
    #plot_elev(X, Y, Z, 'IDW')

    diff = elev - point.elevation
    #print('Point at ({:.9f},{:.9f}), srtm_idw: {:.3f}, diff: {:.3f}'.format(lat, lon, elev, diff))
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return elev
# get_elev_idw

# -----------------------------------------------------------------------------
def get_elev_bilin(point):
  elev = None
  try:
    lat = point.latitude
    lon = point.longitude

    inc = 1/3600

    lat_deg = math.floor(lat)
    lon_deg = math.floor(lon)
    #print('lat = {:.9f}'.format(lat))
    #print('lon = {:.9f}'.format(lon))
    lat_sec = lat % 1
    lon_sec = lon % 1
    #print('lat_sec = {:.9f}'.format(lat_sec))
    #print('lon_sec = {:.9f}'.format(lon_sec))
    lat_idx = int(lat_sec * 3600)
    lon_idx = int(lon_sec * 3600)
    #print('lat_idx = {}'.format(lat_idx))
    #print('lon_idx = {}'.format(lon_idx))

    data = gdata[lat_deg][lon_deg]

    lat1 = lat_deg + lat_idx * inc
    lon1 = lon_deg + lon_idx * inc
    elev1 = data[lat_idx][lon_idx]
    #print('Base point at ({:.9f},{:.9f}), elev: {:.3f}'.format(lat1, lon1, elev1))

    lat2 = lat_deg + lat_idx * inc
    lon2 = lon_deg + (lon_idx + 1) * inc
    elev2 = data[lat_idx][lon_idx+1]

    lat3 = lat_deg + (lat_idx + 1) * inc
    lon3 = lon_deg + lon_idx * inc
    elev3 = data[lat_idx+1][lon_idx]

    lat4 = lat_deg + (lat_idx + 1) * inc
    lon4 = lon_deg + (lon_idx + 1) * inc
    elev4 = data[lat_idx+1][lon_idx+1]

    # bilinear interpolaton
    R1 = (lon2 - lon)/(lon2 - lon1)*elev1 + (lon - lon1)/(lon2 - lon1)*elev2
    R2 = (lon2 - lon)/(lon2 - lon1)*elev3 + (lon - lon1)/(lon2 - lon1)*elev4
    elev = (lat3 - lat)/(lat3 - lat1)*R1 + (lat - lat1)/(lat3 - lat1)*R2

    s = '{:.2f}'.format(elev) # don't use round()
    elev = float(s)

    #X = [lon1, lon2]
    #Y = [lat1, lat3]
    #Z = [elev1, elev2, elev3, elev4]
    #spline = interp2d(X, Y, Z, kind='linear')
    #Xp = [lon]
    #Yp = [lat]
    #Zp = spline(Xp, Yp) # should be the the same as above

    #X = np.array([lon1, lon, lon2])
    #Y = np.array([lat1, lat, lat3])
    #Z = np.array([[elev1, np.nan, elev2], [np.nan, elev, np.nan], [elev3, np.nan, elev4]])
    #plot_elev(X, Y, Z, 'Bilinear')

    diff = elev - point.elevation
    #print('Point at ({:.9f},{:.9f}), srtm_bil: {:.3f}, diff: {:.3f}'.format(lat, lon, elev, diff))
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return elev
# get_elev_bilin

# -----------------------------------------------------------------------------
def procfile(gpxname):
  global gdata, interp, datadir

  rc = 0
  try:
    fd = open(gpxname, 'r')
    gpx = gpxpy.parse(fd)
    fd.close()

    abounds = GPXBounds(90, -90, 180, -180)
    for track in gpx.tracks:
      bounds = track.get_bounds()
      if bounds.min_latitude < abounds.min_latitude:
        abounds.min_latitude = bounds.min_latitude
      if bounds.max_latitude > abounds.max_latitude:
        abounds.max_latitude = bounds.max_latitude
      if bounds.min_longitude < abounds.min_longitude:
        abounds.min_longitude = bounds.min_longitude
      if bounds.max_longitude > abounds.max_longitude:
        abounds.max_longitude = bounds.max_longitude

    min_lat = math.floor(abounds.min_latitude)
    max_lat = math.floor(abounds.max_latitude)
    min_lon = math.floor(abounds.min_longitude)
    max_lon = math.floor(abounds.max_longitude)
    #sys.stderr.write('min_latitude: {} -> {}\n'.format(abounds.min_latitude, min_lat))
    #sys.stderr.write('max_latitude: {} -> {}\n'.format(abounds.max_latitude, max_lat))
    #sys.stderr.write('min_longitude: {} -> {}\n'.format(abounds.min_longitude, min_lon))
    #sys.stderr.write('max_longitude: {} -> {}\n'.format(abounds.max_longitude, max_lon))
    if min_lat == max_lat and min_lon == max_lon:
      sys.stderr.write('Bounding box: {} {} (1 SRTM file)\n'.format(min_lat, min_lon))
    else:
      sys.stderr.write('Bounding box: {} {} -> {} {} ({} SRTM files)\n'.format(min_lat, min_lon, max_lat, max_lon, \
                       (max_lat - min_lat + 1) * (max_lon - min_lon + 1)))

    gdata = {}
    for lat in range(min_lat, max_lat+1):
      gdata[lat] = {}
      for lon in range(min_lon, max_lon+1):
        if lat < 0:
          if lon < 0: fname = 's{:02d}_w{:03d}_1arc_v3.pickle'.format(abs(lat), abs(lon))
          else: fname = 's{:02d}_e{:03d}_1arc_v3.pickle'.format(abs(lat), lon)
        else:
          if lon < 0: fname = 'n{:02d}_w{:03d}_1arc_v3.pickle'.format(lat, abs(lon))
          else: fname = 'n{:02d}_e{:03d}_1arc_v3.pickle'.format(lat, lon)

        if datadir == '<prog>{}data'.format(os.sep):
          fpath = '{}data{}{}'.format(where, os.sep, fname)
        else:
          fpath = '{}{}{}'.format(datadir, os.sep, fname)
        if not os.path.isfile(fpath):
          raise ImportError('missing {}'.format(fname))

        sys.stderr.write('Reading {}\n'.format(fpath))
        fd = open(fpath, 'rb')
        data = pickle.load(fd)
        fd.close()

        if plot:
          plot_tile(data, fname.replace('.pickle', ''))

        gdata[lat][lon] = data

    #pprint(('gdata', gdata))

    if interp == 'bilinear':
      get_elev = get_elev_bilin; ext = 'bil'
    else:
      get_elev = get_elev_idw; ext = 'idw'

    for track in gpx.tracks:
      for segment in track.segments:
        for point in segment.points:
          #print('Point at ({:.9f},{:.9f}), elev: {:.3f}'.format(point.latitude, point.longitude, point.elevation))
          point.elevation = get_elev(point)
          #print()

    gpxname = gpxname.replace('.gpx', '-{}.gpx'.format(ext))
    fd = open(gpxname, 'w')
    sys.stderr.write('Writing {}\n'.format(gpxname))
    fd.write(gpx.to_xml())
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
  global where, prog, interp, datadir, plot

# argv = ['srtm.py', 'sample.gpx']
  where = ntdirname(argv[0])
  prog = ntbasename(argv[0]).replace('.py', '').replace('.PY', '')

  parser = argparse.ArgumentParser(description='Provides elevation for specified geographic coordinates',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inp_files', metavar='<input.gpx>', nargs='+',
                       help='Input GPX file(s)')
  parser.add_argument('-i', metavar='<interp>', dest='interp',
                      choices=['bilinear', 'idw'], default='bilinear',
                      help='Interpolation type: bilinear or idw')
  parser.add_argument('-d', metavar='<datadir>', default='<prog>{}data'.format(os.sep),
                      dest='datadir', help='SRTM data directory')
  parser.add_argument('-p', action='store_true', default=False, dest='plot',
                       help='Plot read SRTM data')

  args = parser.parse_args()
  interp = args.interp
  datadir = args.datadir
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
