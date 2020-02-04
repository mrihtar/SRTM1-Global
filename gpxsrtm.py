#!/usr/bin/env python3
# -*- coding: utf-8 -*-
prog_ver = 'gpxsrtm v1.16 Copyright (c) 2019-2020 Matjaz Rihtar'
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
from scipy.interpolate import griddata, interp2d, Rbf

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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D # needed for 3d projection

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
    if title is not None: plt.title(title) # fontsize=16
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

def distance(lat1, lon1, lat2, lon2):
  delta_lat = math.radians(lat1 - lat2)
  delta_lon = math.radians(lon1 - lon2)
  lat1 = math.radians(lat1)
  lat2 = math.radians(lat2)

  a = pow(math.sin(delta_lat/2), 2) + \
      math.cos(lat1) * math.cos(lat2) * pow(math.sin(delta_lon/2), 2)
  c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
  d = EARTH_RADIUS * c
  return d
# distance

# -----------------------------------------------------------------------------
def get_elev(point, it='bil'):
  elev = None
  try:
    #print('----------------------------------------')
    #print('Processing "{}"'.format(point.name))
    lat = point.latitude
    lon = point.longitude
    #print('lat = {:.7f}'.format(lat))
    #print('lon = {:.7f}'.format(lon))

    inc = 1/3600; rinc = 1 - inc    # 0.0002778, 0.9997223
    inc2 = 1/7200; rinc2 = 1 - inc2 # 0.0001389, 0.9998612

    lat_deg = math.floor(lat)
    lon_deg = math.floor(lon)
    #print('lat_deg = {}'.format(lat_deg))
    #print('lon_deg = {}'.format(lon_deg))

    alat_sec = abs(lat) % 1
    alon_sec = abs(lon) % 1
    #print('alat_sec = {:.7f}'.format(alat_sec))
    #print('alon_sec = {:.7f}'.format(alon_sec))
    if source == 'alos':
      ofs = inc2
      if lat < 0:
        if alat_sec > rinc2:
          lat_deg -= 1
          #print('fixed lat_deg = {}'.format(lat_deg))
          alat_sec = 1 - inc2
        else:
          alat_sec = 1 - alat_sec - inc2
      else:
        if alat_sec < inc2:
          lat_deg -= 1
          #print('fixed lat_deg = {}'.format(lat_deg))
          alat_sec = 1 - inc2
        else:
          alat_sec = alat_sec - inc2

      if lon < 0:
        if alon_sec > rinc2:
          lon_deg -= 1
          #print('fixed lon_deg = {}'.format(lon_deg))
          alon_sec = 1 - inc2
        else:
          alon_sec = 1 - alon_sec - inc2
      else:
        if alon_sec < inc2:
          lon_deg -= 1
          #print('fixed lon_deg = {}'.format(lon_deg))
          alon_sec = 1 - inc2
        else:
          alon_sec = alon_sec - inc2
    else:
      ofs = 0
      if lat < 0:
        alat_sec = 1 - alat_sec
      if lon < 0:
        alon_sec = 1 - alon_sec

    lat_idx = int(alat_sec * 3600)
    lon_idx = int(alon_sec * 3600)
    #print('lat_idx = {}'.format(lat_idx))
    #print('lon_idx = {}'.format(lon_idx))

    data = gdata[lat_deg][lon_deg]

    lat1 = lat_deg + ofs + lat_idx * inc
    lon1 = lon_deg + ofs + lon_idx * inc
    elev1 = data[lat_idx][lon_idx]

    lat2 = lat_deg + ofs + lat_idx * inc
    lon2 = lon_deg + ofs + (lon_idx + 1) * inc
    elev2 = data[lat_idx][lon_idx+1]

    lat3 = lat_deg + ofs + (lat_idx + 1) * inc
    lon3 = lon_deg + ofs + lon_idx * inc
    elev3 = data[lat_idx+1][lon_idx]

    lat4 = lat_deg + ofs + (lat_idx + 1) * inc
    lon4 = lon_deg + ofs + (lon_idx + 1) * inc
    elev4 = data[lat_idx+1][lon_idx+1]

    if it == 'bil': 
      # bilinear interpolaton
      #print('Point1 at ({:.7f},{:.7f}), elev: {:.2f}'.format(lat1, lon1, elev1))
      #print('Point2 at ({:.7f},{:.7f}), elev: {:.2f}'.format(lat2, lon2, elev2))
      #print('Point3 at ({:.7f},{:.7f}), elev: {:.2f}'.format(lat3, lon3, elev3))
      #print('Point4 at ({:.7f},{:.7f}), elev: {:.2f}'.format(lat4, lon4, elev4))

      R1 = (lon2 - lon)/(lon2 - lon1)*elev1 + (lon - lon1)/(lon2 - lon1)*elev2
      R2 = (lon2 - lon)/(lon2 - lon1)*elev3 + (lon - lon1)/(lon2 - lon1)*elev4
      elev = (lat3 - lat)/(lat3 - lat1)*R1 + (lat - lat1)/(lat3 - lat1)*R2

      #X = [lon1, lon2]
      #Y = [lat1, lat3]
      #Z = [elev1, elev2, elev3, elev4]
      #spline = interp2d(X, Y, Z, kind='linear')
      #Xp = [lon]
      #Yp = [lat]
      #Zp = spline(Xp, Yp) # result should be the the same as above

    else: # it == 'idw'
      # Inverse Distance Weighting (IDW) interpolation
      w = 0
      elev = 0

      dist1 = distance(lat, lon, lat1, lon1)
      if dist1 == 0.0: dist1 = 0.01
      #print('Point1 at ({:.7f},{:.7f}), elev: {:.2f}, dist: {:.2f}'.format(lat1, lon1, elev1, dist1))
      w += 1/dist1
      elev += elev1/dist1

      dist2 = distance(lat, lon, lat2, lon2)
      if dist2 == 0.0: dist2 = 0.01
      #print('Point2 at ({:.7f},{:.7f}), elev: {:.2f}, dist: {:.2f}'.format(lat2, lon2, elev2, dist2))
      w += 1/dist2
      elev += elev2/dist2

      dist3 = distance(lat, lon, lat3, lon3)
      if dist3 == 0.0: dist3 = 0.01
      #print('Point3 at ({:.7f},{:.7f}), elev: {:.2f}, dist: {:.2f}'.format(lat3, lon3, elev3, dist3))
      w += 1/dist3
      elev += elev3/dist3

      dist4 = distance(lat, lon, lat4, lon4)
      if dist4 == 0.0: dist4 = 0.01
      #print('Point4 at ({:.7f},{:.7f}), elev: {:.2f}, dist: {:.2f}'.format(lat4, lon4, elev4, dist4))
      w += 1/dist4
      elev += elev4/dist4

      if w > 0:
        elev = elev/w

    s = '{:.2f}'.format(elev) # don't use round()
    elev = float(s)

    #X = np.array([lon1, lon, lon2])
    #Y = np.array([lat1, lat, lat3])
    #Z = np.array([[elev1, np.nan, elev2], [np.nan, elev, np.nan], [elev3, np.nan, elev4]])
    #plot_elev(X, Y, Z, it)

    diff = elev - point.elevation
    #print('Point at ({:.7f},{:.7f}), elev: {:.2f}, elev_{}: {:.2f}, diff: {:.2f}'.format(lat, lon, point.elevation, it, elev, diff))
  except:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    exc = traceback.format_exception_only(exc_type, exc_obj)
    name = sys._getframe().f_code.co_name
    errmsg = '{}({}): {}\n'.format(name, exc_tb.tb_lineno, exc[-1].strip())
    sys.stderr.write(errmsg)
  return elev
# get_elev

# -----------------------------------------------------------------------------
def procfile(gpxname):
  global gdata, interp, datadir

  rc = 0
  try:
    fd = open(gpxname, 'r', encoding='utf-8')
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

    #abounds.min_latitude = 45.99990
    #abounds.max_latitude = 46.00050
    #abounds.min_longitude = 13.00025
    #abounds.max_longitude = 13.00050
    #print('min_latitude: {:.7f}'.format(abounds.min_latitude))
    #print('max_latitude: {:.7f}'.format(abounds.max_latitude))
    #print('min_longitude: {:.7f}'.format(abounds.min_longitude))
    #print('max_longitude: {:.7f}'.format(abounds.max_longitude))

    inc = 1/3600; rinc = 1 - inc    # 0.0002778, 0.9997223
    inc2 = 1/7200; rinc2 = 1 - inc2 # 0.0001389, 0.9998612

    if source == 'alos':
      amin_lat_sec = abs(abounds.min_latitude) % 1
      #print('amin_lat_sec = {:.7f}'.format(amin_lat_sec))
      if abounds.min_latitude < 0:
        if amin_lat_sec > rinc2: abounds.min_latitude -= inc
      else:
        if amin_lat_sec < inc2: abounds.min_latitude -= inc

      amax_lat_sec = abs(abounds.max_latitude) % 1
      #print('amax_lat_sec = {:.7f}'.format(amax_lat_sec))
      if abounds.min_latitude < 0:
        if amax_lat_sec > rinc2: abounds.max_latitude -= inc
      else:
        if amax_lat_sec < inc2: abounds.max_latitude -= inc

      amin_lon_sec = abs(abounds.min_longitude) % 1
      #print('amin_lon_sec = {:.7f}'.format(amin_lon_sec))
      if abounds.min_longitude < 0:
        if amin_lon_sec > rinc2: abounds.min_longitude -= inc
      else:
        if amin_lon_sec < inc2: abounds.min_longitude -= inc

      amax_lon_sec = abs(abounds.max_longitude) % 1
      #print('amax_lon_sec = {:.7f}'.format(amax_lon_sec))
      if abounds.min_longitude < 0:
        if amax_lon_sec > rinc2: abounds.max_longitude -= inc
      else:
        if amax_lon_sec < inc2: abounds.max_longitude -= inc

    min_lat = math.floor(abounds.min_latitude)
    max_lat = math.floor(abounds.max_latitude)
    min_lon = math.floor(abounds.min_longitude)
    max_lon = math.floor(abounds.max_longitude)
    #print('min_latitude: {:.7f} -> {}'.format(abounds.min_latitude, min_lat))
    #print('max_latitude: {:.7f} -> {}'.format(abounds.max_latitude, max_lat))
    #print('min_longitude: {:.7f} -> {}'.format(abounds.min_longitude, min_lon))
    #print('max_longitude: {:.7f} -> {}'.format(abounds.max_longitude, max_lon))
    if min_lat == max_lat and min_lon == max_lon:
      sys.stderr.write('Bounding box: {} {} (1 data file)\n'.format(min_lat, min_lon))
    else:
      sys.stderr.write('Bounding box: {} {} -> {} {} ({} data files)\n'.format(min_lat, min_lon, max_lat, max_lon, \
                       (abs(max_lat - min_lat) + 1) * (abs(max_lon - min_lon) + 1)))

    gdata = {}
    for lat in range(min_lat, max_lat+1):
      alat = abs(lat)
      gdata[lat] = {}
      for lon in range(min_lon, max_lon+1):
        alon = abs(lon)
        if lat < 0:
          if lon < 0:
            if source == 'srtm':
              fname = 's{:02d}_w{:03d}_1arc_v3.pickle'.format(alat, alon)
            else: # alos
              fname = 'S{:03d}W{:03d}_AVE_EXT.pickle'.format(alat, alon)
          else:
            if source == 'srtm':
              fname = 's{:02d}_e{:03d}_1arc_v3.pickle'.format(alat, lon)
            else: # alos
              fname = 'S{:03d}E{:03d}_AVE_EXT.pickle'.format(alat, lon)
        else:
          if lon < 0:
            if source == 'srtm':
              fname = 'n{:02d}_w{:03d}_1arc_v3.pickle'.format(lat, alon)
            else: # alos
              fname = 'N{:03d}W{:03d}_AVE_EXT.pickle'.format(lat, alon)
          else:
            if source == 'srtm':
              fname = 'n{:02d}_e{:03d}_1arc_v3.pickle'.format(lat, lon)
            else: # alos
              fname = 'N{:03d}E{:03d}_AVE_EXT.pickle'.format(lat, lon)

        if datadir == '<prog>{}data'.format(os.sep):
          fpath = '{}data{}{}'.format(where, os.sep, fname)
        else:
          fpath = '{}{}{}'.format(datadir, os.sep, fname)
        if not os.path.isfile(fpath):
          raise ImportError('missing {}'.format(fpath))

        sys.stderr.write('Reading {}\n'.format(fpath))
        fd = open(fpath, 'rb')
        data = pickle.load(fd)
        fd.close()

        if plot:
          plot_tile(data, lon, 1/3600, lat, 1/3600, fname.replace('.pickle', ''))

        gdata[lat][lon] = data

    #pprint(('gdata', gdata))

    if interp == 'bilinear': it = 'bil'
    else: it = 'idw'

    for waypoint in gpx.waypoints:
      #print('Waypoint at ({:.7f},{:.7f}), elev: {:.2f}'.format(waypoint.latitude, waypoint.longitude, waypoint.elevation))
      waypoint.elevation = get_elev(waypoint, it)
      #print()

    for track in gpx.tracks:
      for segment in track.segments:
        for point in segment.points:
          #print('Point at ({:.7f},{:.7f}), elev: {:.2f}'.format(point.latitude, point.longitude, point.elevation))
          point.elevation = get_elev(point, it)
          #print()

    gpxname = gpxname.replace('.gpx', '-{}-{}.gpx'.format(source, it))
    fd = open(gpxname, 'w', encoding='utf-8')
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
  global where, prog, source, interp, datadir, plot

# argv = ['gpxsrtm.py', 'test.gpx']
  where = ntdirname(argv[0])
  prog = ntbasename(argv[0]).replace('.py', '').replace('.PY', '')

  parser = argparse.ArgumentParser(description='Provides elevation for specified geographic coordinates',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inp_files', metavar='<input.gpx>', nargs='+',
                       help='Input file(s) in GPX format')
  parser.add_argument('-s', metavar='<source>', dest='source',
                      choices=['srtm', 'alos'], default='srtm',
                      help='Data source: srtm or alos')
  parser.add_argument('-i', metavar='<interp>', dest='interp',
                      choices=['bilinear', 'idw'], default='bilinear',
                      help='Interpolation type: bilinear or idw')
  parser.add_argument('-d', metavar='<datadir>', default='<prog>{}data'.format(os.sep),
                      dest='datadir', help='SRTM/ALOS data directory')
  parser.add_argument('-p', action='store_true', default=False, dest='plot',
                       help='Plot read SRTM/ALOS data')

  args = parser.parse_args()
  source = args.source
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
