#!/usr/bin/env python
# Author: Pietro Milillo - Oct 26 2013 
# This Script create a .grd file given an unwrapped interferogram and a DEM( the .grd file is ready to be plotted with GMT then  )
# modified to use filtered coherence for mask and adjusted coordinates EJF 2014/8/29
# modified to also convert masked.unw.geo if present and correctly detrend EJF 2014/12/10
# added conversion of los.rdr.geo to .grd for loading into CSI EJF 2014/2/16
# this version requires GIAnT so has to run in Python 2.x
# renamed to "unw2grd.py" EJF 2015/03/05
# modified to check Python version and ISCE version EJF 2015/03/09

from __future__ import print_function

import numpy as np
import os
import sys

if sys.version_info[0] > 2:
    sys.stderr.write('This version only works with Python 2\n')
    exit(1)

import isce
from iscesys.Parsers.FileParserFactory import createFileParser
try:
    from iscesys.ImageApi import DataAccessor as DA
except:
    sys.stderr.write('This version requires Python 2 ISCE DataAcessor ImageApi\n')
    exit(1)

from isceobj.Util import key_of_same_content
import math
import cmath as c
import tsinsar as ts 

import lxml.objectify as OB


# CSK wavelength hard-coded, but updated below from insarProc.xml if available
wavelength = 0.0312283810417
Pi = 3.141592653589
# unfiltered cohTH= 0.19
cohTH= 0.4

DEMfile = 'dem.crop'
name='filt_topophase.unw.geo'
nameM='masked.unw.geo' # assume same sizes
#coh='topophase.cor.geo'
coh='phsig.cor.geo'
cohxml= coh + '.xml'
los = 'los.rdr.geo'

# reading parameters from the xml file, first find file
if os.path.exists(name):
    imagexml= name + '.xml'
elif os.path.exists(nameM):
    imagexml= nameM + '.xml'

startLat = []
endLat = []
deltaLat = []
startLon = []
endLon = []
deltaLon = []
length = []
width = []
names = []

ext = None
dataType = None
width = None
PA = createFileParser('xml')
dictNow, dictFact, dictMisc = PA.parse(imagexml)
numBands = 0

numBands = key_of_same_content('number_bands', dictNow)[1]
bandsType = key_of_same_content('data_type', dictNow)[1]

coordinate1  =  key_of_same_content('coordinate1',dictNow)[1]
width = (int(key_of_same_content('size',coordinate1)[1]))
startLon=(float(key_of_same_content('startingValue',coordinate1)[1]))
deltaLon=(float(key_of_same_content('delta',coordinate1)[1]))


coordinate2  =  key_of_same_content('coordinate2',dictNow)[1]
length = (int(key_of_same_content('size',coordinate2)[1]))

startLat=(float(key_of_same_content('startingValue',coordinate2)[1]))
deltaLat=(float(key_of_same_content('delta',coordinate2)[1]))

# ISCE uses grid-centered coordinates, so need to subtract one
endLat= startLat + (deltaLat * (length-1) )
endLon= startLon + (deltaLon * (width-1) )

#print('LAMBDA: ' + str(wavelength)
print('WIDTH:  ' + str(width))
print('LENGTH: ' + str(length))
#print('LAT:    ' + str(startLat) + ' ' + str(endLat))
#print('LON:    ' + str(startLon) + ' ' + str(endLon))

# get other metadata from insarProc.xml, if available (modified from prepGiant.py
if os.path.exists('insarProc.xml'):
     exampleXML = 'insarProc.xml'
     fid = open(exampleXML,'r')
     xObj = OB.fromstring(fid.read())
     fid.close()

     wavelength = float(xObj.master.wavelength)
     heading = float(xObj.runGeocode.inputs.PEG_HEADING) * 180.0 / np.pi	

print('radar wavelength ', wavelength)

if os.path.exists(name):

    # Set to nan pixels with low coherence for unmasked phase

    phs = ts.load_mmap(name, width, length, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)
    phs = np.array(phs)

    COH = ts.load_mmap(coh, width, length, quiet=True, map='BIL', nchannels=1, channel=1, conv=False)
    # phsig is only one channel
    COH = np.array(COH)
    COH[COH > cohTH ] = 1
    COH[COH <= cohTH ] = 0

    phs = COH*phs
    phs[phs == 0]= np.nan
    phs.tofile(name + '.phs')

    command = 'xyz2grd ' + name + '.phs ' + '-R' + str(startLon) + '/' + str(endLon) + '/' + str(endLat) + '/' + str(startLat) + ' -I' + str(width) + '+/' + str(length) + '+ -ZTLf -N0.0 -G' + name + '.grd -V'
    print (command)
    os.system(command)
    command = 'grdmath ' + name + '.grd ' + str((wavelength/(4 *Pi))) + ' MUL = ' + name + '.math.grd'
    print (command)
    os.system(command)
    command = 'grdtrend ' + name + '.math.grd' + ' -N1 -D'+ name + '.math-N1.grd '
    print (command)
    os.system(command)

# los vector angles
# first band is incidence in degrees
los1 = ts.load_mmap(los, width, length, quiet=True, map='BIL', nchannels=2, channel=1, conv=False)
#los1 = los1[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]
#los1[los1 == 0]= np.nan
los1.tofile('los_inc.flt')
command = 'xyz2grd ' + 'los_inc.flt ' + '-R' + str(startLon) + '/' + str(endLon) + '/' + str(endLat) + '/' + str(startLat) + ' -I' + str(width) + '+/' + str(length) + '+ -ZTLf -N0.0 -G' + 'los_inc.grd -V'
print (command)
os.system(command)

# second band is LOS direction in "math" convention degrees
los2 = ts.load_mmap(los, width, length, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)
#los2 = los2[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]
# convert from math convention to geophysics convention for azimuth
azim=180.0-los2  
azim.tofile('los_azim.flt')
head=azim-90.0
head.tofile('los_head.flt')
command = 'xyz2grd ' + 'los_head.flt ' + '-R' + str(startLon) + '/' + str(endLon) + '/' + str(endLat) + '/' + str(startLat) + ' -I' + str(width) + '+/' + str(length) + '+ -ZTLf -N90.0 -G' + 'los_head.grd -V'
print (command)
os.system(command)

# now dealing with the DEM 
Dfile = 'hgt.flt'
    
hdata = np.memmap(DEMfile, dtype=np.int16,shape=(length, width), mode='r')
fdata = np.memmap(Dfile, dtype=np.float32,shape=(length,width), mode='w+')

for kk in xrange(length):
            fdata[kk,:] = hdata[kk,:]



#DEM GMT command
command = 'xyz2grd ' + Dfile + ' ' + '-R' + str(startLon) + '/' + str(endLon) + '/' + str(endLat) + '/' + str(startLat) + ' -I' + str(width) + '+/' + str(length) + '+ -ZTLf -N0.0 -Gdem.grd -V'
print (command)
os.system(command)

# now masked phase
if os.path.exists(nameM):
    phsM = ts.load_mmap(nameM, width, length, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)
    phsM = np.array(phsM)

    # masked phase does not need coherence mask, just convert masked area to NaNs
    phsM[phsM == 0]= np.nan
    phsM.tofile(nameM + '.phs')

    command = 'xyz2grd ' + nameM + '.phs ' + '-R' + str(startLon) + '/' + str(endLon) + '/' + str(endLat) + '/' + str(startLat) + ' -I' + str(width) + '+/' + str(length) + '+ -ZTLf -N0.0 -G' + nameM + '.grd -V'
    print (command)
    os.system(command)
    command = 'grdmath ' + nameM + '.grd ' + str((wavelength/(4 *Pi))) + ' MUL = ' + nameM + '.math.grd'
    print (command)
    os.system(command)
    command = 'grdtrend ' + nameM + '.math.grd' + ' -N1 -D'+ nameM + '.math-N1.grd '
    print (command)
    os.system(command)





# opening the unwrapped interferogram
## Old way for opening and interferogram

##a=np.fromfile(name, np.complex64)
##amp=np.zeros((length*width))
##phs=np.zeros((length*width))

##for i in range(length*width): # you extract both the modulus and the phase
    #amp[i]=abs(a[i])
##   	phs[i]=c.phase(a[i])

##nz = np.where(np.isnan(phs))
##phs[nz] = 0.0
##phs = phs.astype(np.float32)
##phs.tofile(name + '.phs')
