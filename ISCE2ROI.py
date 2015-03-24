#! /usr/bin/env python
# Create ROI_pac format Inputs for Paul's version of Rowena's resamptool script
# modified from Pietro's script "PreparePaul.py" EJF 2014/12/11-29
# uses some GIAnT functions so must run under Python2
# does not yet convert the LOS angles to the ROI_pac convention

import numpy as np 
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import FormatStrFormatter
import sys
import h5py
import datetime as dt
import tsinsar as ts
import argparse
import os
import isce
from iscesys.Parsers.FileParserFactory import createFileParser
from iscesys.ImageApi import DataAccessor as DA
from isceobj.Util import key_of_same_content
import math
import cmath as c
from collections import OrderedDict
from xml.etree import ElementTree as ET

import lxml.objectify as OB

def findCo(startCo, delta, step):
	endCo = startCo + ( delta * step )
	return endCo

# set output file names--hardcoded for now
#outfile = 'geo_DSC_Coseismic_phs.BIL.unw'
outfile = 'geo_filt_phs.BIL.unw'
outfileCOH = 'geo_filt_coh.BIL.cor'
outfileLOS = 'geo_los.BIL.unw'
outfileAZO = 'geo_azo.BIL.unw'
outfileRGO = 'geo_rgo.BIL.unw'

name='filt_topophase.unw.geo'
nameM='masked.unw.geo' # assume same sizes
los = 'los.rdr.geo'

imagexml= name + '.xml'

#coh='topophase.cor.geo'
coh='phsig.cor.geo'
cohxml= coh + '.xml'
cohTH= 0.4   # higher threshold for filtered coherence

# pixel offset tracking output
offFile = 'offset_fields_offsets.bil.geo'
offFileXml = offFile + '.xml'
offSnrFile = 'offset_fields_snr.rdr.geo'
snrTH = 2 # snr threshold for offsets

# default wavelength (will try to get from insarProc.xml) 
# wavelength of CSK
wavelength = 0.0312283810417
# wavelength of TSX
#wavelength = 0.03106657823461874
#ampl='modeling/amp.flt'

# subsetting output--disabled for now
#COLTL = 976
#ROWTL = 1512
#COLBR = 1131
#ROWBR = 1696

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

# get metadata about image from interferogram .xml
ext = None
dataType = None
width = None
PA = createFileParser('xml')
dictNow, dictFact, dictMisc = PA.parse(imagexml)
numBands = 0

numBands = key_of_same_content('number_bands', dictNow)[1]
bandsType = key_of_same_content('data_type', dictNow)[1]
map = key_of_same_content('scheme', dictNow)[1]

coordinate1  =  key_of_same_content('coordinate1',dictNow)[1]
width = (int(key_of_same_content('size',coordinate1)[1]))
startLon=(float(key_of_same_content('startingValue',coordinate1)[1]))
deltaLon=(float(key_of_same_content('delta',coordinate1)[1]))
coordinate2  =  key_of_same_content('coordinate2',dictNow)[1]
length = (int(key_of_same_content('size',coordinate2)[1]))
startLat=(float(key_of_same_content('startingValue',coordinate2)[1]))
deltaLat=(float(key_of_same_content('delta',coordinate2)[1]))
print('WIDTH:  ' + str(width))
print('LENGTH: ' + str(length))

# include whole scene
ROWTL=0
COLTL=0
ROWBR=length
COLBR=width

LatTL=findCo(startLat, deltaLat, ROWTL)
LonTL=findCo(startLon, deltaLon, COLTL)
LatBR=findCo(startLat, deltaLon, ROWBR)
LonBR=findCo(startLon, deltaLat, COLBR)
width = COLBR - COLTL
length = ROWBR - ROWTL
print('subset WIDTH:  ' + str(width))
print('subset LENGTH: ' + str(length))

# get other metadata from insarProc.xml, if available (modified from prepGiant.py
if os.path.exists('insarProc.xml'):
     exampleXML = 'insarProc.xml'
     fid = open(exampleXML,'r')
     xObj = OB.fromstring(fid.read())
     fid.close()

     wavelength = float(xObj.master.wavelength)
     heading = float(xObj.runGeocode.inputs.PEG_HEADING) * 180.0 / np.pi	

print('radar wavelength ', wavelength)

# now get metadata about pixel offset tracking files, if present, in case different
if os.path.exists(offFileXml):
     fid = open(offFileXml,'r')
     xObjOff = OB.fromstring(fid.read())
     fid.close()

     PA2 = createFileParser('xml')
     dictNow2, dictFact2, dictMisc2 = PA2.parse(offFileXml)

     coordinate1  =  key_of_same_content('coordinate1',dictNow2)[1]
     widthOff = (int(key_of_same_content('size',coordinate1)[1]))
     startLonOff =(float(key_of_same_content('startingValue',coordinate1)[1]))
     deltaLonOff =(float(key_of_same_content('delta',coordinate1)[1]))
     coordinate2  =  key_of_same_content('coordinate2',dictNow2)[1]
     lengthOff = (int(key_of_same_content('size',coordinate2)[1]))
     startLatOff =(float(key_of_same_content('startingValue',coordinate2)[1]))
     deltaLatOff =(float(key_of_same_content('delta',coordinate2)[1]))
     print('Offsets width:  ' + str(widthOff))
     print('Offsets length: ' + str(lengthOff))


# End reading XML files 

# Importing requested data
# regular interferogram
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

    # output
    defo = phs	

# los vector angles
los1 = ts.load_mmap(los, width, length, quiet=True, map='BIL', nchannels=2, channel=1, conv=False)
los1 = los1[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]
los2 = ts.load_mmap(los, width, length, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)
los2 = los2[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]

# amplitude from interferogram
amp = ts.load_mmap(name, width, length, quiet=True, map='BIL', nchannels=2, channel=1, conv=False)
amp = amp[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]

#If we use time series we create a coherence map equals 1 except for nan values
#COH = ts.load_mmap(coh, width, length, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)
#COH = np.array(COH)
#COH = COH[min(ROWTL,ROWBR):max(ROWTL,ROWBR),min(COLTL,COLBR):max(COLTL,COLBR)]

#COH=np.ones((amp.shape[0],amp.shape[1]),dtype=np.float32)
#nz = np.where(np.isnan(defo))
#COH[nz] = 0.0

# now read the pixel offsets
if os.path.exists(offFile):
    azo = ts.load_mmap(offFile, widthOff, lengthOff, quiet=True, map='BIL', nchannels=2, channel=1, conv=False)
    rgo = ts.load_mmap(offFile, widthOff, lengthOff, quiet=True, map='BIL', nchannels=2, channel=2, conv=False)

    # first mask out pixels flagged by OffsetField
    azo = np.array(azo)
    azo[azo == -100.0] = np.nan
    rgo = np.array(rgo)
    rgo[rgo == -100.0] = np.nan

    # now mask with SNR
    SNR = ts.load_mmap(offSnrFile, widthOff, lengthOff, quiet=True, map='BIL', nchannels=1, channel=1, conv=False)
    SNR = np.array(SNR)
    # convert to mask
    SNR[SNR <= snrTH ] = 0
    SNR[SNR > snrTH ] = 1

    azom = SNR*azo
    azom[azom == 0]= np.nan
    rgom = SNR*rgo
    rgom[rgom == 0]= np.nan


#End Import requested data


#Start writing outputs

outunw=np.zeros((amp.shape[0],2*amp.shape[1]))
outunw[:,:width] = amp 
outunw[:,width:] = defo
outunw = outunw.astype('f')
outunw.tofile(outfile)

outunw=np.zeros((amp.shape[0],2*amp.shape[1]))
outunw[:,:width] = amp 
outunw[:,width:] = COH
outunw = outunw.astype('f')
outunw.tofile(outfileCOH )

outunw=np.zeros((amp.shape[0],2*amp.shape[1]))
outunw[:,:width] = los1
outunw[:,width:] = los2
outunw = outunw.astype('f')
outunw.tofile(outfileLOS)

if os.path.exists(offFile):
	if width==widthOff: # we can use amplitude from interferogram
		outunw=np.zeros((amp.shape[0],2*amp.shape[1]))
		outunw[:,:width] = amp 
		outunw[:,width:] = azom
		outunw = outunw.astype('f')
		outunw.tofile(outfileAZO)

		outunw=np.zeros((amp.shape[0],2*amp.shape[1]))
		outunw[:,:width] = amp 
		outunw[:,width:] = rgom
		outunw = outunw.astype('f')
		outunw.tofile(outfileRGO)

#create rsc file 

startLon=max(LonTL,LonBR)
endLon=min(LonTL,LonBR)

startLat=max(LatTL,LatBR)
endLat=min(LatTL,LatBR)

rdict = OrderedDict()

rdict['WIDTH'] = width
rdict['FILE_LENGTH'] = length
rdict['XMAX'] = length -1 
rdict['XMIN'] = 0
rdict['YMAX'] = width -1 
rdict['YMIN'] = 0
rdict['X_FIRST'] = LonTL
rdict['Y_FIRST'] = LatTL
rdict['Y_STEP'] = deltaLat
rdict['X_STEP'] = deltaLon
rdict['X_UNIT'] = 'degrees'
rdict['Y_UNIT'] = 'degrees'


rsccoh = outfileCOH + '.rsc'
ts.write_rsc(rdict, rsccoh)

rsclos = outfileLOS + '.rsc'
ts.write_rsc(rdict, rsclos)

rdict['WAVELENGTH'] = wavelength
rdict['HEADING'] = heading
rscname = outfile + '.rsc'
ts.write_rsc(rdict, rscname)

# pixel offsets already in meters, so set wavelength negative
if os.path.exists(offFile):
	rdict['WAVELENGTH'] = -1
	ts.write_rsc(rdict,outfileRGO+'.rsc')
        rdict['HEADING'] = heading + 90.0      # effective heading of azimuth offsets
	ts.write_rsc(rdict,outfileAZO+'.rsc')
