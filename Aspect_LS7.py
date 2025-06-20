#!/usr/bin/env/ python3

#============================================================
# Modules
#============================================================
# import ee & geetools
import ee

# other modules
from pylab import *
import numpy as np
import time
import math
import datetime
import matplotlib.pyplot as plt
import sys

# Initialization of GEE
ee.Initialize()

time1 = time.time()

#============================================================
# Difinitions
#============================================================
RADIX = 2
#ndsithresh = 0.45

# Masking Cloud, Shadow, Scan line
def extractQABits(qaBand, bitStart, bitEnd):
    numBits = bitEnd - bitStart +1
    qaBits = qaBand.rightShift(bitStart).mod(math.pow(RADIX, numBits))
    return qaBits

def CloudMask_LS(image):
    qa = image.select('BQA')
    mask = qa.bitwiseAnd(1 << 4).eq(0)
    return mask

def ShadowMask_LS(image):
    qa = image.select('BQA')
    bitStartShadowConf = 7
    bitStartShadowConf = 8
    qaBitsShadowConf = extractQABits(qa,bitStartShadowConf,bitStartShadowConf)
    mask = qaBitsShadowConf.lt(1)
    return mask

def ScanLineMask(image):
    area = image.select('B1')
    area = area.multiply(0).add(1)
    #scan line detect
    edge = ee.Algorithms.CannyEdgeDetector(area,1,1)
    connected = edge.mask(edge).lt(0.8).connectedPixelCount(35, True)
    edgeLong = connected.gte(35)
    edge = edgeLong
    line = area.mask(edge.focal_max(30, 'square', 'meters')).rename('ScanLine')
    #mask
    mask = area.where(line,0)
    return mask


def applyMask(image):
	Cmask = CloudMask_LS(image)
	Smask = ShadowMask_LS(image)
	Lmask = ScanLineMask(image)
	masked = image.updateMask(Cmask)
	masked = masked.updateMask(Smask)
	masked = masked.updateMask(Lmask)
	masked = masked.updateMask(localRGImask) # glaciers mask
	masked = masked.updateMask(WBmask) # water bodies mask
	return masked


# NDSI
def NDSI(image):
    ndsi = image.normalizedDifference(['B2','B5']).rename('ndsi')
    return image.addBands([ndsi])


# Snow area, Snow line, DEM
def SnowAreaDetect(image):
    ndsi = image.select('ndsi')
    snow = ndsi.gte(ndsithresh)
    cls = ndsi.multiply(0)
    cls = cls.where(snow,1).rename('SnowArea')
    return image.addBands([cls])

def EdgeDetect(image):
    snowarea = image.select('SnowArea')
    edge = ee.Algorithms.CannyEdgeDetector(snowarea,1,1)
    connected = edge.mask(edge).lt(0.8).connectedPixelCount(35, True)
    edgeLong = connected.gte(35)
    edge = edgeLong
    buffer = snowarea.mask(edge.focal_max(30, 'square', 'meters')).rename('SnowLine')
    return image.addBands([buffer])


# get Altitude on Snow Line
def MaskBySnowLine(image):
    dem_line = DEM.select('AVE_DSM')
    dem_line = dem_line.updateMask(image.select('SnowLine'))
    return image.addBands([dem_line])


# get average alt & date
def getAveAltAsp(image):
    #calculate ave alt for 4 direction
    asp1 = image.updateMask(mask4.select('North'))
    asp2 = image.updateMask(mask4.select('East'))
    asp3 = image.updateMask(mask4.select('South'))
    asp4 = image.updateMask(mask4.select('West'))
    val1 = asp1.reduceRegion(ee.Reducer.mean(),target,30).get('AVE_DSM')
    val2 = asp2.reduceRegion(ee.Reducer.mean(),target,30).get('AVE_DSM')
    val3 = asp3.reduceRegion(ee.Reducer.mean(),target,30).get('AVE_DSM')
    val4 = asp4.reduceRegion(ee.Reducer.mean(),target,30).get('AVE_DSM')

    #others & output
    date = image.get('system:time_start')
    val = image.reduceRegion(ee.Reducer.mean(),target,30).get('AVE_DSM')
    num = image.reduceRegion(ee.Reducer.count(),target,30).get('SnowLine')
    ft = ee.Feature(None, {'system:time_start': date,
                        'date'  : ee.Date(date).format('YYYY-MM-dd'),
                        'value' : val,
                        'north' : val1,
                        'east'  : val2,
                        'south' : val3,
                        'west'  : val4,
                        'number': num})
    return ft


# get Julian Day
def datestdtojd (stddate):
    fmt='%Y-%m-%d'
    sdtdate = datetime.datetime.strptime(stddate, fmt)
    sdtdate = sdtdate.timetuple()
    jdate = sdtdate.tm_yday
    jyear = sdtdate.tm_year
    return(jdate, jyear)


#============================================================
# Output file names
#============================================================
# set glacier name & location
args = sys.argv
glacier = '_'+args[1]
REG = [np.float(args[3]), np.float(args[2])]

# set threshold of NDSI
ndsithresh = float(args[4])

# set period
ST = [1999, 1, 1]
ET = [2019, 12, 31]
#ET = [2000, 12, 31]
StartTime = datetime.datetime(ST[0], ST[1], ST[2])
EndTime = datetime.datetime(ET[0], ET[1], ET[2])

firstyear, lastyear = ST[0], ET[0]
nbyears = lastyear-firstyear + 1

# area = 'G' # to consider only the Glacier
area = 'C' # to consider the Catchment

print('ndsi threshold: '+str(ndsithresh))
print(glacier)


#============================================================
# Import satellite data
#============================================================
# set period
#StartTime = datetime.datetime(1999, 1, 1)
#EndTime = datetime.datetime(2013, 12, 31)

# get target area
sheds = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_9")
point = ee.Geometry.Point(REG)
#target = sheds.filter(ee.Filter.eq("HYBAS_ID",HID))
target = sheds.filterBounds(point)

# import Landsat 8
ls7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_TOA")\
        .filterBounds(target)\
        .filterDate(StartTime,EndTime)\
        .filter(ee.Filter.lte('CLOUD_COVER',50))

print('Number of Images (Landsat 7): ', ls7.size().getInfo())

#import DEM & reduce resolution
DEM = ee.Image("JAXA/ALOS/AW3D30/V2_2")
proj=DEM.projection().getInfo()
crs = proj['crs']
DEM100 = DEM.resample('bilinear').reproject(crs=crs,scale=100.0)

# Output file names
fOut = 'ASP'+glacier+'_LS7_'+area

#============================================================
# Main
#============================================================
# 1. Import Water Bodies from Global Surface Water
WBmask = ee.Image("JRC/GSW1_2/GlobalSurfaceWater")\
            .select("max_extent").unmask().focal_max(2).eq(0)

# 2. Import Randolph Glacier Inventory
rgi15 = ee.FeatureCollection("users/orie_sasaki/15_rgi60_SouthAsiaEast")
localRIG15 = rgi15.filterBounds(target)
rgi14 = ee.FeatureCollection("users/orie_sasaki/14_rgi60_SouthAsiaWest")
localRIG14 = rgi14.filterBounds(target)
rgi13 = ee.FeatureCollection("users/orie_sasaki/13_rgi60_CentralAsia")
localRIG13 = rgi13.filterBounds(target)

# 3. Import Scherler2018_global_debris
debris13 = ee.FeatureCollection("users/orie_sasaki/13_rgi60_CentralAsia_S2_DC_2015_2017_NDSI")
localDebris13 = debris13.filterBounds(target)
debris14 = ee.FeatureCollection("users/orie_sasaki/14_rgi60_SouthAsiaWest_S2_DC_2015_2017_NDSI")
localDebris14 = debris14.filterBounds(target)
debris15 = ee.FeatureCollection("users/orie_sasaki/15_rgi60_SouthAsiaEast_S2_DC_2015_2017_NDSI")
localDebris15 = debris15.filterBounds(target)

# 4. Reduce them to the target under study and merge the domains together
localDebris = localDebris13.merge(localDebris14).merge(localDebris15)
localRGI = localRIG15.merge(localRIG13).merge(localRIG14)

debrisUnmask = localDebris.reduceToImage(properties = ['Area'], reducer = ee.Reducer.first())
localRGImask = localRGI.reduceToImage(properties = ['Area'], reducer = ee.Reducer.first())

# 5. Unmask the debris covered part from the glaciers mask
if area == 'C':
    localRGImask = localRGImask.mask(debrisUnmask.unmask().focal_max(2).eq(0))
    localRGImask = localRGImask.unmask().eq(0)
else:
    localRGImask = localRGImask.mask().eq(1)


# Calculate slope aspect
aspect = ee.Terrain.aspect(DEM)
aspect100 = ee.Terrain.aspect(DEM100)

# Make direction mask
Names4 = ['Flat','North','East','South','West']
Thresh4 = [[45.0,135.0],[135.0,225.0],[225.0,315.0]]
Thresh4 = array(Thresh4)

zeros = aspect.multiply(0)
zeros100 = aspect100.multiply(0)
for num in range(0,5):
    if num == 0:
        #flat
        tmp1 = aspect100.lt(0.0)
        mask4 = zeros100.where(tmp1,1).rename('Flat')
    elif num == 1:
        #north
        tmp1 = aspect100.gte(0.0).And(aspect100.lte(45.0))
        tmp2 = aspect100.gt(315.0)
        mask = zeros100.where(tmp1,1)
        mask = mask.where(tmp2,1).rename('North')
        mask4 = mask4.addBands([mask])
    else:
        Asp = Names4[num]
        t1 = Thresh4[num-2,0]
        t2 = Thresh4[num-2,1]
        tmp1 = aspect100.gt(t1).And(aspect100.lte(t2))
        mask = zeros100.where(tmp1,1).rename(Asp)
        mask4 = mask4.addBands([mask])


# NDSI, Snow area, Edge detection
ls7_NDSI = ls7.map(NDSI)
SnowArea_ls7 = ls7_NDSI.map(SnowAreaDetect)
Edges = SnowArea_ls7.map(EdgeDetect)

# Masking cloud & shadow
Edges_masked = Edges.map(applyMask)

# Masking DEM by Snow Line
Edges_masked = Edges_masked.select('SnowLine')
dem_line = Edges_masked.map(MaskBySnowLine)

# get Average Altitude & date
AveAlt = dem_line.map(getAveAltAsp)


# Export to Google Drive
task = ee.batch.Export.table.toDrive(
    collection=AveAlt,
    description=fOut,
    folder='SnowLineAltitude',
    fileFormat='CSV')
task.start()

count=0
while task.active():
    print ('Time:',count*0.5,"min")
    #print (task.status())
    count = count+1
    time.sleep(30)

# last time
time2 = time.time()
TotalTime = (time2-time1) / 60.0
print (f"Total time: {TotalTime:.2f}min")
