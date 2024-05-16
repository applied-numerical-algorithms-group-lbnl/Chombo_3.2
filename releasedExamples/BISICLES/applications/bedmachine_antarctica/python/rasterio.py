
from osgeo import gdal, osr
from osgeo.gdalconst import *

def write_raster(data,path, EPSG=3031):
    """
    write BisiclesData speeed,Tb to a geotiff path
    """
    N_BANDS = 2
    x = data.x
    y = data.y
    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    dx = x[1]-x[0]
    dataset = driver.Create(path, len(x), len(y), N_BANDS, gdal.GDT_Float64)
    dataset.SetGeoTransform( ( np.min(x), dx, 0, np.max(y), 0, -dx) ) 
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    wkt_projection = srs.ExportToWkt()
    dataset.SetProjection(wkt_projection)
    np.ma.set_fill_value(data,-9999.0)
    
    #speed
    d = np.flipud((data.speed)/365.0)
    dataset.GetRasterBand(1).WriteArray(d)
    dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
    
    
    #basal traction
    d = np.flipud(data.Tb/1.0e3)
    dataset.GetRasterBand(2).WriteArray(d)
    dataset.GetRasterBand(2).SetNoDataValue(-9999.0) 

    dataset.FlushCache()  # Write to disk. 

def read_raster(raster_path):
    """        
    Opens a tiff as specified by the user    
    Returns an array of the raster with co-oordinates
    """
    from osgeo import gdal,gdalconst  
    import numpy as np
    
    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    src = gdal.Open(raster_path, gdalconst.GA_ReadOnly)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)

    
    data=src.ReadAsArray()
    print("Opened %s" %(raster_path))
    #print(src.GetMetadata())
    
    tol = 1.0e-6
    x = np.arange(ulx,lrx-tol, +xres)
    y = np.arange(lry,uly-tol, -yres)
    
    return x,y,np.flipud(data[:,:])*365 # m/day -> m/a
