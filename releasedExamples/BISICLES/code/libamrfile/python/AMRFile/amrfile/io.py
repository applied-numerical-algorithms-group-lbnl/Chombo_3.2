#AMRFile hdf5 i/o functions via libamrfile (a C++ shared lib with C interface)

from ctypes import *;
import numpy as np
libamrfile = CDLL("libamrfile.so")

def __error__(status):
    if (status.value != 0):
        raise Exception('libamrfile error',status.value)
    

def load(filename):
    status = c_int(-1)
    amrID = c_int(-1)
    ascii_filename =  filename.encode('ascii')
    libamrfile.amr_read_file(pointer(status), pointer(amrID), ascii_filename)
    __error__(status)
    return amrID

def free(amrID):
    status = c_int(-1)
    libamrfile.amr_free(pointer(status), pointer(amrID)) 
    __error__(status)

def freeAll():
    status = c_int(0)
    libamrfile.amr_free_all(pointer((status)))
    __error__(status)

def queryCompID(amrID, name):
    #look up component ID by name
    status = c_int(-1)
    comp = c_int(-1)
    namelen = c_int(len(name))
    ascii_name =  name.encode('ascii')
    libamrfile.amr_query_comp_id(pointer(status), pointer(comp), pointer(amrID), ascii_name, pointer(namelen))
    __error__(status)
    return comp.value

def queryDomainCorners(amrID, level):
   status = c_int(-1)
   lo=np.intc([0,0])
   hi=np.intc([0,0])
   libamrfile.amr_query_domain_corners(pointer(status),
                                       lo.ctypes.data_as(POINTER(c_int)),
                                       hi.ctypes.data_as(POINTER(c_int)),
                                       pointer(amrID),
                                       pointer(c_int(level)))
   __error__(status)
   return lo,hi

def queryTime(amrID):
    status = c_int(-1)
    time = c_double(-1)
    libamrfile.amr_query_time(pointer(status),
                              pointer(time),
                              pointer(amrID))
    __error__(status)
    return time.value

def readBox2D(amrID, level, lo, hi, component, interpolationOrder = 0):

    compid = c_int(-1)
    if (isinstance(component,str)):
        compid = c_int(queryCompID(amrID, component))
    elif (isinstance(component,int)):
        compid = c_int(component)               
    

    status = c_int(-1)
    nx = hi[0] - lo[0]+1
    ny = hi[1] - lo[1]+1
    x = np.zeros(nx)
    y = np.zeros(ny)
    v = np.asfortranarray(np.zeros((nx,ny)))
    nplo = np.intc(lo)
    nphi = np.intc(hi)

    libamrfile.amr_read_box_data_2d(pointer(status), 
                                    v.ctypes.data_as(POINTER(c_double)), 
                                    x.ctypes.data_as(POINTER(c_double)),
                                    y.ctypes.data_as(POINTER(c_double)), 
                                    pointer(amrID), 
                                    pointer(c_int(level)), 
                                    nplo.ctypes.data_as(POINTER(c_int)), 
                                    nphi.ctypes.data_as(POINTER(c_int)), 
                                    pointer(compid), 
                                    pointer(c_int(interpolationOrder)))

    __error__(status)
    return x,y,v.transpose()

def queryLevelNumber(amrID):
    
    status = c_int(-1)
    n_level = c_int(-1)
    libamrfile.amr_query_n_level(pointer(status), pointer(n_level), pointer(amrID))
    
    return n_level.value
    
def queryFABNumber(amrID, level):
    
    status = c_int(-1)
    n_fab = c_int(-1)
    libamrfile.amr_query_n_fab(pointer(status), pointer(n_fab),
                               pointer(amrID), pointer(c_int(level)))
    
    return n_fab.value
    
def readFAB(amrID, level, fab, comp, ng=0):
    

    
    status = c_int(-1)
    nx = c_int(-1)
    ny = c_int(-1)
    nz = c_int(-1)
    
    libamrfile.amr_query_fab_dimensions_2d(pointer(status), pointer(nx), pointer(ny),
                                           pointer(nz), pointer(amrID),
                                           pointer(c_int(level)), 
                                           pointer(c_int(fab)))
     
    nx = nx.value + 2*ng
    ny = ny.value  + 2*ng                                      
    x = np.zeros(nx)
    y = np.zeros(ny)
    v = np.asfortranarray(np.zeros((nx,ny)))

    libamrfile.amr_read_fab_data_2d(
            pointer(status),
            v.ctypes.data_as(POINTER(c_double)), 
            x.ctypes.data_as(POINTER(c_double)),
            y.ctypes.data_as(POINTER(c_double)), 
            pointer(amrID),        
            pointer(c_int(level)), 
            pointer(c_int(fab)),
            pointer(c_int(comp)),
            pointer(c_int(ng)))
    
    return x,y,v.transpose()