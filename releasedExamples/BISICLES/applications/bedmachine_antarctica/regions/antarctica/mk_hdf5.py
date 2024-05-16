INTERMEDIATE_DATA_PATH='../../intermediate_data'

def name(s,ext):
    return '{}/antarctica_bedmachine_{}.{}'.format(INTERMEDIATE_DATA_PATH,s,ext)

def temp_name(s,n_layer,ext,path='.'):
    return '{}/antarctica_bedmachine_temperature_morlighem_{}_{}.{}'.format(
        path,s,n_layer,ext)

def temp_m10_name(s,n_layer,ext,path='.'):
    return '{}/antarctica_bedmachine_temperature_morlighem_{}_{}_m10.{}'.format(
        path,s,n_layer,ext)


NCTOAMR='~/Development/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'

import os
def system(cmd):
    print(cmd)
    os.system(cmd)


system(NCTOAMR + ' {} antarctica_bedmachine_inverse_500m.2d.hdf5 umod umodc btrc'.format(name('500m','nc')))
system(NCTOAMR + ' {} antarctica_bedmachine_geometry_500m.2d.hdf5 thk topg'.format(name('500m','nc')))

args = ''
for i in range(0,24):
    args += 'temp{:06d} '.format(i)
    
system('{} {} {} {}'.format(NCTOAMR,
                               temp_name('4km','24','nc', path=INTERMEDIATE_DATA_PATH), 
                               temp_name('4km','24','2d.hdf5'),
                               args))

system('{} {} {} {}'.format(NCTOAMR,
                               temp_m10_name('4km','24','nc', path=INTERMEDIATE_DATA_PATH), 
                               temp_m10_name('4km','24','2d.hdf5'),
                               args))
