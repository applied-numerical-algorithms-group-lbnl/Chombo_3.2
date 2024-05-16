

EXTERNAL_DATA_PATH='../external_data'
#BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineGreenland_2019-11-05_v01.nc')
BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineGreenland-2021-04-20.nc')
MEASURES_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'Greenland_ice_speed_v2017.nc')
INTERMEDIATE_DATA_PATH='.'

def name(s):
    return '{}/greenland_bedmachine_{}.nc'.format(INTERMEDIATE_DATA_PATH,s)

OUTPUT_NC = name('150m')
N_LAYER=24

import os, time

def create_new_only(file_name, fun):
    if os.path.isfile(file_name):
        mb = int(os.path.getsize(file_name)/1024/1024)
        mt = time.ctime(os.path.getsize(file_name))
        print ('{} ({} Mb {}) exists - delete if you want to recreate it'.format(file_name, mb, mt ))
    else:
        print ('creating {} ...'.format(file_name))
        fun(file_name)
        print ('...done')

def preprocess_150m(output_nc):        
    from preprocess_thk_bed_btrc import preprocess
    preprocess(output_nc, BEDMACHINE_NC, MEASURES_NC)
    return None

def preprocess_therm(output_nc):
    from preprocess_therm_bc import preprocess    
    preprocess(output_nc)
    return None
    
    
                        
create_new_only(OUTPUT_NC, preprocess_150m)

#create coarse grid versions for later convenience
suffix = ['150m','300m','600m','1200m','2400m','4800m']

for i in range(1,len(suffix)):
    fine = name(suffix[i-1])
    coarse = name(suffix[i])
    def f(file_name):
        from coarsen import coarsenc
        coarsenc(file_name , fine)
        return(None)
    create_new_only(coarse, f)

create_new_only('greenland_bedmachine_therm_bc_4800m.nc', preprocess_therm)