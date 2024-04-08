#!/usr/bin/python
from argparse import *
import os
import glob
import platform
from datetime import date
import signal

def signal_handler(signum, frame):
    raise Exception("Timed out!")


today = date.today()

print("Today's date:", today)

parser = ArgumentParser()

parser.add_argument('--input_template_dir', type=str, help='Directory of input file templates [../_input_templates].' ,default="../_input_templates")
parser.add_argument('--sbatch_instead_of_source', type=bool, help='Whether run.sh calls sbatch instead of source for each case[True].' ,default=True)
parser.add_argument('--batch', type=str, help='batch file template'   ,default="../_batch_templates/saul.host.batch")
parser.add_argument('--max_num_proc', type=int, help='max number of processors for each run'   ,default='8')
parser.add_argument('--max_num_proc_2d', type=int, help='max number of processors for 2d runs'   ,default='4')
parser.add_argument('--prefix', type=str, help='name of test["prch_compare"]',default="prch_saul_on_host")

args = parser.parse_args()
print(args)
homestr = os.getcwd();
print ("homedir = " + homestr)
use_sbatch = args.sbatch_instead_of_source
neartop_directory = homestr + "/_" +args.prefix
strtoday =str(today.month) + "_" + str(today.day) + "_" + str(today.year)
top_directory = neartop_directory + "_" + strtoday
print ("top_directory = " + top_directory)
if not os.path.exists(top_directory):
    print_str = "making directory " + top_directory
    print (print_str)
    os.mkdir(top_directory)

batch_template = homestr + "/" + args.batch
run_all_file_name = top_directory + "/run_all.sh"
f_run_all       = open( run_all_file_name,'w')
f_run_all.write('#/usr/bin/csh\n')
batch_root = "batch_4586.sh"
#make and save executables because redundant compilation is silly and slow
all_exec_dir_name = top_directory + "/executables"
if not os.path.exists(all_exec_dir_name):
    print_str = "making directory " + all_exec_dir_name
    print (print_str)
    os.mkdir(all_exec_dir_name)

print_str = "Creating and caching executables into " +  all_exec_dir_name
print( print_str )
i_dim = 2
# loop through dimensions
#should be <=
#  set to < for short test
#while i_dim <= 2:
while i_dim <= 3:
    opt_status = "opt_high"
    opt_gmake  = " OPT=HIGH "
    deb_status = "debug_true"
    deb_gmake  = " DEBUG=TRUE "
    dim_status = "dim_"  + str(i_dim)
    dim_gmake  = " DIM=" + str(i_dim) + " "
    mpi_gmake  = " MPI=TRUE "
    mpi_status  = "mpi_true"
    str_compile = deb_gmake + dim_gmake + opt_gmake + mpi_gmake
    i_max_iter = 0
    i_opera = 0
    #while i_opera <= 0:
    while i_opera <= 3:
        op_str = "4586"
        if(i_opera == 0):
            op_str = "helmholtz"
        if(i_opera == 1):
            op_str = "conductivity"
        if(i_opera == 2):
            op_str = "viscous_tensor"
        if(i_opera == 3):
            op_str = "resistivity"

        #loop through old/amr
        i_old_amr = 0
        while i_old_amr <= 1:
                
            exec_prefix = "old"
            exec_directory   = homestr + "/../../_old_" + op_str
            if (i_old_amr == 1):
                exec_prefix = "amr"
                exec_directory   = homestr + "/../../amr_" + op_str 
            
            opera_status =  exec_prefix + "_" + op_str + "_op" 
            option_string =  op_str + "." + dim_status + "."  + deb_status + "." + mpi_status
            executable_name = exec_prefix + "." + option_string + ".exe"
            full_exec_name = all_exec_dir_name + "/" + executable_name 
            command_str = "cd " + exec_directory + "; make main " + str_compile + "; mv main.exe " + full_exec_name + ";"
            print_str  = "about to compile  " + full_exec_name
            print(print_str)
            os.system(command_str)
            print_str  = "compiled/archived " + full_exec_name
            print(print_str)
            #exit()
            #i_max_iter = i_max_iter + 1
            #if(i_max_iter > 1000):
            #    exit()

            i_min_case = 0
            i_max_case = 3
            #  set to max_case to make a short test
            i_case = i_min_case
            #i_case = i_max_case
            while i_case <= i_max_case:
                case_status= "case_"  + str(i_case)
                str_config = opera_status + "." + case_status + "." + dim_status+ "."  + opt_status+ "."  +  deb_status 
                print_str = "configuration string = " + str_config
                print (print_str)

                config_directory = top_directory + "/" + str_config
                print_str =  "configuration directory = " + config_directory
                print(print_str)
                if not os.path.exists(config_directory):
                    print_str = "making directory " + config_directory
                    print (print_str)
                    os.mkdir(config_directory)

                #(inner loop)
                # loop through num processors
                i_max_proc = args.max_num_proc
                if(i_dim == 2):
                    i_max_proc = args.max_num_proc_2d
                i_num_proc = 1
                while i_num_proc <= i_max_proc:
                    proc_status = str(i_num_proc) + "_procs"
                    rundir_name = config_directory + "/" + proc_status
                    print_str = " rundir_name = " + rundir_name
                    print (print_str)
                    if not os.path.exists(rundir_name):
                        print_str = "making directory " + rundir_name
                        print (print_str)
                        os.mkdir(rundir_name)
                    #soft link executable to rundir_name/main.exe
                    command_str = "ln -s "   + full_exec_name + " "     + rundir_name + "/main.exe"
                    print_str   = "linking " + full_exec_name + " to "  + rundir_name + "/main.exe"
                    print(print_str)
                    os.system(command_str)
                    # copy the right copy_n.inputs over
                    input_name = "case_" + str(i_case) + ".inputs"
                    command_str = "cd " + exec_directory + "; cp _inputs/" + input_name + " " + rundir_name + ";"
                    printstr  = "copying input file " + input_name + " to " + rundir_name
                    print(printstr)
                    os.system(command_str)
                    f_batchtemplate = open(batch_template,'r')
                    batch_file_name = rundir_name + "/" +  batch_root
                    printstr  = "creating batch file " + batch_file_name  + " from " + batch_template 
                    print(printstr)
                    f_batch = open(batch_file_name, 'w')
                    for batchster in f_batchtemplate:
                        t1str = batchster;
                        t2str = t1str.replace("NUM_NODE", str(i_num_proc))
                        t3str = t2str.replace("EXECUTABLE_FILE", "main.exe")
                        t4str = t3str.replace("INPUT_FILE", input_name)
                        f_batch.write(t4str)

# srun directly. can only be done interactively
                    batch_command = "\n pushd " +  rundir_name + "; source " + batch_root + "; popd \n"
                    if(use_sbatch):
                        batch_command = "\n pushd " +  rundir_name + "; sbatch " + batch_root + "; popd \n"
                        
                    f_run_all.write( batch_command)
                    printstr  = "closing batch file"
                    print(printstr)
                    f_batch.close()
                    f_batchtemplate.close()
                    i_num_proc = 2*i_num_proc
                    
                print_str = "end loop over procs"
                print( print_str )
                i_case = i_case + 1
            print_str  = "end loop over cases "
            print( print_str )
            i_old_amr = i_old_amr +1

        print_str  = "end loop over i_old_amr "
        print( print_str )
        i_opera = i_opera +1

    print_str  = "end loop over operators " 
    print(print_str)
    i_dim = i_dim+ 1

print_str  = "end loop over dimensions "
print(print_str)
print_str = "Closing run all script file and exiting" 
print( print_str )
f_run_all.close()
        
