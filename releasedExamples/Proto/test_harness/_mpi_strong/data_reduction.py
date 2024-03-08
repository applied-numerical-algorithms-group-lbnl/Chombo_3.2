#!/usr/bin/python
from argparse import *
import os
import glob
import platform
from datetime import date
import signal
import subprocess
import re

def signal_handler(signum, frame):
    raise Exception("Timed out!")
#./_archive_no_mains_no_hdf/_prch_strong_3_5_2024/old_helmholtz_op.case_3.dim_2.opt_high.debug_false/4_procs/pout.0: RMultiGrid::solveNoInitResid:   Final residual norm  = 3.135270e-13, Multigrid iterations = 13
today = date.today()

print("Today's date:", today)

parser = ArgumentParser()

parser.add_argument('--data_directory', type=str, help='Location of prch_strong directory' ,default="./_archive_no_mains_no_hdf/_prch_strong_3_5_2024/")
parser.add_argument('--output_prefix', type=str, help='Prefix for output file',default="summary")
parser.add_argument('--max_num_proc', type=int, help='max number of processors for each run'   ,default='8')
parser.add_argument('--max_num_proc_2d', type=int, help='max number of processors for 2d runs'   ,default='4')

args = parser.parse_args()
print(args)
homestr = os.getcwd();
print ("homedir = " + homestr)
strtoday =str(today.month) + "_" + str(today.day) + "_" + str(today.year)
summary_file_name  = args.output_prefix + "_" + strtoday + ".tex"
print ("summary_file_name = " + summary_file_name)

f_summary       = open( summary_file_name,'w')
f_summary.write("\\begin{small} \n ")
f_summary.write("\\begin{table} \n ")
f_summary.write("\\begin{center}\n ")
f_summary.write("\\begin{tabular}{|c|c|c|c|c|c||c|} \\hline \n")
f_summary.write("Op & D & Case & $N_p$ & Final $|R|$  &  Iter & Main Time \\\\  \n")
top_dir  = args.data_directory 
case_label     = "case_label"
case_caption   = "case_caption"
i_dim = 2
while i_dim < 3:
    i_opera = 0
    #make this <= 3 to turn on resistivity
    while i_opera < 1:
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
        i_min_case = 0
        i_max_case = 3
        #  set to max_case to make a short test
        #i_case = i_min_case
        i_case = i_max_case
        while i_case <= i_max_case:
            i_old_amr = 0
            while i_old_amr <= 1:
                operator_name = "amr_" + op_str
                op_entry = "op_str + (Proto)"

                if(i_old_amr == 1):
                    operator_name = "old_" + op_str
                    op_entry = op_str +  "(ChF)"

                    i_max_proc = args.max_num_proc
                    if(i_dim == 2):
                        i_max_proc = args.max_num_proc_2d
                    i_num_proc = 1
                    while i_num_proc <= i_max_proc:
                        main_str  = "main_time"
                        resid_str = "norm_res"

                        comm_str_resid = "grep Final_residual  " + top_dir + "*/*/pout.0        | grep "  + operator_name + " | grep dim_" + str(i_dim) + " | grep case_" + str(i_case) + " | grep " + str(i_num_proc) + "_procs"
                        comm_str_iter  = "grep Final_iteration " + top_dir + "*/*/pout.0        | grep "  + operator_name + " | grep dim_" + str(i_dim) + " | grep case_" + str(i_case) + " | grep " + str(i_num_proc) + "_procs"
                        comm_str_time  = "grep Total_Time "      + top_dir + "*/*/time.table.0  | grep "  + operator_name + " | grep dim_" + str(i_dim) + " | grep case_" + str(i_case) + " | grep " + str(i_num_proc) + "_procs"
                        comm_str = comm_str_time
                        comm_str_bash = "set -o  pipefail;" + comm_str
                        output = subprocess.check_output(comm_str, shell=True)
#                        if ("Final" not in output):
                        if ("main" not in output):
                            print_str("grep found nothing here")
                            print(print_str)
                        else:
                            print_str = "Op=" + op_str +  ", case = " + str(i_case) +  ", Np = " + str(i_num_proc)
                            print(print_str )
                            print_str = "full command = " + comm_str
                            print(print_str)
                            print_str = "output       = " + output
                            print(print_str)



#                        #print("output = " + output)
#                        file_entry = op_entry + " & " + str(i_dim) + " & " + str(i_case) + " & " + str(i_num_proc) + "& " + resid_str + " & " + "max_iter"  + " & " + main_str + "\\\\";
#                        f_summary.write(file_entry + "\n");
                        i_num_proc = 2*i_num_proc
                       
                print_str = "end loop over procs"
                print( print_str )
                i_old_amr = i_old_amr +1

            print_str  = "end loop over i_old_amr "
            f_summary.write("\\hline \n")
            print( print_str )
            i_case = i_case + 1
        print_str  = "end loop over cases "
        print( print_str )
        i_opera = i_opera +1

    print_str  = "end loop over operators " 
    print(print_str)
    i_dim = i_dim+ 1

print_str  = "end loop over dimensions " 
print(print_str)


f_summary.write("\\end{tabular} \n")
f_summary.write("\\end{center}   \n")
f_summary.write("\\label{"   + case_label   + "} \n") 
f_summary.write("\\caption{" + case_caption + "} \n" )
f_summary.write("\\end{table} \n")
f_summary.write("\\end{small} \n")

print("closing summary and exiting")
f_summary.close()
#########i_dim = 2
########## loop through dimensions
##########should be <=
##########  set to < for short test
##########while i_dim <= 2:
#########while i_dim <= 3:
#########    opt_status = "opt_high"
#########    opt_gmake  = " OPT=HIGH "
#########    deb_status = "debug_false"
#########    deb_gmake  = " DEBUG=FALSE "
#########    dim_status = "dim_"  + str(i_dim)
#########    dim_gmake  = " DIM=" + str(i_dim) + " "
#########    mpi_gmake  = " MPI=TRUE "
#########    mpi_status  = "mpi_true"
#########    str_compile = deb_gmake + dim_gmake + opt_gmake + mpi_gmake
#########    i_max_iter = 0
#########    i_opera = 0
#########    #while i_opera <= 0:
#########    while i_opera <= 3:
#########        op_str = "4586"
#########        if(i_opera == 0):
#########            op_str = "helmholtz"
#########        if(i_opera == 1):
#########            op_str = "conductivity"
#########        if(i_opera == 2):
#########            op_str = "viscous_tensor"
#########        if(i_opera == 3):
#########            op_str = "resistivity"
#########
#########        #loop through old/amr
#########        i_old_amr = 0
#########        while i_old_amr <= 1:
#########                
#########            exec_prefix = "old"
#########            exec_directory   = homestr + "/../../_old_" + op_str
#########            if (i_old_amr == 1):
#########                exec_prefix = "amr"
#########                exec_directory   = homestr + "/../../amr_" + op_str 
#########            
#########            opera_status =  exec_prefix + "_" + op_str + "_op" 
#########            option_string =  op_str + "." + dim_status + "."  + deb_status + "." + mpi_status
#########            executable_name = exec_prefix + "." + option_string + ".exe"
#########            full_exec_name = all_exec_dir_name + "/" + executable_name 
#########            command_str = "cd " + exec_directory + "; make main " + str_compile + "; mv main.exe " + full_exec_name + ";"
#########            print_str  = "about to compile  " + full_exec_name
#########            print(print_str)
#########            os.system(command_str)
#########            print_str  = "compiled/archived " + full_exec_name
#########            print(print_str)
#########            #i_max_iter = i_max_iter + 1
#########            #if(i_max_iter > 1000):
#########            #    exit()
#########
#########            i_min_case = 0
#########            i_max_case = 3
#########            #  set to max_case to make a short test
#########            i_case = i_min_case
#########            #i_case = i_max_case
#########            while i_case <= i_max_case:
#########                case_status= "case_"  + str(i_case)
#########                str_config = opera_status + "." + case_status + "." + dim_status+ "."  + opt_status+ "."  +  deb_status 
#########                print_str = "configuration string = " + str_config
#########                print (print_str)
#########
#########                config_directory = top_directory + "/" + str_config
#########                print_str =  "configuration directory = " + config_directory
#########                print(print_str)
#########                if not os.path.exists(config_directory):
#########                    print_str = "making directory " + config_directory
#########                    print (print_str)
#########                    os.mkdir(config_directory)
#########
#########                #(inner loop)
#########                # loop through num processors
#########                i_max_proc = args.max_num_proc
#########                if(i_dim == 2):
#########                    i_max_proc = args.max_num_proc_2d
#########                i_num_proc = 1
#########                while i_num_proc <= i_max_proc:
#########                    proc_status = str(i_num_proc) + "_procs"
#########                    rundir_name = config_directory + "/" + proc_status
#########                    print_str = " rundir_name = " + rundir_name
#########                    print (print_str)
#########                    if not os.path.exists(rundir_name):
#########                        print_str = "making directory " + rundir_name
#########                        print (print_str)
#########                        os.mkdir(rundir_name)
#########                    #soft link executable to rundir_name/main.exe
#########                    command_str = "ln -s "   + full_exec_name + " "     + rundir_name + "/main.exe"
#########                    print_str   = "linking " + full_exec_name + " to "  + rundir_name + "/main.exe"
#########                    print(print_str)
#########                    os.system(command_str)
#########                    # copy the right copy_n.inputs over
#########                    input_name = "case_" + str(i_case) + ".inputs"
#########                    command_str = "cd " + exec_directory + "; cp _inputs/" + input_name + " " + rundir_name + ";"
#########                    printstr  = "copying input file " + input_name + " to " + rundir_name
#########                    print(printstr)
#########                    os.system(command_str)
#########                    f_batchtemplate = open(batch_template,'r')
#########                    batch_file_name = rundir_name + "/" +  batch_root
#########                    printstr  = "creating batch file " + batch_file_name  + " from " + batch_template 
#########                    print(printstr)
#########                    f_batch = open(batch_file_name, 'w')
#########                    for batchster in f_batchtemplate:
#########                        t1str = batchster;
#########                        t2str = t1str.replace("NUM_NODE", str(i_num_proc))
#########                        t3str = t2str.replace("EXECUTABLE_FILE", "main.exe")
#########                        t4str = t3str.replace("INPUT_FILE", input_name)
#########                        f_batch.write(t4str)
#########
#########                    batch_command = "\n pushd " +  rundir_name + "; source " + batch_root + "; popd \n;"
#########                    f_run_all.write( batch_command)
#########                    printstr  = "closing batch file"
#########                    print(printstr)
#########                    f_batch.close()
#########                    f_batchtemplate.close()
#########                    i_num_proc = 2*i_num_proc
#########                    
#########                print_str = "end loop over procs"
#########                print( print_str )
#########                i_case = i_case + 1
#########            print_str  = "end loop over cases "
#########            print( print_str )
#########            i_old_amr = i_old_amr +1
#########
#########        print_str  = "end loop over i_old_amr "
#########        print( print_str )
#########        i_opera = i_opera +1
#########
#########    print_str  = "end loop over operators " 
#########    print(print_str)
#########    i_dim = i_dim+ 1
#########
#########print_str  = "end loop over dimensions "
#########print(print_str)
#########print_str = "Closing run all script file and exiting" 
#########print( print_str )
#########f_run_all.close()
        
