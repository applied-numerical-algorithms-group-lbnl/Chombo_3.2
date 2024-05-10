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

parser.add_argument('--data_directory', type=str, help='directory where data lives' ,default="_visc_compare_5_10_2024")
parser.add_argument('--output_prefix', type=str, help='Prefix for output file[summary_of_]',default="summary_of")
parser.add_argument('--max_num_proc', type=int, help='max number of processors for each run'   ,default='8')
parser.add_argument('--max_num_proc_2d', type=int, help='max number of processors for 2d runs'   ,default='4')


args = parser.parse_args()
print(args)
home_str = os.getcwd();
print ("homedir = " + home_str)
today_str =str(today.month) + "_" + str(today.day) + "_" + str(today.year)
summary_file_name  = args.output_prefix + "_" + args.data_directory + ".tex"
print ("summary_file_name = " + summary_file_name)

top_directory = home_str + "/" + args.data_directory
print("top_directory=" + top_directory)
label_str = "tab::data_reduction_table_" + today_str
caption_str = "Performance data for " + args.data_directory
f_summary       = open( summary_file_name,'w')
f_summary.write("\\begin{small} \n ")
f_summary.write("\\begin{table} \n ")
f_summary.write("\\begin{center}\n ")
f_summary.write("\\begin{tabular}{|c|c|c|c|c|c||c|} \\hline \n")
f_summary.write("Op & D & Case & $N_p$ & Final $|R|$  &  Iter & Main Time \\\\  \n")
opt_status = "opt_high"
deb_status = "debug_false"
mpi_status  = "mpi_true"
main_str  = "main_time"
resi_str = "norm_res"

i_dim = 2
while i_dim <= 3:
    dim_status = "dim_"  + str(i_dim)
    i_max_iter = 0
    i_opera = 0
    #while i_opera <= 0:
    while i_opera <= 2:
        opera_str      = "4586"
        subdir_name    = "4586"
        operator_entry = "4586"
        pout_base_name      = "4586"
        if(i_opera == 0):
            opera_str   = "chf_viscous_tensor_op"
            subdir_name = "_old_viscous_tensor"
            operator_entry = "ChF Viscous Op"
            pout_base_name      = "pout_old_vto.0"
        if(i_opera == 1):
            opera_str       = "proto_viscous_tensor_op"
            subdir_name     =   "amr_viscous_tensor"
            operator_entry  = "Proto Slow VTO"
            pout_base_name  = "pout.0"
        if(i_opera == 2):
            opera_str       = "proto_fast_vto"
            subdir_name     =   "amr_fast_vto"
            operator_entry  = "Proto Fast VTO"
            pout_base_name  = "pout_fastvto.0"

        operator_name = opera_str
        #exit()
        #i_max_iter = i_max_iter + 1
        #if(i_max_iter > 1000):
        #    exit()

        i_min_case = 0
        i_max_case = 4
        i_case = i_min_case
        while i_case <= i_max_case:
            case_status= "case_"  + str(i_case)
            str_config = opera_str + "." + case_status + "." + dim_status+ "."  + opt_status+ "."  +  deb_status 
            print_str = "configuration string = " + str_config
            print (print_str)

            config_directory = top_directory + "/" + str_config
            print_str =  "configuration directory = " + config_directory
            print(print_str)

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
                    print_str = "directory " + rundir_name + " not found"
                    print (print_str)
                    exit()
                pout_name = rundir_name + "/" + pout_base_name
                time_name = rundir_name + "/time.table.0"
                comm_str_resi = "grep Final_residual  " + pout_name
                comm_str_iter  = "grep Final_iteration " + pout_name
                comm_str_time  = "grep Total_Time "      + time_name
                has_pout =  os.path.exists(pout_name) 
                has_time =  os.path.exists(time_name) 
                if( has_pout ):
                    print_str = "YES pout file  found in "  + config_directory
                    print(print_str)            
                else:                           
                    printstr = pout_name
                    print_str = "NO  pout file  found in "  + config_directory
                    print(print_str)            
                if( has_time ):                 
                    print_str = "YES time file  found in "  + config_directory
                    print(print_str)            
                else:                           
                    print_str = "NO  time file  found in "  + config_directory
                    print(print_str)
                if (has_time and has_pout):
                    print_str = "YES both files found in "  + config_directory
                    print(print_str)

                        
                    output_str_resi  = "4586"
                    output_str_iter  = "4586"
                    output_str_time  = "4586"
                    bork_flag = 0
                    try:
                        output_str_resi  = subprocess.check_output(comm_str_resi , shell=True)
                    except subprocess.CalledProcessError as err_resi:
                        bork_flag = 1
                        print(err_resi.output)
                            
                    try:
                        output_str_iter  = subprocess.check_output(comm_str_iter , shell=True)
                    except subprocess.CalledProcessError as err_iter:
                        bork_flag = 1
                        print (err_iter.output)

                    try:
                        output_str_time  = subprocess.check_output(comm_str_time , shell=True)
                    except subprocess.CalledProcessError as err_time:
                        bork_flag = 1
                        print (err_time.output)

                    if(bork_flag == 0):
                        resi_list =  output_str_resi.split()
                        iter_list =  output_str_iter.split()
                        time_list =  output_str_time.split()
                        resi_len  = len(resi_list)
                        iter_len  = len(iter_list)
                        time_len  = len(time_list)
                        resi_field = resi_list[resi_len-1]
                        iter_field = iter_list[iter_len-1]
                        time_field = time_list[time_len-1]
                        file_entry = operator_entry + " & " + str(i_dim) + " & " + str(i_case) + " & " + str(i_num_proc) + "& " + resi_field + " & " + iter_field  + " & " + time_field + "\\\\";
                        f_summary.write(file_entry + "\n");
                    else:
                        print_str = "borkflag causes this case to be skipped "
                            
                else:
                    print_str = "skipping this case for lack of data"
                    print( print_str)
                print("end inner loop")
                i_num_proc = 2*i_num_proc
                    
            print_str = "end loop over procs"
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
f_summary.write("\\label{" +     label_str + "} \n") 
f_summary.write("\\caption{" + caption_str + "} \n" )
f_summary.write("\\end{table} \n")
f_summary.write("\\end{small} \n")

print("closing summary and exiting")
f_summary.close()
        
