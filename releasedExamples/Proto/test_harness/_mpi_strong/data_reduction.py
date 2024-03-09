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

parser.add_argument('--data_directory', type=str, help='Location of prch_strong directory' ,default="./_prch_strong_3_8_2024")
parser.add_argument('--output_prefix', type=str, help='Prefix for output file',default="summary")
parser.add_argument('--max_num_proc', type=int, help='max number of processors for each run'   ,default='8')
parser.add_argument('--max_num_proc_2d', type=int, help='max number of processors for 2d runs'   ,default='4')
parser.add_argument('--prefix', type=str, help='name of test["prch_compare"]',default="prch_strong")

args = parser.parse_args()
print(args)
home_str = os.getcwd();
print ("homedir = " + home_str)
today_str =str(today.month) + "_" + str(today.day) + "_" + str(today.year)
summary_file_name  = args.output_prefix + "_" + today_str + ".tex"
print ("summary_file_name = " + summary_file_name)
neartop_directory = home_str + "/_" +args.prefix
top_directory = args.data_directory

label_str = "tab::data_reduction_table_" + today_str
caption_str = "Performance data for " + args.data_directory
f_summary       = open( summary_file_name,'w')
f_summary.write("\\begin{small} \n ")
f_summary.write("\\begin{table} \n ")
f_summary.write("\\begin{center}\n ")
f_summary.write("\\begin{tabular}{|c|c|c|c|c|c||c|} \\hline \n")
f_summary.write("Op & D & Case & $N_p$ & Final $|R|$  &  Iter & Main Time \\\\  \n")
top_dir  = args.data_directory 

i_dim = 2
while i_dim < 3:
    dim_status = "dim_"  + str(i_dim)
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
            case_status= "case_"  + str(i_case)
            i_old_amr = 0
            while i_old_amr <= 1:
                operator_name = "old_" + op_str
                op_entry = op_str +  "(ChF)"
                exec_prefix = "old"
                if(i_old_amr == 1):
                    exec_prefix = "amr"
                    operator_name = "amr_" + op_str
                op_entry = op_str + "(Proto)"

                i_max_proc = args.max_num_proc
                if(i_dim == 2):
                    i_max_proc = args.max_num_proc_2d
                i_num_proc = 1
                while i_num_proc <= i_max_proc:
                    proc_status = str(i_num_proc) + "_procs"
                    main_str  = "main_time"
                    resid_str = "norm_res"
                    opt_status = "opt_high"
                    deb_status = "debug_false"

                    opera_status =  exec_prefix + "_" + op_str + "_op" 
                    config_str = opera_status + "." + case_status + "." + dim_status+ "."  + opt_status+ "."  +  deb_status 
                    config_directory = top_directory + "/" + config_str
                    rundir_name = config_directory + "/" + proc_status
                    pout_name = rundir_name + "/pout.0"
                    time_name = rundir_name + "/time.table.0"
                    comm_str_resid = "grep Final_residual  " + pout_name
                    comm_str_iter  = "grep Final_iteration " + pout_name
                    comm_str_time  = "grep Total_Time "      + time_name
                    has_pout =  os.path.exists(pout_name) 
                    has_time =  os.path.exists(time_name) 
                    if( has_pout ):
                        print_str = "YES pout file  found in "  + config_directory
                        print(print_str)            
                    else:                           
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
                        output_str_time  = subprocess.check_output(comm_str_time , shell=True)
                        output_str_resid = subprocess.check_output(comm_str_resid, shell=True)
                        output_str_iter  = subprocess.check_output(comm_str_iter , shell=True)
                        print_str_time  = "time_output  = " + output_str_time
                        print(print_str_time)
                        print_str_resid = "resid_output = "  + output_str_resid
                        print_str_iter  = "iter_output  = "  + output_str_iter
                        print(print_str_resid)
                        print(print_str_iter)
                        exit()
                    else:
                        print_str = "skipping this case for lack of data"
                        print( print_str)
                    i_num_proc = 2*i_num_proc
                    exit()
                    #end of inner loop
                       
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
f_summary.write("\\label{" +     label_str + "} \n") 
f_summary.write("\\caption{" + caption_str + "} \n" )
f_summary.write("\\end{table} \n")
f_summary.write("\\end{small} \n")

print("closing summary and exiting")
f_summary.close()
        
