#!/usr/bin/python
from argparse import *
import os
import glob
import platform
from datetime import date

today = date.today()

print("Today's date:", today)

parser = ArgumentParser()

parser.add_argument('--input_template_dir', type=str, help='Directory of input file templates [../_input_templates].' ,default="../_input_templates")
parser.add_argument('--batch', type=str, help='batch file template'   ,default="../_batch_templates/spencer.batch")
parser.add_argument('--max_num_proc', type=int, help='max number of processors for each run'   ,default='8')
parser.add_argument('--prefix', type=str, help='name of test["prch_compare"]',default="prch_strong_scale")

args = parser.parse_args()
print(args)
homestr = os.getcwd();
print ("homedir = " + homestr)
neartop_directory = homestr + "/_" +args.prefix
strtoday =str(today.month) + "_" + str(today.day) + "_" + str(today.year)
top_directory = neartop_directory + "_" + strtoday
print ("top_directory = " + top_directory)
if not os.path.exists(top_directory):
    printstr = "making directory " + top_directory
    print (printstr)
    os.mkdir(top_directory)

batch_template = homestr + "/" + args.batch
configure_log_name = top_directory + "/configure_test.log"
run_all_file_name = top_directory + "/run_all.sh"
f_run_all       = open( run_all_file_name,'w')
f_configure_log = open(configure_log_name,'w')
f_run_all.write('#/usr/bin/csh\n')

i_min_case = 0
i_max_case = 3

batch_root = "batch_4586.sh"

# loop through cases to run
#should be min case
#  set to max_case to make a short test
i_case = i_min_case
while i_case <= i_max_case:
    i_dim = 2
    # loop through dimensions
    #should be <=
    #  set to < for short test
    while i_dim <= 3:

        #I hesitate to use the char '=' in directory names
        #This also allows the directory name to be less shouty
        opt_status = "_opt_high_"
        opt_gmake  = " OPT=HIGH "
        deb_status = "_debug_false_"
        deb_gmake  = " DEBUG=FALSE "
        dim_status = "_dim_"  + str(i_dim)
        dim_gmake  = " DIM=" + str(i_dim) + " "
        case_status= "_case_"  + str(i_case)
        case_gmake = "Unneeded though we could work it in."
        mpi_gmake  = " MPI=TRUE "
            
        str_compile = deb_gmake + dim_gmake + opt_gmake + mpi_gmake
        str_config = case_status + dim_status + opt_status +  deb_status 
        printstr = "configuration string = " + str_config
        print (printstr)
        f_configure_log.write(printstr + "\n")
        config_directory = top_directory + "/" + str_config
        printstr =  "configuration directory = " + config_directory
        print(printstr)
        f_configure_log.write(printstr + "\n")
        if not os.path.exists(config_directory):
            printstr = "making directory " + config_directory
            print (printstr)
            f_configure_log.write(printstr + "\n")
            os.mkdir(config_directory)

        i_opera = 0
        #set to <= if you want to run resistivity
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
                directory_prefix = "/_old_"
                exec_directory   = homestr + "/../../_old_" + op_str
                if (i_old_amr == 1):
                    directory_prefix = "/amr_"
                    exec_directory   = homestr + "/../../amr_" + op_str
                        
                #(inner loop)
                # loop through num processors
                i_max_proc = args.max_num_proc
                i_num_proc = 1
                while i_num_proc <= i_max_proc:
                    mpi_status = "_" + str(i_num_proc) + "_procs_"  
                    rundir_name = config_directory + directory_prefix + op_str + mpi_status
                    printstr = "rundir_name = " + rundir_name
                    print (printstr)
                    f_configure_log.write(printstr + "\n")
                    if not os.path.exists(rundir_name):
                        os.mkdir(rundir_name)
                            
                    # create executables and copy them over
                    # copy the right copy_n.inputs over
                    # make a batch file 
                    input_name = "case_" + str(i_case) + ".inputs"
                    command_str = "cd " + exec_directory + "; make main " + str_compile + "; mv main.exe " + rundir_name + "; cp _inputs/" + input_name + " " + rundir_name + ";"
                    printstr  = "command_str = " + command_str
                    print(printstr)
                    f_configure_log.write(printstr + "\n")

                    os.system(command_str)
                    f_batchtemplate = open(batch_template,'r')
                    batch_file_name = rundir_name + "/" +  batch_root
                    printstr = "batch_file_name = " + batch_file_name
                    print(printstr)
                    f_configure_log.write(printstr + "\n")
                    #exit()

                    f_batch = open(batch_file_name, 'w')
                    for batchster in f_batchtemplate:
                        t1str = batchster;
                        t2str = t1str.replace("NUM_NODE", str(i_num_proc))
                        t3str = t2str.replace("EXECUTABLE_FILE", "main.exe")
                        t4str = t3str.replace("INPUT_FILE", input_name)
                        f_batch.write(t4str)

                    batch_command = "\n pushd " +  rundir_name + "; source " + batch_root + "; popd \n;"
                    f_run_all.write( batch_command)
                    f_batch.close()
                    f_batchtemplate.close()
                    i_num_proc = 2*i_num_proc
                    # end loop over processor counts

                printstr =     "End of the loops over proc counts."
                print(printstr)
                f_configure_log.write(printstr + "\n")
                i_old_amr = i_old_amr + 1
                #loop through old/amr

            printstr =     "End of the old/amr iteration."
            print(printstr)
            f_configure_log.write(printstr + "\n")
            i_opera = i_opera + 1
            # end loop over operators


        printstr =     "End of operator list iteration."
        print(printstr)
        f_configure_log.write(printstr + "\n")
        i_dim = i_dim + 1
        # end loop over values of i_dim
            
    printstr =     "End of the loops over dimensions."
    print(printstr)
    f_configure_log.write(printstr + "\n")
    i_case = i_case + 1
    # end loop over values of i_case

printstr =     "End of the loop over cases."
print(printstr)
f_configure_log.write(printstr + "\n")
printstr =     "Closing run all."
print(printstr)
f_configure_log.write(printstr + "\n")
f_run_all.close()
printstr =     "Closing configuration log and leaving."
print(printstr)
f_configure_log.write(printstr + "\n")
f_configure_log.close();

