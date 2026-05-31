import argparse
import os 
import sys
import pandas as pd
import numpy as np

# Running together with the script:
# ================================== #

#     #!/bin/bash

#     #SBATCH --job-name=ggp
#     #SBATCH --cpus-per-task=1
#     #SBATCH --mem-per-cpu=4G

#     #SBATCH --time=23:00:00
#     #SBATCH --qos=1day

#     #SBATCH --output=std.out
#     #SBATCH --mail-type=END,FAIL,TIME_LIMIT

#     ${COMMAND}


def get_input_files(directory, keyword=None):
    if os.path.isfile(directory):
        return [directory]
    entries = os.listdir(directory)
    final_files = []
    if keyword == None:
        for e in entries:
            if e.endswith(".csv"):
                final_files.append(os.path.join(directory,e))
    else:
        for e in entries:
            if e.endswith(".csv") and keyword in e:
                final_files.append(os.path.join(directory,e))   
    return sorted(final_files)

         

def get_arg_list(args):
    s = ''
    for a in args:
        s += ' ' + a
    return s + ' '


def get_parameter_file(file, suffix):
    parameter_file = file[:-4] + suffix
    # directory = "/".join(file.split("/")[:-2]+ ["parameters"])
    # sample = "_".join(file.split("/")[-1].split("_")[:2])
    # parameter_file = os.path.join(directory, sample) + suffix
    return parameter_file

def get_col_from_csv_config(csv_config, col):
    with open(csv_config, 'r') as file:
        for line in file:
            if line.startswith(col):
                return line.split("=")[-1].strip()
    return None
    

########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Run RealTrace for all files in directory.",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-d',
                        dest='dir',
                        help='Directory with input files or path to input file',
                        required=True)
    
    parser.add_argument('-bin',
                        dest="bin",
                        help="Location of RealTrace binary",
                        default="./../../RealTrace/bin/RealTrace", 
                        required=False)
        
        
    parser.add_argument('-suffix',
                        dest="suffix",
                        help="Suffix of parameter file (for reused parameter files)",
                        default="", 
                        required=False)
    
    parser.add_argument('--dryrun', help="Shows what will be done", action='store_true')
    parser.add_argument('--cluster', help="Submit job to cluster (do not run directly)", action='store_true')
    
    
    parser.add_argument('-o',
                        dest='out' ,
                        help='Output directory (passed on to RealTrace)',
                        default=None, 
                        required=False)


    parser.add_argument('-space',
                        dest='space' ,
                        help='search space, log or linear (passed on to RealTrace)',
                        default='log', 
                        required=False)

    parser.add_argument('-t',
                        dest="tol",
                        help="Tolerance of maximization (passed on to RealTrace)",
                        default="1e-30", 
                        required=False)
    
    parser.add_argument('-m', help="Run maximization (passed on to RealTrace)", action='store_true')
    parser.add_argument('-p', help="Run prediction (passed on to RealTrace)", action='store_true')
    parser.add_argument('-j', help="Run joints (passed on to RealTrace)", action='store_true')

    args = parser.parse_args()

    # ======================================== #
    # ======================================== #


    if os.path.isfile(args.dir):
        print(args.dir, "is file")
        input_files = [args.dir]
        input_directory = os.path.join(*args.dir.split("/")[:-1])
    else:
        input_files = get_input_files(args.dir)
        print(len(input_files), "input files in", args.dir)
        input_directory = args.dir

    for infile in input_files:

        data = pd.read_csv(infile, skiprows=0, dtype=str)
        
        csv_config = input_directory + "/csv_config" + args.suffix +".txt"
        
        filter_col = get_col_from_csv_config(csv_config, "filter_col")
        if filter_col!=None:
            data = data[data[filter_col].astype(bool)]
            
        segment_col = get_col_from_csv_config(csv_config, "segment_col")
        if segment_col!=None:
            segment_idxs = sorted(np.unique(data[segment_col]))
            
        
        realtrace_arg = args.bin +\
                    " -c "      + csv_config + \
                    " -t "      + args.tol + \
                    " -space "  + args.space + \
                    " -i "      + infile 
        if args.m:
            realtrace_arg += ' -m '
        if args.p:
            realtrace_arg += ' -p '
        if args.j:
            realtrace_arg += ' -j '
            
        if segment_col!=None:
            parameter_file = " ".join([infile[:-4] + "_segment" + str(segment) + "_parameters" + args.suffix + ".txt" for segment in segment_idxs])
        else:
            parameter_file = infile[:-4] + "_parameters" + args.suffix + ".txt"

        realtrace_arg +=  " -b " + parameter_file

        if args.out != None:
            realtrace_arg +=  " -o " + args.out

        # ============ run! ============ #
        if args.cluster:
            com = "sbatch  --export=COMMAND='" + realtrace_arg +"'" + " submit_ggp_run.sl"
        else:
            com = realtrace_arg
       
        print(com)    
        if not args.dryrun:  
            os.system(com)


# ================================================================================ #
if __name__ == "__main__":
    main()