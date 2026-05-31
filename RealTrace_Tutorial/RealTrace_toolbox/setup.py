import argparse
import os 
import sys
import pandas as pd
import numpy as np

## CLASS FOR READING ##
class Raw_cell:
    def __init__(self, cell_id = 0, parent_id=-1):
        self.parent_id = parent_id
        self.cell_id = cell_id
        self.length = []
        self.log_length = []
        self.gfp = []
        self.time = []

def df2raw_cells(dataset, 
            time="time", 
            length="length", gfp="fp", 
            cell_id="cell_id", 
            parent_id="parent_id"):
    """ 
    dataset (pandas data frame as read from csv file)
    """
    cell_list = []
    last_cell = ""
    
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            if parent_id==None:
                new_cell = Raw_cell(cell_id=row[cell_id], parent_id=None)
            else:
                new_cell = Raw_cell(cell_id=row[cell_id], parent_id=row[parent_id])
            cell_list.append(new_cell)

        cell_list[-1].length.append(row[length])
        cell_list[-1].log_length.append(np.log(row[length]))
        cell_list[-1].gfp.append(row[gfp])
        cell_list[-1].time.append(row[time])

        last_cell = row[cell_id]
    return cell_list



## FILE MANAGEMENT ##

def get_input_files(directory, keyword=None):
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


def tags_str(tags):
    return ", ".join([tag for tag in tags])
        



def get_paramter_file(directory, infile, seg_idx):
    entries = os.listdir(directory)
    for e in entries:
        if infile.split("/")[-1][:-4] in e  and "_parameter_file.txt" in e and "segment" + str(seg_idx):
            return os.path.join(directory, e)



def estimate_mean_lambda(cells):

    growth_rate = []
    
    for idx, _ in enumerate(cells):
        gr = np.diff(cells[idx].log_length)/np.diff(cells[idx].time)

        growth_rate.append(gr)
    return np.mean(np.concatenate(growth_rate))

def estimate_mean_q(cells, beta=0):
    
    production_rate = []

    for idx, _ in enumerate(cells):
        log_length_interpl = np.exp(cells[idx].log_length[:-1] + np.diff(cells[idx].log_length)/2.) 
        gfp_interpl = cells[idx].gfp[:-1] + np.diff(cells[idx].gfp)/2.
        
        pr = np.diff(cells[idx].gfp)/np.diff(cells[idx].time)/log_length_interpl + gfp_interpl*beta/log_length_interpl
        production_rate.append(pr)
    return np.mean(np.concatenate(production_rate))


def estimate_var_x(cells, mean_lambda):
    deviations = []
    for idx, _ in enumerate(cells):
        deviations.append((cells[idx].log_length[:-1] + mean_lambda*np.diff(cells[idx].time)) - cells[idx].log_length[1:])
    return np.var(np.concatenate(deviations))

    
def estimate_var_g(cells, mean_q, beta):
    deviations = []
    for idx, _ in enumerate(cells):        
        abs_dev = cells[idx].gfp[:-1] + mean_q*np.diff(cells[idx].time)*cells[idx].log_length[:-1] \
                    - beta*np.diff(cells[idx].time)*cells[idx].gfp[:-1] \
                    - cells[idx].gfp[1:]
        deviations.append(abs_dev/np.sqrt(cells[idx].gfp[:-1]))
    return np.var(np.concatenate(deviations))


def estimate_parameters(data_df, args, beta):
    cells_data = df2raw_cells(  data_df, 
                                time="__temp_time__", 
                                length=args.length_col, 
                                gfp=args.fp_col, 
                                cell_id="__temp_cell_id__", 
                                parent_id="__temp_parent_id__")
    
    cc_time = np.mean([cell.time[-1] - cell.time[0] for cell in cells_data])
    
    
    # OUs
    # Growth
    mean_lambda = estimate_mean_lambda(cells_data)
    if args.min_mean_lambda != None:
        mean_lambda = np.max([mean_lambda, args.min_mean_lambda])
        
    gamma_lambda = 1/(cc_time * args.tau_lambda)
    if args.min_gamma_lambda != None:
        gamma_lambda = np.max([gamma_lambda, args.min_gamma_lambda])
        
    var_lambda = (args.cv_lambda*mean_lambda)**2 *2* gamma_lambda
    if args.min_var_lambda != None:
        var_lambda = np.max([var_lambda, args.min_var_lambda])
    
    # Porduction
    mean_q = estimate_mean_q(cells_data, beta)
    if args.min_mean_q != None:
        mean_q = np.max([mean_q, args.min_mean_q])
        
    gamma_q = 1/(cc_time * args.tau_q)
    if args.min_gamma_q != None:
        gamma_q = np.max([gamma_q, args.min_gamma_q])
        
    var_q= (args.cv_q*mean_q)**2 *2* gamma_q
    if args.min_var_q != None:
        var_q = np.max([var_q, args.min_var_q])
    
    
    params = {}
    params["mean_lambda"] = ["bound", mean_lambda, mean_lambda/2., mean_lambda*2]
    params["gamma_lambda"] = ["bound", gamma_lambda, gamma_lambda/100., gamma_lambda*100.] 
    params["var_lambda"] = ["bound",  var_lambda, var_lambda/100., var_lambda*100.]


    params["mean_q"] = ["bound", mean_q, mean_q/2., mean_q*2]
    params["gamma_q"] = ["bound", gamma_q, gamma_q/100., gamma_q*100.] 
    params["var_q"] = ["bound",  var_q, var_q/100., var_q*100.]

    if args.opt_beta:
        params["beta"] =  ["bound", beta, beta/100, beta*100]

    else:
        params["beta"] =  ["fixed", beta]

    # measurment noise
    var_x = estimate_var_x(cells_data, mean_lambda)
    var_g = estimate_var_g(cells_data, mean_q, beta)
    
    params["var_x"] = ["bound", var_x, var_x/100., var_x*100.]
    params["var_g"] = ["bound", var_g, var_g/100., var_g*100.]

    # cell division
    params["var_dx"] = ["free", 1e-3]
    params["var_dg"] = ["free", 1]
    return params


        
            
def write_params2file(params, filename, step_scale=0.1):
    with open(filename, 'w') as fout:
        fout.write("# Automatically estimated parameter for initialzing MLE search\n")
        for param in params:
            if params[param][0] == "bound":
                fout.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format(param, params[param][1],
                                                                 params[param][1]*step_scale,
                                                                 params[param][2],
                                                                 params[param][3]) )
            elif params[param][0] == "free":
                fout.write("{:s} = {:.2E}, {:.2E}\n".format(param, params[param][1], params[param][1]*step_scale))
            elif params[param][0] == "fixed":
                fout.write("{:s} = {:.2E}\n".format(param, params[param][1]))
                

def write_config2file(args, filename):
    with open(filename, 'w') as fout:
        fout.write("# Automatically generated configuration file\n")
        fout.write("time_col = {:s}\n".format(args.time_col))
        fout.write("length_col = {:s}\n".format(args.length_col))
        fout.write("fp_col = {:s}\n".format(args.fp_col))
        
            
        if args.rescale_time!=None:
            fout.write("rescale_time = {:.3E}\n".format(args.rescale_time))
            
        if args.segment_col!=None:
            fout.write("segment_col = {:s}\n".format(args.segment_col))
            
        if args.filter_col!=None:
            fout.write("filter_col = {:s}\n".format(args.filter_col)) 
        
        fout.write("cell_tags = " + tags_str(args.cell_tags) +"\n")
        fout.write("parent_tags = " + tags_str(args.parent_tags)+"\n")
#################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Estimate init parameters from raw input files and write csv_config file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-d',
                        dest='dir',
                        help='Directory with input files. It will be searched for .csv files.',
                        required=True)

    # COLUMNS #
    parser.add_argument('-time_col',
                        dest='time_col',
                        help='Time column in the input file(s).',
                        default="time_sec",
                        required=False)
    
    parser.add_argument('-length_col',
                        dest='length_col',
                        help='Length column in the input file(s).',
                        default="length_um",
                        required=False)

    parser.add_argument('-fp_col',
                        dest='fp_col',
                        help='Fluorescence column in the input file(s).',
                        default="gfp_nb",
                        required=False)
    
    parser.add_argument('-filter_col',
                        dest='filter_col',
                        help='Filter column in the input file(s).',
                        default=None,
                        required=False)  
    
    parser.add_argument('-segment_col',
                        dest='segment_col',
                        help='Segment column in the input file(s).',
                        default=None,
                        required=False)  
    
    parser.add_argument('-cell_tags',
                        dest='cell_tags',
                        help='Cell tags specifies the columns that will be used to compose unique cell id.',
                        nargs='+',
                        default=["lane_ID", "id"],
                        required=False)  
    
    parser.add_argument('-parent_tags',
                        dest='parent_tags',
                        help='Parent tags specifies the columns that will be used to compose unique parent id. Must be consistent with cell tags.',
                        nargs='+',
                        default=["lane_ID", "parent_id"],
                        required=False)  
    
    
    # PARAMS #
    parser.add_argument('-tau_lambda',
                        dest='tau_lambda',
                        help='Time scale of growth rate noise (relative to cell cycle time)',
                        type=float,
                        default=0.5,
                        required=False)
    
    parser.add_argument('-cv_lambda',
                        dest='cv_lambda',
                        help='CV (noise level) growth rate',
                        type=float,
                        default=0.5,
                        required=False)    

    parser.add_argument('-tau_q',
                        dest='tau_q',
                        help='Time scale of production rate noise (relative to cell cycle time)',
                        nargs='+',
                        type=float,
                        default=0.2,
                        required=False)         
    
    parser.add_argument('-cv_q',
                        dest='cv_q',
                        help='CV (noise level) production',
                        type=float,
                        default=0.5,
                        required=False)    
    
    
    parser.add_argument('-beta',
                        dest='beta',
                        help='Bleaching rate',
                        nargs='+',
                        type=float,
                        default=[0.0],
                        required=False)   
    
    parser.add_argument('--opt_beta', help="Set beta to a bound parameter, i.e. RealTrace optimizes it.", action='store_true')

    # MIN PARAMETERS
    
    parser.add_argument('-min_mean_lambda',
                        dest='min_mean_lambda',
                        help='Minimal mean growth rate (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False)   
    
    parser.add_argument('-min_gamma_lambda',
                        dest='min_gamma_lambda',
                        help='Minimal gamma lambda (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False)  
    
    parser.add_argument('-min_var_lambda',
                        dest='min_var_lambda',
                        help='Minimal var_lambda (kick size squared) (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False) 
    
    parser.add_argument('-min_mean_q',
                        dest='min_mean_q',
                        help='Minimal mean production rate (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False)   
    
    parser.add_argument('-min_gamma_q',
                        dest='min_gamma_q',
                        help='Minimal gamma q (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False)  
    
    parser.add_argument('-min_var_q',
                        dest='min_var_q',
                        help='Minimal var_q (kick size squared) (mainly needed for starvation conditions)',
                        type=float,
                        default=None,
                        required=False)  

    # MISC#
    parser.add_argument('-rescale_time',
                        dest='rescale_time',
                        help='Rescale time to change unit. E.g. -rescale_time 60 will result in t->t/60',
                        type=float,
                        default=None,
                        required=False)    
    
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
        data[args.length_col] = data[args.length_col].astype(float)
        data[args.fp_col]    = data[args.fp_col].astype(float)
        data[args.time_col]   = data[args.time_col].astype(float)
    
          ####     
            
        if args.filter_col!=None:
            data = data[data[args.filter_col].astype(bool)]
        if args.segment_col!=None:
            segment_idxs = sorted(np.unique(data[args.segment_col]))
            
        ### Create tags
        data["__temp_cell_id__"] = ""
        for tag in args.cell_tags:
            data["__temp_cell_id__"] = data["__temp_cell_id__"] + data[tag] + "."
    
            
        data["__temp_parent_id__"] = ""
        for tag in args.parent_tags:
            data["__temp_parent_id__"] = data["__temp_parent_id__"] + data[tag] + "."
            

            
        if args.rescale_time!=None:
            data["__temp_time__"] = data[args.time_col].astype(float)/args.rescale_time
        else:
            data["__temp_time__"] = data[args.time_col].astype(float)
            
        if args.segment_col!=None:
            for segment in segment_idxs:
                if len(args.beta)==0:
                    beta = args.beta
                else:
                    beta = args.beta[int(segment)]
                params = estimate_parameters(data[data[args.segment_col]==segment], args, beta)
                
                # write to file named as expected by "run_*.py"
                parameter_file = infile[:-4] + "_segment" + str(segment) + "_parameters.txt"
                write_params2file(params, parameter_file)
        
        else:
            params = estimate_parameters(data, args, args.beta[0])
            
            # write to file named as expected by "run_*.py"
            parameter_file = infile[:-4] + "_parameters.txt"
            write_params2file(params, parameter_file)

    write_config2file(args, input_directory + "/csv_config.txt")
        
    return
       

if __name__ == "__main__":
    main()