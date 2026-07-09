import argparse
import os 
import sys


import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
# from matplotlib.patches import Patch
import matplotlib as mpl
# from scipy import stats

import glob

# ---------------------------------------------------------------------------------------------------------- #
# GGP_cell class to read prediction_files easily
# ---------------------------------------------------------------------------------------------------------- #
class GGP_cell:
    def __init__(self, cell_id = 0, parent_id=-1):
        self.parent_id = parent_id
        self.cell_id = cell_id
        self.log_length = []
        self.gfp = []
        self.time = []

        self.mean_x = []
        self.mean_g = []
        self.mean_l = []
        self.mean_q = []

        self.cov_xx = []
        self.cov_gg = []
        self.cov_ll = []
        self.cov_qq = []

        self.x = []
        self.g = []
        self.l = []
        self.q = []

def df2ggp_cells(dataset, 
            time="time", 
            log_length="log_length", gfp="fp", 
            mean_x="mean_x", mean_g="mean_g", 
            mean_l="mean_l", mean_q="mean_q",
            cov_xx="cov_xx",
            cov_gg="cov_gg",
            cov_ll="cov_ll",
            cov_qq="cov_qq",
            cell_id="cell_id", 
            parent_id="parent_id"):
    """ 
    dataset (pandas data frame as read from csv file) to list of GGP_cell instances, m
    written for ggp output
    """
    cell_list = []
    last_cell = ""
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            new_cell = GGP_cell(
                        cell_id=row[cell_id], 
                        parent_id=row[parent_id])
            cell_list.append(new_cell)

        cell_list[-1].log_length.append(row[log_length])
        cell_list[-1].gfp.append(row[gfp])
        cell_list[-1].time.append(row[time])

        cell_list[-1].mean_x.append(row[mean_x])
        cell_list[-1].mean_g.append(row[mean_g])
        cell_list[-1].mean_l.append(row[mean_l])
        cell_list[-1].mean_q.append(row[mean_q])

        cell_list[-1].cov_xx.append(row[cov_xx])
        cell_list[-1].cov_gg.append(row[cov_gg])
        cell_list[-1].cov_ll.append(row[cov_ll])
        cell_list[-1].cov_qq.append(row[cov_qq])

        last_cell = row[cell_id]
    return cell_list

def df2cells(dataset, 
            time="time", 
            x="x", g="g", 
            l="l", q="q",
            cell_id="cell_id", 
            parent_id="parent_id"):
    """ 
    dataset (pandas data frame as read from csv file) to list of GGP_cell instances, m
    written for ggp output
    """
    cell_list = []
    last_cell = ""
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            new_cell = GGP_cell(
                        cell_id=row[cell_id], 
                        parent_id=row[parent_id])
            cell_list.append(new_cell)

        cell_list[-1].time.append(row[time])

        cell_list[-1].x.append(row[x])
        cell_list[-1].g.append(row[g])
        cell_list[-1].l.append(row[l])
        cell_list[-1].q.append(row[q])

        last_cell = row[cell_id]
    return cell_list



def header_lines(filename, until="cell_id"):
    with open(filename,'r') as fin:
        for i, line in enumerate(fin):
            if line.startswith(until):
                return i
    return None

def plot_predictions(filename, 
                     start=None, stop=None, step=None, 
                     n_random_cells=None,
                     time_unit=("min", 60), 
                     xlim=[None, None], 
                     outfile=None, 
                     show=True, 
                     color_scheme="tab"):
    
    """ needs a prediction file, start, stop, step refers to cells """
    fig, axes = plt.subplots(4, 1, figsize=(8,10), sharex=True)
    ax = axes.ravel()

    for a in ax:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

        a.spines['right'].set_color('none')
        a.yaxis.tick_left()

        a.spines['top'].set_color('none')
        a.xaxis.tick_bottom()
        
    data = pd.read_csv(filename, skiprows= header_lines(filename, until="cell_id"))
    all_cells = df2ggp_cells(data)

    np.random.seed(0)
    if n_random_cells != None:
        cells_data = np.random.choice(all_cells, size=n_random_cells, replace=False)
    else:
        cells_data = all_cells[start: stop: step]
     

    norm = mpl.colors.Normalize(vmin=-len(cells_data)/2, vmax=len(cells_data))
    if len(cells_data)==1:
        norm = mpl.colors.Normalize(vmin=-10*len(cells_data), vmax=len(cells_data))

    cmap_data = mpl.cm.ScalarMappable(cmap='Oranges', norm=norm)
    cmap_data.set_array([])

    cmap_prediction = mpl.cm.ScalarMappable(cmap='Blues', norm=norm)
    cmap_prediction.set_array([])
    colors= ['tab:blue',
            'tab:orange',
            'tab:green',
            'tab:red',
            'tab:purple',
            'tab:brown',
            'tab:pink',
            'tab:gray',
            'tab:olive',
            'tab:cyan']*int(len(cells_data)/10+1)

    s = 0.5
    lw = 0.7
    for i, cell in enumerate(cells_data):
        if color_scheme=="tab":
            data_color =colors[i]
            prediction_color =colors[i]
        else:
            data_color = cmap_data.to_rgba(i)
            prediction_color = cmap_prediction.to_rgba(i)

        time = np.array(cell.time) / time_unit[1]

        ax[0].scatter(time, cell.log_length, color=data_color, s=s)
        ax[0].plot(time, cell.mean_x, color=prediction_color, lw=lw)
        ax[0].fill_between(time, cell.mean_x-np.sqrt(cell.cov_xx), cell.mean_x+np.sqrt(cell.cov_xx), 
                    color=prediction_color, alpha=0.2)

        ax[1].scatter(time, cell.gfp, color=data_color, s=s)
        ax[1].plot(time, cell.mean_g, color=prediction_color, lw=lw)
        ax[1].fill_between(time, cell.mean_g-np.sqrt(cell.cov_gg), cell.mean_g+np.sqrt(cell.cov_gg), 
                    color=prediction_color, alpha=0.2)

        ax[2].plot(time, np.array(cell.mean_l)*time_unit[1], color=prediction_color, lw=lw)
        ax[2].fill_between(time, 
                           (cell.mean_l-np.sqrt(cell.cov_ll))*time_unit[1], 
                           (cell.mean_l+np.sqrt(cell.cov_ll))*time_unit[1], 
                            color=prediction_color, alpha=0.2)

        ax[3].plot(time, np.array(cell.mean_q)*time_unit[1], color=prediction_color, lw=lw)
        ax[3].fill_between(time, 
                           (cell.mean_q-np.sqrt(cell.cov_qq))*time_unit[1], 
                           (cell.mean_q+np.sqrt(cell.cov_qq))*time_unit[1], 
                            color=prediction_color, alpha=0.2)

    ax[0].set_ylabel("log cell size")   
    ax[1].set_ylabel("total GFP")   
    ax[2].set_ylabel(r"growth rate $\lambda$")   
    ax[3].set_ylabel(r"production rate $q$")   


    ax[3].set_xlabel("time ({:s})".format(time_unit[0]))  

    for i in range(4):
        ax[i].set_xlim(xlim)
        
    plt.tight_layout()

    if outfile != None:
        fig.savefig(outfile, dpi=300)
    elif show:
        plt.show()
    fig.clear()
    plt.close(fig)
    return ax

# ==================================================== #

def mk_missing_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory) 
    return directory

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

    return final_files


def get_prediction_files(path):
    if os.path.isfile(path):
        return [path]
    candidates = sorted(glob.glob(path + '/**/*prediction.csv', recursive=True))
    # ignore directories that end with '.czi' and only select files!
    return [c for c in candidates if os.path.isfile(c)]

########################################################################################################################
########################################################################################################################
########################################################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Plot predicted traces of example cell cycles.',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-d',
                        dest='dir',
                        help='Directory(-ies) that will be searched for prediction files (recursively)',
                        nargs='+',
                        type=str,
                        required=True)

    parser.add_argument('-time_unit',
                        dest='time_unit',
                        help='Time unit and rescaling of the time eg t/60',
                        nargs='+',
                        default=["min", "1"],
                        type=str,
                        required=False)
    
    parser.add_argument('-n_random_cells',
                        dest='n_random_cells',
                        help='Plot n randoms cells (overwrites range parameters)',
                        type=int,
                        default=20,
                        required=False)
    
    parser.add_argument('-r',
                        dest='range',
                        help='Range for cells that will be plotted (start, stop, step) as used by np.arange()',
                        nargs='+',
                        type=int,
                        default=[None, None, 1],
                        required=False)
    
    parser.add_argument('-o',
                        dest='out',
                        help='Output directory, rather than same as prediction file',
                        nargs='+',
                        type=str,
                        default=[],
                        required=False)

    parser.add_argument('--replot', help="Replot plots for files, that already exist", action='store_true')


    args = parser.parse_args()
    
    if len(args.out)>0:
        if len(args.out)!=len(args.dir):
            print("Different number of input and output directories. Stop")
            return
    else:
        for o in args.out:
            mk_missing_dir(o)
            
            
    # ======================================== #
    # ======================================== #
    for i, directory in enumerate(args.dir):
        input_files = get_prediction_files(directory)

        for infile in input_files:
            print(infile)

            if len(args.out)>0:
                out_file = os.path.join(args.out[i], infile.split("/")[-1][:-4]) + ".png"
            else:
                out_file = infile[:-4] + ".png"

            if args.replot or not os.path.exists(out_file):
                plot_predictions(infile, 
                                start=args.range[0], stop=args.range[1], step=args.range[2], 
                                n_random_cells=args.n_random_cells,
                                time_unit=(args.time_unit[0], float(args.time_unit[1])), 
                                xlim=[None, None], 
                                outfile=out_file, 
                                show=False)


# ================================================================================ #
if __name__ == "__main__":
    main()