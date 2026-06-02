# Readme for the RealTrace toolbox

## Setup
`setup.py` prepares everything for running RealRrace. This can be run on a single input file or a directory with many files (both specified via `-d`).

First, it outputs a csv_config file that will be used by RealTrace to interpret the the input file. 

Second, it estimates initial parameters, based on heuristics and rough estimates. That is:
1. The means of the OU processes (mean growth and productionr rate) are estimated naively (i.e. ignoring the measuremnt errors) and are typically very acurate.
2. The `gamma_lambda/q` and `var_lambda/q` are estimated based on the set CV and the time scale of the respective process (see options below). The defaults are  set to typical values we have seen in the past for mother machine data. The time scale is set relative to the mean cell cycle time which is estimated fromt the data. 
3. Measurement noise parameters are estimated (very crudely) by assuming the biol. dynamics is noise-free. This will overestimate the parameters, however, the order of magnitude should be correct. 
4. The bleaching rate needs to be set.
5. Parameters for cell division noise are set to `var_dx = 1e-3` and `var_dg = 1`.

Based on the initial values p0 for each parameter the initial step size for the optimization, and the bounds are set. That is `step = p0/10`, `lower bound = p0/100`, `upper bound = p0*100`.

Bounds (even very wide ones) have proven to be very useful in dealing with numerical issues during the optimization.    

The full list of options is:
```
usage: setup.py [-h] -d DIR [-time_col TIME_COL] [-length_col LENGTH_COL]
                [-fp_col FP_COL] [-filter_col FILTER_COL]
                [-segment_col SEGMENT_COL]
                [-cell_tags CELL_TAGS [CELL_TAGS ...]]
                [-parent_tags PARENT_TAGS [PARENT_TAGS ...]]
                [-tau_lambda TAU_LAMBDA] [-cv_lambda CV_LAMBDA]
                [-tau_q TAU_Q [TAU_Q ...]] [-cv_q CV_Q]
                [-beta BETA [BETA ...]] [--opt_beta]
                [-min_mean_lambda MIN_MEAN_LAMBDA]
                [-min_gamma_lambda MIN_GAMMA_LAMBDA]
                [-min_var_lambda MIN_VAR_LAMBDA] [-min_mean_q MIN_MEAN_Q]
                [-min_gamma_q MIN_GAMMA_Q] [-min_var_q MIN_VAR_Q]
                [-rescale_time RESCALE_TIME]

Estimate init parameters from raw input files and write csv_config file.

optional arguments:
  -h, --help            show this help message and exit
  -d DIR                Directory with input files. It will be searched for
                        .csv files. (default: None)
  -time_col TIME_COL    Time column in the input file(s). (default: time_sec)
  -length_col LENGTH_COL
                        Length column in the input file(s). (default:
                        length_um)
  -fp_col FP_COL        Fluorescence column in the input file(s). (default:
                        gfp_nb)
  -filter_col FILTER_COL
                        Filter column in the input file(s). (default: None)
  -segment_col SEGMENT_COL
                        Segment column in the input file(s). (default: None)
  -cell_tags CELL_TAGS [CELL_TAGS ...]
                        Cell tags specifies the columns that will be used to
                        compose unique cell id. (default: ['lane_ID', 'id'])
  -parent_tags PARENT_TAGS [PARENT_TAGS ...]
                        Parent tags specifies the columns that will be used to
                        compose unique parent id. Must be consistent with cell
                        tags. (default: ['lane_ID', 'parent_id'])
  -tau_lambda TAU_LAMBDA
                        Time scale of growth rate noise (relative to cell
                        cycle time) (default: 0.5)
  -cv_lambda CV_LAMBDA  CV (noise level) growth rate (default: 0.5)
  -tau_q TAU_Q [TAU_Q ...]
                        Time scale of production rate noise (relative to cell
                        cycle time) (default: 0.2)
  -cv_q CV_Q            CV (noise level) production (default: 0.5)
  -beta BETA [BETA ...]
                        Bleaching rate (default: [0.0])
  --opt_beta            Set beta to a bound parameter, i.e. RealTrace
                        optimizes it. (default: False)
  -min_mean_lambda MIN_MEAN_LAMBDA
                        Minimal mean growth rate (mainly needed for starvation
                        conditions) (default: None)
  -min_gamma_lambda MIN_GAMMA_LAMBDA
                        Minimal gamma lambda (mainly needed for starvation
                        conditions) (default: None)
  -min_var_lambda MIN_VAR_LAMBDA
                        Minimal var_lambda (kick size squared) (mainly needed
                        for starvation conditions) (default: None)
  -min_mean_q MIN_MEAN_Q
                        Minimal mean production rate (mainly needed for
                        starvation conditions) (default: None)
  -min_gamma_q MIN_GAMMA_Q
                        Minimal gamma q (mainly needed for starvation
                        conditions) (default: None)
  -min_var_q MIN_VAR_Q  Minimal var_q (kick size squared) (mainly needed for
                        starvation conditions) (default: None)
  -rescale_time RESCALE_TIME
                        Rescale time to change unit. E.g. -rescale_time 60
                        will result in t->t/60 (default: None)
```

Rescale the time to change the time unit (see also RealTrace Readme). Note that this changes the time unit of the initial parameters, too. E.g. let's say the input file uses seconds for the time column, then `-rescale_time 60` will transform this into minutes and all parameters (including the bleaching rate!) will be in units of minutes.

                        
                        

## Plot
`plot_results.py` plots the inferred RealTrace traces fetching the prediction files in the directory specified by `-d`. This code also includes some useful functions and classes that can be used for further analysis. 
```
usage: plot_results.py [-h] -d DIR [DIR ...]
                       [-time_unit TIME_UNIT [TIME_UNIT ...]]
                       [-n_random_cells N_RANDOM_CELLS] [-r RANGE [RANGE ...]]
                       [-o OUT [OUT ...]] [--replot]

Plot predicted traces of example cell cycles.

optional arguments:
  -h, --help            show this help message and exit
  -d DIR [DIR ...]      Directory(-ies) that will be searched for prediction
                        files (recursively) (default: None)
  -time_unit TIME_UNIT [TIME_UNIT ...]
                        Time unit and rescaling of the time eg t/60 (default:
                        ['min', '1'])
  -n_random_cells N_RANDOM_CELLS
                        Plot n randoms cells (overwrites range parameters)
                        (default: 20)
  -r RANGE [RANGE ...]  Range for cells that will be plotted (start, stop,
                        step) as used by np.arange() (default: [None, None,
                        1])
  -o OUT [OUT ...]      Output directory, rather than same as prediction file
                        (default: [])
  --replot              Replot plots for files, that already exist (default:
                        False)
```

Requirements: `setup.py` needs numpy and pandas (tested with `np.__version__ '1.19.1'` and `pd.__version__
'1.3.5'`). In additon, `plot_result.py` needs maplotlib (tested with `mpl.__version__ '3.5.3'`).