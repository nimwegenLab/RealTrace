# Readme for the RealTrace toolbox

## 1. Setup
The `setup.py` script prepares everything for running RealRrace. 

Input directory:
  -d DIR                directory with input files
  
Column setting specifying which columns will be used here and for RealTrace (see also RealTrace Readme):
  -time_col TIME_COL    time column
  -length_col LENGTH_COL
                        length column
  -fp_col FP_COL        fluorescence column
  -filter_col FILTER_COL
                        filter column
  -segment_col SEGMENT_COL
                        segment column
  -cell_tags CELL_TAGS [CELL_TAGS ...]
                        cell_tags
  -parent_tags PARENT_TAGS [PARENT_TAGS ...]
                        parent_tags
                        
Options for the initial parameter estimates:                 
  -tau_lambda TAU_LAMBDA
                        time scale growth (rel. to cell cycle time)
  -cv_lambda CV_LAMBDA  CV growth rate
  -tau_q TAU_Q [TAU_Q ...]
                        time scale prod. (rel. to cell cycle time)
  -cv_q CV_Q            CV production
  
Bleaching rate:
  -beta BETA [BETA ...]
                        bleaching rate
  --free_beta           Set beta to 'free'
  
Set minimum for the initial parameters (this is mainly relavant for starvation conditions where growth vanishes):
  -min_mean_lambda MIN_MEAN_LAMBDA
                        minimal mean growth rate
  -min_gamma_lambda MIN_GAMMA_LAMBDA
                        minimal gamma lambda
  -min_var_lambda MIN_VAR_LAMBDA
                        minimal var_lambda (kick size squared)
  -min_mean_q MIN_MEAN_Q
                        minimal mean production rate
  -min_gamma_q MIN_GAMMA_Q
                        minimal gamma q
  -min_var_q MIN_VAR_Q  minimal var_q (kick size squared)
  
Rescale the time to change the time unit (see also RealTrace Readme). Note that this changes the time unit of the initial parameters too. E.g. let's say the input file uses seconds for the time column, then `-rescale_time 60` will transform this into minutes and all parameters (including the bleaching rate!) will be in units of minutes:
  -rescale_time RESCALE_TIME
                        rescale time
                        
                        
                        
## 2. Run
Run RealTrace for a single input file or all input files in input directory:

  -d DIR          directory with input files or path to input file
  -bin BIN        location of RealTrace binary
  -suffix SUFFIX  suffix of parameter file (for reused parameter files)
  --dryrun        Shows what will be done
  --cluster       Submit job to cluster (do not run directly)
  
RealTrace options (see also RealTrace Readme):
  -o OUT          output dir (None)
  -space SPACE    search space, log or linear (log)
  -t TOL          Tolerance of maximization (1e-15)
  -m              Run maximization
  -p              Run prediction
  -j              Run joints

