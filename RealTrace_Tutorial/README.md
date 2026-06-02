# Readme for the RealTrace Tutorial

This is a tutorial illustrating the usage of RealTrace. For initial installation, please see the [installation wiki](https://github.com/nimwegenLab/RealTrace/wiki/2.-Installation).

## Jupyter notebook
The [Jupyter notebook](https://github.com/nimwegenLab/RealTrace/tree/main/RealTrace_Tutorial/RealTrace_Tutorial_notebook.ipynb) is a step-by-step guide for the usage of RealTrace. It also illustrates the use of RealTrace's output and performs some basic anaysis. This includes:
1. Plotting of the inferred cell dynamics 
2. Calculation of the correlation functions
3. Plotting sampled cell trajectories and illustrates their use by recalculating the correlation functions


## Minimal Example
If you want to run a minimal example, the `Minimal_example_data` folder contains all files to run RealTrace directly on a synthetic data set. To run RealTrace on this data run from this directory (expected runtime is around 2min depending on the machine):
```
../bin/RealTrace -i Minimal_example_data/input.csv -b Minimal_example_data/parameter_file.txt -c Minimal_example_data/csv_config.txt -m -p
```