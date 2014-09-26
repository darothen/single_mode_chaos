Title: PCE DAKOTA Utility and Analysis Scripts
Author: Daniel Rothenberg
Date: September 13, 2014
Comment: This document summarizes the scripts used to perform the probabilistic collocation method with DAKOTA to compute a polynomial chaos expansion of the parcel model.
Base Header Level: 3

## Overview

This repository archives the scripts used to perform polynomial chaos expansions using the probabilistic collocation method with the DAKOTA software. Please see the [File Summary] for details on each particular file, including the output with which each is associated.

## File Summary 

---

### PCE Scripts

1. ``dakota_poly.py``: Representation of polynomial output by DAKOTA PCE routines in order to analyze coefficients, save for future use, or evaluate directly in Python. 
2. ``dist_tools.py``: Collection of functions mapping various distributions to one another; *design_lhs_exp()* for performing LHS experiments.
3. ``gen_scripts.py``: Read a pre-configured PCE experiment from the data saved by the notebook and generate the DAKOTA and Python scripts for running the model.
4. ``process_outputs.py``: Read the output from a DAKOTA PCE experiment in order to extract fitted polynomial.
5. ``pce_experiment.ipynb``: Main driver for customizing/tweaking PCE experiment. 
    -- **OUTPUT**
        1. ``${exp_name}_exp.dict``: dictionary with the setup of the PCE experiment, including the experiment name, the variables (and their parameter-space definitions), the body/name of the function, and the directive passed to DAKOTA.
        2. ``config.p``: similar to above, but for a single particular PCE computation to be saved in its archive directory
        2. ``${exp_name}_results.dict``: mapping of the iterated experiment results from DAKOTA, indicating which timestamped folder in the overall **save/** directory archives the simulation results.
6. **Templates**
    - ``model_pce_parallel.template``: Run DAKOTA in batch mode with asynchronous parcel model evaluations, driving the model with file I/O mode
    - ``model_pce.template``: Use a direct Python interface to sequentially perform model evaluations
    - ``model_run_script.template``: Run the model by reading an input file and writing an output file in order to interface with DAKOTA
    - ``model_run_linked.template``: Run the model directly through a Python function call

### Analysis Scripts 

1. ``pce_sample.py``: Perform an LHS experiment over the parameter space defined for a simulation, running the model **n_samples** times. Saves the design (sample) points, their projection in *z*_space for evaluating with the PCE, and the results of evaluating hte model on the design points.
    - **OUTPUT** 
        1. ``${exp_name}_LHS_sample.npz``
2. ``pce_vis.py``: Create some visualizations detailing the performance of a given PCE against the LHS study. *Requires that ``pce_sample.py`` has been run on the same __${exp_name}!__* 
3. ``pcm_param.py``: Generate a portable HDF5 file containing the data necessary to evaluate a polynomial chaos expansion; also implements the logic to quickly retrieve and evaluate the parameterization
4. ``proc_sobol.py``: Process the output DAKOTA file for a set of experiments and produce DataFrames summarizing the Sobol indices for all the terms.
    - **OUTPUT**
        1. ``${exp_name}_sobol.dict``
        2. ``sobol.df`` in experiment save directory

### PCE Output Directory

1. ``config.p``:
2. ``model_run.py``:
3. ``model_pce.in``:
4. ``dakota.rst``:
5. ``model_pce.data``:
6. ``model_pce.out``:
7. ``coeffs_all``:
8. ``pce.p``:
9. ``sobol.df``: