Title: Single-mode Activation Parameterization using PCM / DAKOTA
Author: Daniel Rothenberg
Date: June 8, 2015

This repository archives the scripts used to perform polynomial chaos expansions using the probabilistic collocation method with the DAKOTA software. Please see the [File Summary] for details on each particular file, including the output with which each is associated.

For more information, please contact Daniel Rothenberg (darothen@mit.edu).

## Contents

The contents of this repository belong to four different categories:

1. Drivers
2. Post-processors
3. Experiment Scripts
4. Utilities

The **Drivers** are scripts which automate the process of performing a series of chaos expansions via DAKOTA, as well as the utility packages which aid it. **Post-processors** analyze the DAKOTA output and package its results into forms that can be studied or re-used elsewhere. A set of **Experiment Scripts** help to conduct, study, and visualize sampling experiments which aim to validate the derived chaos expansion. Finally, several **Utilities** are available to help clean up much of the superfluous results and otuputs which aren't critical to the experiment.

## Workflow

- Modify the notebook `pce_experiment.ipynb`; you should only need to modify the contents of the first cell underneath "Experiment Setup."
    + The second cell will record the experiment configuration in a `dict` serialized to disk
    + The next cell (underneath "Generate the model...") will analyze the saved `dict` and generate both the DAKOTA runtime configuration script and the Python interface to the model being studied, using templates saved in this repository; it will then conduct the requested chaos expansions
    + The cell under "Process output" will extract the chaos expansion from the DAKOTA output into aform that is readily re-usable through Python. It will also extract the computed Sobol` indices for later analysis
- Sample from the parameter space used to derive the chaos expansion
    + This is automated in `pce_experiment.ipynb` by calling the script `pce_sample.py`, although it can be run separately
    + On future iterations (e.g. if you change something high-level in the chaos expansion mechanics), rather than re-compute the sampling design, follow these steps:
        * copy *SM_OLS_LHS_design.csv* to *{exp_name}_LHS_design.csv*
        * run the script `pce_sample.py` with an explicit reference "-r" to this new design 
- Run the `postprocess.ipynb` notebook. This will archive all the results found for the list of experiments enumerated in the second cell of the notebook. It will also run `pcm_param.py` to create the python module and datastore as well as the ASCII files for evaluating the chaos expansions in Fortran.
- *Optionally*, use the command line utility `gen_nc.py` to create a self-contained netCDF file which encapsulates all the expansion details and can be used with the CESM-MARC.

## Output

The notebook `postprocess.ipynb` was designed to consolidate all the output into a single archive containing all the necessary utility scripts to evaluate a given chaos expansion and study it versus the original model and other parameerizations in a clean context. Originally, all output is cached with timestamps in the ``save/`` directory, and mapped through the saved Python `dict`s. For a given experiment `exp` with expansion orders derived up through *n*, this script will produce a folder structure like this:

```
    | exp/
    | --->/sample.csv/                
    | --->/--------->/expansion_order_1/
    | ...
    | --->/--------->/expansion_order_n/
    | --->/--------->/----------------->config.p
    | --->/--------->/----------------->model_run.py
    | --->/--------->/----------------->model_pce.in
    | --->/--------->/----------------->dakota.rst
    | --->/--------->/----------------->model_pce.dat
    | --->/--------->/----------------->model_pce.out
    | --->/--------->/----------------->resp_coeffs_all
    | --->/--------->/----------------->pce.p
    | --->/--------->/----------------->sobol.p
    | --->/--------->/----------------->exp_n.ascii
```

which contains all the intermediate output from DAKOTA and the sampling experiment. 

The finishing touch of `postprocessing.ipynb` is to package multiple experiments together in a single archive. This is a common task, say for when one feature of the chaos expansion (such as regression method) is varied. The end result of this operation is an archive with folders for each experiment like above, but with several additional files:

```
    | exp/
    | --->/exp_sampling_results.csv
    | --->/exp_sampling_tidy.csv
    | --->/exp_stats.p
    | --->/dakota_poly.p
    | --->/dist_tools.py
    | --->/pcm_param.py
    | --->/pcm_param.h5
```

The sampling results files here are pre-arranged in forms that are suitable for different plots / analyses (e.g. "tidy" datasets for use with Seaborn).

---

## File Summary

Note that deprecated or non-maintained scripts are denoted with *italics*. Outputs enumerated here are those directly created by various scripts, which may be intercepted later on by `postprocessing.ipynb` (see [above][Output])

### PCE Scripts

1. ``dakota_poly.py``: Representation of polynomial output by DAKOTA PCE routines in order to analyze coefficients, save for future use, or evaluate directly in Python. 
2. ``dist_tools.py``: Collection of functions mapping various distributions to one another; *design_lhs_exp()* for performing LHS experiments.
3. ``gen_scripts.py``: Read a pre-configured PCE experiment from the data saved by the notebook and generate the DAKOTA and Python scripts for running the model.

    **OUTPUT**
    - ``model_run.py``: Python interface to the model being analyzed via chaos expansion
    - ``model_pce.in``: Driver script for DAKOTA chaos expansion

4. ``process_outputs.py``: Read the output from a DAKOTA PCE experiment in order to extract fitted polynomial.
5. ``pce_experiment.ipynb``: Main driver for customizing/tweaking PCE experiment. 
    
    **OUTPUT**
    - ``${exp_name}_exp.dict``: dictionary with the setup of the PCE experiment, including the experiment name, the variables (and their parameter-space definitions), the body/name of the function, and the directive passed to DAKOTA
    - ``config.p``: similar to above, but for a single particular PCE computation to be saved in its archive directory
    - ``${exp_name}_results.dict``: mapping of the iterated experiment results from DAKOTA, indicating which timestamped folder in the overall **save/** directory archives the simulation results
   
6. ``gen_nc.py``: Convert chaos expansion output into a format readable by the CESM-MARC initialization routines
    **OUTPUT**
    - ``${exp_name}.nc``
7. **Templates**
    - ``model_pce_parallel.template``: Run DAKOTA in batch mode with asynchronous parcel model evaluations, driving the model with file I/O mode
    - ``model_pce.template``: Use a direct Python interface to sequentially perform model evaluations
    - ``model_run_script.template``: Run the model by reading an input file and writing an output file in order to interface with DAKOTA
    - ``model_run_linked.template``: Run the model directly through a Python function call

### Analysis Scripts 

1. ``pce_sample.py``: Perform an LHS experiment over the parameter space defined for a simulation, running the model **n_samples** times. Saves the design (sample) points, their projection in *z*_space for evaluating with the PCE, and the results of evaluating hte model on the design points.

    **OUTPUT** 
    - ``${exp_name}_LHS_design.csv``: A CSV containing the sample dataset points for all variables evaluating as part of the sampling study
    - ``${exp_name}_LHS_design_result.csv``: A CSV mapping the parameter samples from the design dataset to output for the model, chaos expansions, and alternative parameterizations

2. ``pce_vis.py``/``pce_vis_old.py``: Create some visualizations detailing the performance of a given PCE against the LHS study. *Requires that ``pce_sample.py`` has been run on the same __${exp_name}!__* 
3. ``pcm_param.py``: Generate a portable HDF5 file containing the data necessary to evaluate a polynomial chaos expansion; also implements the logic to quickly retrieve and evaluate the parameterization

    **OUTPUT**
    - ``pcm_param.h5``: an HDF5 file containing the coefficient and term order vectors/matrices for a given sets of experiments with given expansion orders

4. *proc_sobol.py* Process the output DAKOTA file for a set of experiments and produce DataFrames summarizing the Sobol indices for all the terms.

    **OUTPUT**
    - ``${exp_name}_sobol.dict``
    - ``sobol.df`` 

