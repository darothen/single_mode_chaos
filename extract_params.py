#!/usr/bin/env python
"""
Extract parameters from a CESM/MARC simulation to use to compare a 
chaos expansion to parcel model output.
"""

import os, pickle
from itertools import product

import numpy as np
import pandas as pd

import netCDF4 as nc

import argparse
parser = argparse.ArgumentParser(description="Extract aerosol activation calculation parameters" +
                                             " from a CESM/MARC simulation")
parser.add_argument("exp_name", type=str, help="name of PCM experiment to reference")
parser.add_argument("-s", "--subset", action='store_true',
                    help="Subset CESM data based on PCM bounds")
parser.add_argument("-f", "--cesm_file", type=str, required=True,
                    help="CESM/MARC output data to extract from")
parser.add_argument("-n", "--n", type=int, default=10000, 
                    help="Number of samples to extract from the dataset")
parser.add_argument("-o", "--output", type=str, default="CESM_extract.csv",
                    help="Name of output file")

#: Aerosol modes to extract
MARC_modes = ["ACC", "MOS", "MBS", ]
DST_modes = ["DST01", "DST02", "DST03", "DST04"]
SSLT_modes = ["SSLT01", "SSLT02", "SSLT03", "SSLT04"]

#: Meteorology fields to extract
meteo_fields = ["Q", "P", "T", "WSUB"]

#: water vapor threshold
QCSMALL = 1e-8

#: Mapping of PCM var shorthand to CESM fields
mapping = {
    "V": "WSUB",
}

if __name__ == "__main__":

    args = parser.parse_args()
    print "Referencing PCM %s" % args.exp_name
    print "Extracting (%d) samples from %s" % (args.n, args.cesm_file)
    print "--"*20 + "\n"

    ###########################################################

    ## Load in the experiment dictionary, which contains information on 
    ## the parameter space structure and contents
    exp_dict = pickle.load(open("%s_exp.dict" % args.exp_name, 'r'))

    directive_base = exp_dict['directive_base']
    variables = exp_dict['variables']
    responses = exp_dict['responses']

    print "Reading in activation parameter fields"

    aer_prefixes = ["n", "m", "mu"]

    d = nc.Dataset(args.cesm_file, 'r')
    data = {}
    for (prefix, mode) in product(aer_prefixes, MARC_modes+DST_modes+SSLT_modes):
        if prefix == "mu" and mode in DST_modes + SSLT_modes: continue
        varname = prefix+mode
        print "   ", varname
        data[varname] = d.variables[varname][:]
    if "MOS" in MARC_modes:
        data["mOIM"] = d.variables["mOIM"][:]
        data['rhoMOS'] = data['mMOS']/(data['mOIM']/2e3 + (data['mMOS'] - data['mOIM'])/1.8e3)
        data["kappaMOS"] = (1e-9*data["mOIM"]/2e3 + 0.507*(data['mMOS'] - data['mOIM'])/1.8e3)/ \
                           (data['mMOS']/data['rhoMOS'])
    for field in meteo_fields:
        data[field] = d.variables[field][:]
        if field == "P": 
            data[field] = np.swapaxes(data[field], 1, 3)
            data[field] = np.swapaxes(data[field], 2, 3)
    d.close()

    ## Subset the variables to match the PCM setup
    print "Mapping the PCM vars to the CESM data"
    pcm_vars = [v[0] for v in variables]
    data_mapping = {}
    logvars = {}
    for var in pcm_vars:            
        logvar = "log" in var
        logvars[var] = logvar
        print "   ", var, " -> ", 

        ## Some logic to map the PCM varnaming scheme to the CESM modes
        if len(var) == 1:
            var_short = var
        else:
            bits = var.split("_")
            if len(bits) == 1 and logvar:
                var_short = bits[0][3:]
            elif len(bits) > 1 and logvar:
                prefix, var_short = var.split("_")
                prefix = prefix[3:]
                var_short = prefix.lower()+var_short
            else:
                var_short = var.replace("_", "")
        if var_short in mapping:
            var_short = mapping[var_short]

        if logvar: 
            data[var_short] = np.log10(data[var_short])

        print var_short
        data_mapping[var] = var_short

    ## Mask where QC is too small, kappa is invalid
    mask = (data['Q'] > QCSMALL) #& (data['kappaMOS'] > 0.1)
    n_tot = len(data['Q'][mask])

    row = "   {:<14s} {:<12d}"
    print row.format("valid data points", n_tot)
    if args.subset: print "Subsetting data based on PCM bounds"
    for v in variables:
        pcm_var = v[0]
        _, lo, hi = v[3]
        mask = mask & (data[data_mapping[pcm_var]] > lo) & \
                      (data[data_mapping[pcm_var]] < hi)
        print row.format(pcm_var, len(data['Q'][mask]))
    n_tot = len(data['Q'][mask])


    ## Sample the data
    print "Sampling the data"
    inds = np.random.choice(range(n_tot), size=args.n)
    data_subset = {}
    for key, field in data.iteritems():
        print "  ",key
        assert field.shape == mask.shape
        data_subset[key] = field[mask][inds]

    ####
    print "Extracting data from samples based on mapping"
    data_extract = {}
    for pcm_var, cesm_var in data_mapping.iteritems():
        data_extract[pcm_var] = data_subset[cesm_var]

    data_extract = pd.DataFrame(data_extract)

    ## Finally print a table describing the ranges of the extracted vars
    print "\n" + "--"*30  
    print "Distributions of extracted data" + "\n"
    header = " "*14 + " {:<10s}"*4
    print header.format('min', 'median', 'mean', 'max')
    
    row = "{:<14s}" + " {:<10.2f}"*4
    for var in data_extract:
        var_dat = data_extract[var]
        row_prnt = [var, np.min(var_dat), np.median(var_dat),
                    np.mean(var_dat), np.max(var_dat)]
        print row.format(*row_prnt)
    
    print ""
    
    print "Saving to file %s" % args.output
    data_extract.to_csv(args.output, )


