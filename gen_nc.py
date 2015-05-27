#!/usr/bin/env python
""" Condense post-processed chaos expansion output into a single netCDF 
output file """

import h5py
import netCDF4 as nc
import numpy as np
import pandas as pd

import argparse
import os
import pickle
import time
import subprocess

parser = argparse.ArgumentParser(
    description="Utility script to save chaos expansion output in netCDF form"
)
parser.add_argument("exp_name", type=str, help="name of PCM experiment to save")
parser.add_argument("order", type=int, help="order of expansion from PCM experiment to save"),

def get_git_versioning():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

if __name__ == "__main__":

    args = parser.parse_args()
    print "Processing %s - order %d" % (args.exp_name, args.order)

    ## Get the dictionary with meta-data from the experiment, and
    ## extract the variable names, bounds, and whether they're log
    config_filepath = os.path.join(args.exp_name, 
                                   "expansion_order_%d" % args.order, 
                                   "config.p")
    with open(config_filepath, 'r') as f:
        config = pickle.load(f)
    variables = config['variables']
    names, bounds, logs = [], [], []
    for v in variables:
        name = v[0]
        dist = v[3]
        _, low, hi = dist # unpack dist tuple

        ## Collect this variable's data
        names.append(name)
        bounds.append([low, hi])
        logs.append(name.startswith("log"))
    names = np.array(names, dtype=str)
    bounds = np.array(bounds, dtype=float)
    logs = np.array(logs, dtype=bool)

    ## Get the vector of PCE coefficients and the order matrix
    store = h5py.File("pcm_param.h5", "r")
    coeffs = store[args.exp_name]['expansion_order_%d' % args.order]['coeffs'][:]
    orders = store[args.exp_name]['expansion_order_%d' % args.order]['orders'][:]
    store.close()

    nterms, nvars = orders.shape

    ## Using a structured DataArray, save this information into a netCDF file
    with nc.Dataset(args.exp_name + "_%d.nc" % args.order, 'w', format='NETCDF4') as root:
        print "   writing to file...",

        ## Create and populate the variables
        variable = root.createDimension('variable', nvars)
        term = root.createDimension('term', nterms)

        ## Add and populate the fields of the netCDF file
        varnames = root.createVariable('varname', names.dtype, ('variable',))
        varnames[:] = names
        varnames.description = "name of variable used in expansion"

        lower_bnds = root.createVariable('lower_bnd', bounds.dtype, ('variable', ))
        lower_bnds[:] = bounds[:, 0]   
        lower_bnds.description = "lower bound of variable uniform dist"     

        upper_bnds = root.createVariable('upper_bnd', bounds.dtype, ('variable', ))
        upper_bnds[:] = bounds[:, 1]
        upper_bnds.description = "upper bound of variable uniform dist"  

        log_vars = root.createVariable('log', 'i1', ('variable', )) # promote to NC_BYTE
        log_vars[:] = logs
        log_vars.description = "log10(variable) used in expansion"

        orders_var = root.createVariable('order', orders.dtype, ('term', 'variable'))
        orders_var[:] = orders
        orders_var.description = "chaos term/variable expansion order"

        coeffs_var = root.createVariable('coeff', coeffs.dtype, ('term', ))
        coeffs_var[:] = coeffs
        coeffs_var.descripption = "chaos term coefficient"

        ## Assign some global attributes
        git_commit = get_git_versioning().rstrip()
        root.source = "gen_nc.py | git@%s" % git_commit
        root.history = "Created " + time.ctime(time.time())
        root.description = \
            "Concise order %1d chaos expansion output for experiment %s" % (args.order, args.exp_name)

        print "done."
