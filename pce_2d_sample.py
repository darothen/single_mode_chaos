
import os, pickle
import numpy as np
import pandas as pd
import h5py

from functools import partial
from itertools import combinations

import model_run

from dist_tools import map_transfer_fcns

exp_name = "pcm_ols_parcel"
n_samples_1d = int(50. + 2)

PARALLEL = True
recompute = True

###########################################################

## Unload the configured experiment
print "Generating samples for %s" % exp_name
exp_dict = pickle.load(open("%s_exp.dict" % exp_name, 'r'))

directive_base = exp_dict['directive_base']

variables = exp_dict['variables']
var_names = [v[0] for v in variables]
dists = [v[3][0] for v in variables]
offsets = [v[4] for v in variables]

mapped = map_transfer_fcns(dists, directive_base, True)
maps = { var: m for (var, m) in zip(var_names, mapped) }

mean_dict = {}
range_dict = {}
dist_dict = {}
for v in variables:
    vn = v[0]
    a, b = v[3][1:]
    dist_dict[vn] = [a, b]
    v_mean = np.mean([a, b])
    print vn, a, v_mean, b
    mean_dict[vn] = v_mean
    if vn.startswith("ln"):
        range_dict[vn] = np.log(np.logspace(a, b, num=n_samples_1d, base=np.e))[1:-1]
    else:
        range_dict[vn] = np.linspace(a, b, num=n_samples_1d)[1:-1]

## Compute the analytical maps if desired
sample_fn = "%s_2D_sens.h5" % exp_name
if not os.path.exists(sample_fn) or recompute:
    RUN_SIMS = True
    print "writing new data to", sample_fn
    store = pd.HDFStore(sample_fn, "w")
    #store = h5py.File(sample_fn, "w")
else:
    RUN_SIMS = False
    print "reading data from", sample_fn
    store = pd.HDFStore(sample_fn, 'r')
    #store = h5py.File(sample_fn, "r")

fn = model_run.__dict__[exp_dict['function_name']]

if PARALLEL:
    from IPython.parallel import Client, require
    client = Client(profile='ssh')
    dv = client[:]
    dv['fn'] = fn
    cwd = os.getcwd()
    dv.execute("import sys; sys.path.append('%s'); import model_run" % cwd)


for i, combo in enumerate(combinations(var_names, 2)):
    #if i > 1: break
    var_a, var_b = combo
    print var_a, var_b

    range_a, range_b = range_dict[var_a], range_dict[var_b]
    
    other_means = {}
    for v in var_names:
        if v not in combo:
            other_means[v] = mean_dict[v]

    grid_a, grid_b = np.meshgrid(range_a, range_b)
    grid_a_lin = grid_a.ravel()
    grid_b_lin = grid_b.ravel()

    param_sets = [ { var_a: a, var_b: b } \
                   for a, b in zip(grid_a_lin, grid_b_lin) ]

    for p in param_sets:
        p.update(other_means)

    if RUN_SIMS:
        if PARALLEL:
            print "Executing in parallel"
            view = client.load_balanced_view()
            fn_par = lambda z : fn(**z)
            results = view.map_async(fn_par, param_sets, ordered=True)
            results.wait_interactive()

        else:
            results = [fn(**z) for z in param_sets]

        fn_nact = lambda z, smax: fn(fn_toggle=smax, **z)
        zipped = zip(param_sets, results)
        Nacts = np.array([fn_nact(z, np.exp(smax)) for z, smax in zipped])
        Ns = np.exp(np.array([p['lnN'] for p in param_sets]))
        act_fracs = Nacts/Ns

        ## Saving the full parameter set matrix
        parameters = pd.DataFrame(param_sets)
        store.put("%s/parameters" % ("_".join(combo), ), parameters)

        z_parameters = parameters.copy()
        for key in mean_dict:
            col = z_parameters[key]
            z_vals = maps[key](col, *dist_dict[key])
            z_parameters[key] = z_vals[:]
        store.put("%s/z_parameters" % ("_".join(combo), ), z_parameters)

        results = pd.Series(results, index=parameters.index)
        store.put("%s/smaxes" % ("_".join(combo), ), results)

        act_fracs = pd.Series(act_fracs, index=parameters.index)
        store.put("%s/act_fracs" % ("_".join(combo), ), act_fracs)

store.close()
