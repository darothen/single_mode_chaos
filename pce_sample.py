import os, pickle
import numpy as np
from functools import partial

import model_run
from dist_tools import design_lhs_exp, map_transfer_fcns

exp_name = "tri_modal_ols"
CESM_SAMPLE = True
CESM_FILENAME = "pce_params.csv"
n_samples = int(1e4)

PARALLEL = True
recompute = True

z_func = lambda z: np.exp(z)
#z_func = lambda z: z

###########################################################

## Unload the configured experiment
print "Generating samples for %s" % exp_name

exp_dict = pickle.load(open("%s_exp.dict" % exp_name, 'r'))

directive_base = exp_dict['directive_base']

variables = exp_dict['variables']
print variables

dists = [v[3][0] for v in variables]
offsets = [v[4] for v in variables]


if CESM_SAMPLE:
    import pandas as pd
    df = pd.read_csv("pce_params.csv")
    design = df[[v[0] for v in variables]].values.T
    z_design = np.empty_like(design)
    for i, v in enumerate(variables):        
        dist, a, b = v[3]
        z_design[i, :] = maps[i](design[i, :], a, b)
    #z_design = z_design.T
else:
	maps = map_transfer_fcns(dists, directive_base, True)
    design, z_design = design_lhs_exp(variables, maps, offsets, int(n_samples))

for i in xrange(len(variables)):
    print variables[i][0], variables[i][3][1], np.min(design[i, :]), \
          np.max(design[i, :]), variables[i][3][2]

## Compute the analytical maps if desired
sample_fn = "%s_LHS_sample.npz" % exp_name
if not os.path.exists(sample_fn) or recompute:
    print "writing new data to", sample_fn

    fn = model_run.__dict__[exp_dict['function_name']]
    print fn

    if PARALLEL:
        
        from IPython.parallel import Client, require
        client = Client(profile='ssh')
        dv = client[:]
        dv['z_func'] = z_func
        dv['fn'] = fn
        fn_par = lambda z : fn(*z)
        dv['fn_par'] = fn_par

        cwd = os.getcwd()
        dv.execute("import sys; sys.path.append('%s'); import model_run" % cwd)

        print "Executing Smax calculation in parallel"
        view = client.load_balanced_view()
        results = view.map_async(fn_par, design.T, ordered=True)
        results.wait_interactive()

    else:
        results = [fn(*z) for z in design.T]

    results = np.array(results[:])

    ## Drop bad simulations
    good_design = design.T[results != -9999.0]
    good_z_design = z_design.T[results != -9999.0]
    good_results = results[results != -9999.0]
    zipped = zip(good_design, z_func(good_results))

    fn_nact = lambda z, smax : fn(*z, fn_toggle=smax)
    if PARALLEL:
        print "Executing Nact calculation in parallel"
        dv['fn_nact'] = fn_nact
        Nacts = view.map_async(fn_nact, zipped, ordered=True)
        Nacts.wait_interactive()

    else:
        Nacts = np.array([z, z_func(smax))) for z, smax in zipped])

    Nacts = np.array(Nacts[:])
    np.savez(sample_fn, results=good_results, design=good_design.T, z_design=good_z_design.T,
             Nacts=Nacts)

