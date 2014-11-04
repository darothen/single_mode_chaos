import os, pickle
import numpy as np
import pandas as pd
from functools import partial

import model_run
from dist_tools import design_lhs_exp, map_transfer_fcns

import argparse

parser = argparse.ArgumentParser(description="Conduct sampling study of parcel model and chaos expansion metamodels.")
parser.add_argument("exp_name", type=str, help="name of PCM experiment to reference")
parser.add_argument("-r", "--reference", type=str,  
                    help="Reference set of design points to use")
parser.add_argument("--params", action='store_true', 
                    help="Include sampling of ARG/MBN parameterizations")
parser.add_argument("--n", "-n", type=int, default=10000,
                    help="Number of samples to perform; either total for LHS or " + \
                         "limit for CESM")
parser.add_argument("--parallel", action='store_true', 
                    help="Perform evaluations in parallel")
parser.add_argument("--recompute", action='store_true', 
                    help="Force recompute/overwrite of existing results")

z_func = lambda z: 10.**z

if __name__ == "__main__":

    args = parser.parse_args()
    print "Studying %s" % args.exp_name
    print "--"*20
    if args.reference:
        print " referencing sample %s (max n - %d)" % (args.reference, args.n)
    else:
        print " computing new LHS design (max n - %d)" % args.n

    if args.params:
        print " analyzing MBN/ARG parameterizations"
    if args.parallel:
        print "running in parallel"
        from IPython.parallel import Client
        client = Client(profile='ssh')
        dv = client[:]

    print "\n"

    ###########################################################

    ## Load in the experiment dictionary, which contains information on 
    ## the parameter space structure and contents
    exp_dict = pickle.load(open("%s_exp.dict" % args.exp_name, 'r'))

    directive_base = exp_dict['directive_base']
    variables = exp_dict['variables']
    responses = exp_dict['responses']
    dists = [v[3][0] for v in variables]
    offsets = [v[4] for v in variables]

    ## z-space mapping functions
    maps = map_transfer_fcns(dists, directive_base, True)

    ## Was there a reference design set? If so, use that 
    if args.reference:
        ref_design_name = args.references
        stored_design = pd.read_csv(args.reference)
        design = stored_design[[v[0] for v in variables]].values.T
        z_design = np.empty_like(design)
        for i, v in enumerate(variables):
            dist, a, b = v[3]
            z_design[i, :] = maps[i](design[i, :], a, b)
    else:
        design, z_design = design_lhs_exp(variables, maps, offsets, args.n)
        design_df = pd.DataFrame(design.T, columns=[v[0] for v in variables])
        ref_design_name = "%s_LHS_design.csv" % args.exp_name
        design_df.to_csv(ref_design_name)

    print "                Mins             Maxes"
    print "           range    design  design    range"
    for i in xrange(len(variables)):
        print variables[i][0], variables[i][3][1], np.min(design[i, :]), \
              np.max(design[i, :]), variables[i][3][2]

    ###########################################################

    results_base, _ = os.path.splitext(ref_design_name)
    sample_fn = "%s_results.csv" % results_base
    if not os.path.exists(sample_fn) or args.recompute:
        print "Writing results to", sample_fn

        fn = model_run.__dict__[exp_dict['function_name']]

        ## Analytical model
        print "Detailed model..."
        if args.parallel:
            dv['z_func'] = z_func
            dv['fn'] = fn
            fn_par = lambda z : fn(*z)
            dv['fn_par'] = fn_par

            cwd = os.getcwd()
            dv.execute("import sys; sys.path.append('%s'); import model_run" % cwd)

            view = client.load_balanced_view()
            results = view.map_async(fn_par, design.T, ordered=True)
            results.wait_interactive()
        else:
            results = [fn(*z) for z in design.T]
        results = np.array(results[:])

        results_df = pd.DataFrame(results, columns=["Smax_parcel", 
                                                    "Neq_parcel",
                                                    "Nkn_parcel"])
        print "done."

        ## Chaos expansions
        # Load into memory
        results_dict = pickle.load(open("%s_results.dict" % args.exp_name, 'r'))
        n_runs = len(results_dict.keys())
        for run_name, folder in results_dict.iteritems():        
            all_pces = pickle.load(open(os.path.join("save", folder, "pce.p"), 'rb'))
            results_dict[run_name] = { "%s" % r: all_pces[r] for r in responses }

        print "Chaos Expansions..."
        for i, run_name in enumerate(sorted(results_dict.keys())):
            run_dict = results_dict[run_name]
            print "   ", run_name

            for r in responses:
                print "      ", r
                pce = run_dict[r]

                if args.parallel:
                    dv['pce'] = pce
                    pce_results = view.map_async(pce, z_design.T, ordered=True)
                    pce_results.wait_interactive()
                else:
                    pce_results = pce(z_design)
                pce_results = np.array(pce_results)
                results_df['%s_%s' % (r, run_name)] = pce_results

                ## If this is an Smax response, go ahead and compute 
                ## a derived Nact (Nderiv)
                if r == 'Smax':
                    fn_nderiv = lambda z, smax: fn(*z, fn_toggle=smax)
                    zipped = zip(design.T, pce_results)
                    nderiv = np.array([fn_nderiv(z, z_func(smax)) \
                                      for z, smax in zipped])
                    results_df['Nderiv_%s' % run_name] = nderiv[:, -1]

        print "done."

        ## Parameterizations
        if args.params:
            print "Parameterizations"
            fn_param = lambda z: fn(*z, fn_toggle=True)
            if args.parallel:
                dv['fn_param'] = fn_param
                param_results = view.map_async(fn_param, design.T, ordered=True)
                param_results.wait_interactive()
            else:
                param_results = [fn_param(z) for z in design.T]
            param_results = np.array(param_results[:])

            results_df["Smax_ARG"] = param_results[:, 0]
            results_df["Neq_ARG"] = param_results[:, 1]
            results_df["Smax_MBN"] = param_results[:, 2]
            results_df["Neq_MBN"] = param_results[:, 3]
            print "done."

        ## Save final results
        print "Saving to disk at %s" % sample_fn
        results_df.to_csv(sample_fn)

