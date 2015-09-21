import os, pickle
import numpy as np
import pandas as pd

import model_run
from dist_tools import design_lhs_exp, map_transfer_fcns

import argparse

parser = argparse.ArgumentParser(description="Conduct sampling study of parcel model and chaos expansion metamodels.")
parser.add_argument("exp_name", type=str, help="name of PCM experiment to reference")
parser.add_argument("-r", "--reference", type=str,
                    help="Reference set of design points to use")
parser.add_argument("--params", action='store_true',
                    help="Include sampling of ARG/MBN parameterizations")
parser.add_argument("--parcel", action='store_true',
                    help="Include sampling of numerical parcel model.")
parser.add_argument("--n", "-n", type=int, default=10000,
                    help="Number of samples to perform; either total for LHS or " + \
                         "limit for CESM")
parser.add_argument("--parallel", action='store_true',
                    help="Perform evaluations in parallel")
parser.add_argument("--recompute", action='store_true',
                    help="Force recompute/overwrite of existing results")
parser.add_argument("--project", action='store_true',
                    help="Re-project logarithmic variables to linear space for sampling")

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
        client = Client(profile='legion')
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
        ref_design_name = args.reference
        stored_design = pd.read_csv(args.reference)
        design = stored_design[[v[0] for v in variables]].values.T
        z_design = np.empty_like(design)
        for i, v in enumerate(variables):
            dist, a, b = v[3]
            z_design[i, :] = maps[i](design[i, :], a, b)

        design = design[:, :args.n]
        z_design = z_design[:, :args.n]
        design_df = stored_design

    else:
        design, z_design = design_lhs_exp(variables, maps, offsets, args.n,
                                          project_linear=args.project)
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
        results_df = pd.DataFrame(index=design_df.index)

        ## Set up some parallel processing arguments if they're going
        ## to be used
        if args.parallel:

            ## Set import path on remote engines
            LOCAL_PATH = os.getcwd()
            REMOTE_PATH = "/net/legion" + LOCAL_PATH

            pth = sys.path
            sys.path.append(LOCAL_PATH)
            sys.path.append(REMOTE_PATH)

            with dv.sync_imports():
                import sys
            dv['sys.path'] = pth

            ## Send vital functions to remote engines
            dv['fn'] = fn
            fn_par = lambda z : fn(*z)
            dv['fn_par'] = fn_par
            dv['exp_dict'] = exp_dict
            dv['z_func'] = z_func

            # Generate load-balancing view
            view = client.load_balanced_view()

        ## Analytical model
        if args.parcel:
            print "Detailed model..."
            if args.parallel:
                results = view.map_async(fn_par, design.T, ordered=True)
                results.wait_interactive()
            else:
                results = []
                for i, z in enumerate(design.T):
                    result = fn(*z)
                    results.append(result)
            results = np.array(results[:])

            #results_df = pd.DataFrame(results, columns=["Smax_parcel",
            #                                            "Neq_parcel",
            #                                            "Nkn_parcel"])
            results_df['Smax_parcel'] = results[:, 0]
            results_df['Neq_parcel'] = results[:, 1]
            results_df['Nkn_parcel'] = results[:, 2]

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
                    print "      Computing Nderiv from PCE Smax"
                    fn_nderiv = lambda z: fn(*z[0], fn_toggle=z_func(z[1]))
                    zipped = zip(design.T, pce_results)

                    if args.parallel:
                        dv['fn_nderiv'] = fn_nderiv
                        dv['z_func'] = z_func
                        nderiv = view.map_async(fn_nderiv, zipped, ordered=True)
                        nderiv.wait_interactive()
                    else:
                        nderiv = np.array([fn_nderiv(z) \
                                          for z in zipped])
                    nderiv = np.array(nderiv)
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
            results_df["Nderiv_ARG"] = param_results[:, 1]
            results_df["Smax_MBN"] = param_results[:, 2]
            results_df["Nderiv_MBN"] = param_results[:, 3]
            print "done."

        ## Save final results
        print "Saving to disk at %s" % sample_fn
        if args.reference:
            results_df.index = stored_design.index
        results_df.to_csv(sample_fn)

