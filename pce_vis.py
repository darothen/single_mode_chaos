import os, pickle
import numpy as np
import pandas as pd

from statsmodels.distributions import ECDF

import model_run

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import sklearn.metrics as skm
from statsmodels.graphics.gofplots import qqplot

import argparse

parser = argparse.ArgumentParser(description="Visualize some PCE evaluation metrics")
parser.add_argument("exp_name", type=str, help="name of PCM experiment to reference")
parser.add_argument("-r", "--reference", type=str, required=True,
                    help="Reference sample design to plot")
parser.add_argument("-i", "--interact", action="store_true",
                    help="Enable interactive mode (draw plots to screen)")
parser.add_argument("--params", action='store_true', 
                    help="Include plots for ARG/MBN parameterizations")
parser.add_argument("--plots", nargs="+", type=str,
                    help="Control which plots to produce")

def compute_stats(obs, act):
    mae = skm.mean_absolute_error(act, obs)
    r2  = skm.r2_score(act, obs)
    rmse  = np.sqrt(np.sum((obs-act)**2.)/len(act))
    nrmse = rmse/np.sqrt(np.sum((act**2.)/len(act)))
    

    rel_err = 100.*(obs - act)/act
    ## Mask egregiously high values (1000% error) which screw up the spread
    rel_err = rel_err[np.abs(rel_err) <= 1000.]
    mre = np.mean(rel_err)
    mre_std = np.std(rel_err)

    stats = {
        'mae': mae, 'r2': r2, 'rmse': rmse, 'nrmse': nrmse,
        'mre': mre, 'mre_std': mre_std,
    }

    return stats

def plot_dists_base(param, parcel, var_name, param_name, lims, exp_name = '', ref_name='', savefig=False, **kwargs):
    pct_levs = np.linspace(0., 100., 11)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    plt.subplots_adjust(wspace=0.25, bottom=0.15)
    ax_cdf, ax_pdf = axs

    ## Compute empirical CDFs
    param_cdf = ECDF(param.ravel())
    param_percentiles = [np.percentile(param.ravel(), x) for x in pct_levs]
    parcel_cdf = ECDF(parcel.ravel())
    parcel_percentiles = [np.percentile(parcel.ravel(), x) for x in pct_levs]

    ax_cdf.plot(parcel_percentiles, pct_levs/100., color='k', lw=5,
                   label="parcel model")
    ax_cdf.plot(param_percentiles, pct_levs/100., "--o", ms=8, label=param_name)

    ax_cdf.legend(loc='best')
    ax_cdf.set_xlim(*lims)
    ax_cdf.set_xlabel(var_name)
    ax_cdf.set_ylim(0, 1)
    ax_cdf.set_ylabel("Cumulative Probability")

    ## PDFs
    ax_pdf = sns.distplot(parcel, hist=False,
                          color='k', label="parcel model", ax=ax_pdf,
                          kde_kws={'lw': 5})
    ax_pdf = sns.distplot(param, hist=False, kde_kws={'linestyle': 'dashed'},
                          label=param_name, ax=ax_pdf)
    ax_pdf.set_xlim(*lims)
    ax_pdf.set_xlabel(var_name)
    ax_pdf.set_ylabel("Probability Density")
    ax_pdf.legend(loc="best")

    if savefig:
        var_fn = fn_out_fix(var_name)
        plt.savefig(os.path.join(plot_dir, "%s_%s_%s_%s_cdfs.pdf" % (ref_name, exp_name, var_fn, param_name)))

    return ax_cdf, ax_pdf

def plot_one_one_base(param, parcel, var_name, param_name, lims, coloring='k', loglog=False, ax=None, error_pct=0.5, exp_name='', ref_name='', savefig=False, **kwargs):
    if not ax:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)

    ## Mask for infinites
    mask = np.isfinite(parcel)
    parcel = parcel[mask]
    param  = param[mask]
    if isinstance(coloring, (pd.Series, np.ndarray, list)):
        coloring = coloring[mask]


    ax.scatter(parcel, param, marker='.', s=30, c=coloring,
               edgecolor='none', alpha=0.8, cmap=plt.get_cmap("OrRd"), **kwargs)
    oo = np.linspace(lims[0], lims[1], 100)
    ax.plot(oo, oo, color='grey', lw=3)
    ax.plot(oo, oo*error_pct, color='k', lw=1, alpha=0.8)
    ax.plot(oo, oo*(1.+error_pct), color='k', lw=1, alpha=0.8)

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    if loglog:
        ax.loglog()

    ax.set_xlabel("%s, parcel model" % var_name)
    ax.set_ylabel("%s, %s" % (var_name, param_name))

    stats = compute_stats(param, parcel)
    ax.text(0.05, 0.775, stat_label.format(**stats),
            transform=ax.transAxes, fontsize=12)

    if savefig:
        var_fn = fn_out_fix(var_name)
        plt.savefig(os.path.join(plot_dir, "%s_%s_%s_%s_oneone.pdf" % (ref_name, exp_name, var_fn, param_name)))

    return ax

###########################################################
## Plot aesthetics and configuration

sns.set(style="darkgrid", context="talk", palette="Dark2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def fn_out_fix(s, chars="_${}"):
    for char in chars:
        s = s.replace(char, "")
    return s

fn_fix = lambda s: s.replace("_", "\_")
def param_name(s):
    if s in ["ARG", "MBN"]:
        return s 
    else:
        bits = s.split("_")
        order = int(bits[-1])
        return "PCE order %d" % order

stat_label = "RMSE: {rmse:1.2f} ({nrmse:1.2f})\n" + \
             " MAE: {mae:1.2f}\n"  + \
             " R$^2$: {r2:1.2f}\n" + \
             " MRE: {mre:2.2f}$\%$ ({mre_std:2.2f}$\%$)"

z_func = lambda z: 10.**z

plot_dir = "figs/"

N_lims = 1, 5
S_lims = -5, -1

###########################################################

if __name__ == "__main__":

    args = parser.parse_args()
    print "Creating plots for %s" % args.exp_name
    print "   Sample data:", args.reference
    if args.interact:
        print "   Interactive mode"
        plt.ion()
    else:
        plt.ioff()
    if args.params:
        print "   Plotting ARG/MBN figures"
    if not args.plots:
        ALL_PLOTS = True
        print "   Plotting all plots"
    else:
        ALL_PLOTS = False
        print "   Plotting", ", ".join(args.plots)

    print "\n"

    ## Over-write plotting functions with exp_name
    def plot_one_one(*v, **kw):
        return plot_one_one_base(*v, exp_name=args.exp_name, ref_name=args.reference,
                                 savefig=(not args.interact), **kw)
    def plot_dists(*v, **kw):
        return plot_dists_base(*v, exp_name=args.exp_name, ref_name=args.reference,
                               savefig=(not args.interact),**kw)

    ## Unload the configured experiment
    exp_dict = pickle.load(open("%s_exp.dict" % args.exp_name, 'r'))
    results_dict = pickle.load(open("%s_results.dict" % args.exp_name, 'r'))
    design_df = pd.read_csv("%s.csv" % args.reference, index_col=0)
    results_df = pd.read_csv("%s_results.csv" % args.reference, index_col=0)

    pce_keys = results_dict.keys()
    if args.params: 
        param_keys = pce_keys + ["ARG", "MBN"]
    else:
        param_keys = pce_keys

    all_plots = {}


    ###########################################################
    ## SET 1) one-one plots
    if ALL_PLOTS or "oneone" in args.plots:
        print "One-one plots..."

        oo_kwargs = { 
            #'coloring': 10.**design_df['logV'], 
        }
        
        # a) log10(Smax)
        var_key  = "Smax"
        var_name = r"log10(S_max)"
        lims     = [-4, -1]
        parcel = results_df['%s_parcel' % var_key]
        print var_name
        for key in param_keys:
            print "   ", key
            ax = plot_one_one(results_df['%s_%s' % (var_key, key)], parcel,
                              var_name, param_name(key), lims, **oo_kwargs)
            fig = ax.get_figure()
            all_plots[var_name, key] = (fig, ax)
            if args.interact: plt.show()

        # b) Smax
        var_key  = "Smax"
        var_name = r"S_max"
        lims     = [1e-4, 5e-1]
        parcel = z_func(results_df['%s_parcel' % var_key])
        print var_name
        for key in param_keys:
            print "   ", key
            ax = plot_one_one(z_func(results_df['%s_%s' % (var_key, key)]), parcel,
                              var_name, param_name(key), lims, loglog=True, **oo_kwargs)
            fig = ax.get_figure()
            all_plots[var_name, key] = (fig, ax)
            if args.interact: plt.show()

        # c) Neq
        var_key  = "Neq"
        var_name = r"log10(N_eq)"
        lims     = [1, 4]
        parcel = results_df['%s_parcel' % var_key]
        print var_name
        for key in param_keys:
            print "   ", key
            ax = plot_one_one(results_df['%s_%s' % (var_key, key)], parcel,
                              var_name, param_name(key), lims, **oo_kwargs)
            fig = ax.get_figure()
            all_plots[var_name, key] = (fig, ax)
            if args.interact: plt.show()

        # e) Nderiv
        var_key  = "Nderiv"
        var_name = r"log10(N_d)"
        lims     = [1, 4]
        parcel = results_df['Neq_parcel']
        print var_name
        for key in pce_keys:
            print "   ", key
            ax = plot_one_one(results_df['%s_%s' % (var_key, key)], parcel,
                              var_name, param_name(key), lims, **oo_kwargs)
            fig = ax.get_figure()
            all_plots[var_name, key] = (fig, ax)
            if args.interact: plt.show()

        # c) Nderiv
        var_key  = "Nderiv"
        var_name = r"N$_{d}$"
        lims     = map(lambda x: 10.**x, N_lims)
        parcel = 10.**(results_df['Neq_parcel'])
        print var_name
        for key in pce_keys:
            print "   ", key
            ax  = plot_one_one(10.**(results_df['%s_%s' % (var_key, key)]), parcel,
                               var_name, param_name(key), lims, loglog=True, **oo_kwargs)
            fig = ax.get_figure()
            all_plots[var_name, key] = (fig, ax)
            if args.interact: plt.show()


    ###########################################################
    ## SET 2) one-one plots
    if ALL_PLOTS or "pdf" in args.plots:
        print "CDFs/PDFs..."

        pdf_kwargs = { }

        # a) log10(Smax)
        var_key  = "Smax"
        var_name = r"log10(S_max)"
        lims     = [-5, 0]
        parcel = results_df['%s_parcel' % var_key]
        print var_name
        for key in param_keys:
            print "   ", key
            axs  = plot_dists(results_df['%s_%s' % (var_key, key)], parcel,
                            var_name, param_name(key), lims, **pdf_kwargs)
            fig = axs[0].get_figure()
            all_plots[var_name, key] = (fig, axs)
            if args.interact: plt.show()

        # b) log10(Neq)
        var_key  = "Neq"
        var_name = r"log10(N_eq)"
        lims     = [0, 4]
        parcel = results_df['%s_parcel' % var_key]
        print var_name
        for key in param_keys:
            print "   ", key
            axs  = plot_dists(results_df['%s_%s' % (var_key, key)], parcel,
                            var_name, param_name(key), lims, **pdf_kwargs)
            fig = axs[0].get_figure()
            all_plots[var_name, key] = (fig, axs)
            if args.interact: plt.show()

        # c) log10(Nderiv)
        var_key  = "Nderiv"
        var_name = r"log10(N_d)"
        lims     = [0, 4]
        parcel = results_df['Neq_parcel']
        print var_name
        for key in pce_keys:
            print "   ", key
            axs  = plot_dists(results_df['%s_%s' % (var_key, key)], parcel,
                            var_name, param_name(key), lims, **pdf_kwargs)
            fig = axs[0].get_figure()
            all_plots[var_name, key] = (fig, axs)
            if args.interact: plt.show()

    ## Clean up
    #plt.close('all')


'''
###############################################################
## Loop over all the experiments
run_names = sorted(results_dict.keys())
for i, run_name in enumerate(run_names):
    #if run_name != "expansion_order_4": continue
    run_dict = results_dict[run_name]
        
    n_terms.append(run_dict['pce'].nterms)

    pce_lhs_results = run_dict['pce'](z_design)
    np.save("%s_%s_pce_smax.npy" % (exp_name, run_name), pce_lhs_results)
    pce_lhs_results_positive = np.ma.masked_less_equal(pce_lhs_results, 0)

    pce_ecdf = ECDF(pce_lhs_results.ravel())
    pce_percentiles = [np.percentile(pce_lhs_results.ravel(), x) for x in pct_levs]

    n_pos = pce_lhs_results_positive.count()
    n_tot = len(pce_lhs_results)

    if use_log:
        ln_pce = np.log(pce_lhs_results[pce_lhs_results > 0.])
        ln_ana = np.log(lhs_results[pce_lhs_results > 0.])
        rmse = np.sqrt(np.sum((ln_pce - ln_ana)**2.)/n_good)
        RMSEs.append(rmse)
    else:
        rmse = np.sqrt(np.sum(pce_lhs_results - lhs_results)**2.)/n_tot
        RMSEs.append(rmse)

    mae = skm.mean_absolute_error(lhs_results, pce_lhs_results)
    r2  = skm.r2_score(lhs_results, pce_lhs_results)
    rel_err = 100.*(pce_lhs_results - lhs_results)/lhs_results
    mre = np.mean(rel_err)
    mre_std = np.std(rel_err)

    print run_name
    print """
response fn eval
----------------
    (min, max): {min}, {max}
    # of terms: {n_terms}
          RMSE: {rmse}
           MAE: {mae}
           R^2: {r2}
  mean rel err: {mre:2.2f}% ({mre_std:2.2f}%)
""".format(n_terms=run_dict['pce'].nterms, rmse=rmse, mae=mae, r2=r2,
           min=pce_lhs_results.min(), max=pce_lhs_results.max(), 
           mre=mre, mre_std=mre_std)

    #if not run_name == "sparse_grid_level_6": continue
    
    if PLOT_1 or PLOT_2 or PLOT_3 or PLOT_5:
        fig, [[ax_cdf, ax_pdf, ax_bins], 
              [ax_oo, ax_nact, ax_af], 
              [ax_oo_bin, ax_nact_bin, ax_af_bin]] = \
            plt.subplots(3, 3, num=run_name, figsize=(15, 12))

    ###########################################################
    ## 1) CDF/PDF for true vs each PCE generated in the experiment
    print "Figure 1 - CDF/PDF for each PCE from experiment"
    if PLOT_1:
        ax_cdf.plot(lhs_percentiles, pct_levs/100., color='k', lw=5,
                       label="analytical")
        ax_cdf.plot(pce_percentiles, pct_levs/100., "-", ms=2, label=fn_fix(run_name))

        ax_cdf.legend(loc='lower right')
        ax_cdf.set_xlim(res_min, res_max)
        if use_log:
            ax_cdf.semilogx()
        ax_cdf.set_xlabel('Response')
        ax_cdf.set_ylim(0, 1)
        ax_cdf.set_ylabel("Cumulative Probability")

        if use_log:
            mask = pce_lhs_results > 0
            ax_pdf = sns.distplot(lhs_results[mask], hist=False,
                                  color='k', label="analytical", ax=ax_pdf,
                                  kde_kws={'lw': 5})
            ax_pdf = sns.distplot(pce_lhs_results[mask], hist=False,
                                  color='r', label=fn_fix(run_name), ax=ax_pdf)
        else:
            ax_pdf = sns.distplot(lhs_results, hist=False,
                                  color='k', label="analytical", ax=ax_pdf,
                                  kde_kws={'lw': 5})
            ax_pdf = sns.distplot(pce_lhs_results, hist=False,
                                  color='r', label=fn_fix(run_name), ax=ax_pdf)
        ax_pdf.set_xlabel('Response')
        ax_pdf.set_ylabel("Probability Density")
        ax_pdf.legend(loc="upper left")

    else: print "...skipping"

    ###########################################################
    ## 2) one-one plot for analytical vs each PCE
    print "Figure 2 - One-to-one plots for each PCE from experiment"
    if PLOT_2:

        if use_log:
            ss = np.logspace(np.log10(res_min), np.log10(res_max), 100)
            mask = pce_lhs_results > 0
            ax_oo.scatter(lhs_results[mask], pce_lhs_results[mask], marker='.', s=12, 
                color='k', alpha=0.5, label=run_name, edgecolor='none')
            ax_oo.loglog()
            ax_oo.plot(ss, ss, color='grey', lw=3)#, zorder=100)
        else:
            ss = np.linspace(res_min, res_max, 100)
            ax_oo.scatter(lhs_results, pce_lhs_results, color='k', alpha=0.5, 
                marker='.', s=12, label=run_name,
                edgecolor='none')
            #bins = np.linspace(res_min, res_max, 21)
            #ax_oo.hist2d(lhs_results, pce_lhs_results, bins=bins,
            #             norm=LogNorm())
            ax_oo.plot(ss, ss, color='grey', lw=3)#, zorder=100)
            ax_oo.plot(ss, ss*.5, color='k', lw=1, alpha=0.8)#, zorder=100)
            ax_oo.plot(ss, ss*2., color='k', lw=1, alpha=0.8)#, zorder=100)
        ax_oo.set_xlim(res_min, res_max)
        ax_oo.set_ylim(res_min, res_max)
        ax_oo.set_xlabel("Analytical")
        ax_oo.set_ylabel("Emulator")
        ax_oo.set_title("Modeled response function", loc='left')
        ax_oo.text(0.05, 0.775, stat_label.format(rmse=rmse, mae=mae, 
                                                  r2=r2, mre=mre, mre_std=mre_std), 
                   transform=ax_oo.transAxes)

    else: print "...skipping"

    ###########################################################
    ## 2) one-one plot for analytical vs each PCE, but number activated
    print "Figure 3 - One-to-one plots for Nacts from experiment"
    if PLOT_3:

        fn = model_run.__dict__[exp_dict['function_name']]
        zipped = zip(design.T, pce_lhs_results)

        fn_nact = lambda z, smax : fn(*z, fn_toggle=smax)
        pce_nacts = np.array([fn_nact(z, z_func(smax)) for z, smax in zipped])

        if READ_CACHED_NACTS:
            pce_nacts = np.load("%s_%s_pce_nact.npy" % (exp_name, run_name))
        else:
            pce_nacts = np.array([fn_nact(z, 10.**(smax)) for z, smax in zipped])
            np.save("%s_%s_pce_nact.npy" % (exp_name, run_name), pce_nacts)

        ss = np.linspace(10, 40000, 100)

        ax_nact.set_xlabel("Analytical")
        ax_nact.set_ylabel("Emulator")
        ax_nact.set_title("Computed CDNC", loc='left')

        mask = lhs_results > np.log10(0.01/100) # 0.05 % Smax

        ax_nact.scatter(lhs_nacts[mask], pce_nacts[mask], 
                        c=lhs_results[mask], marker='.', s=12, label=run_name,
                        edgecolor='none', cmap=plt.get_cmap("OrRd"))
        #bins = np.logspace(1, 4.4, 21)
        #ax_nact.hist2d(lhs_nacts, pce_nacts)
        ax_nact.plot(ss, ss, color='grey', lw=3)#, zorder=100)
        ax_nact.plot(ss, ss*.5, color='k', lw=1, alpha=0.8)#, zorder=100)
        ax_nact.plot(ss, ss*2., color='k', lw=1, alpha=0.8)#, zorder=100)
        ax_nact.semilogx()
        ax_nact.semilogy()
        ax_nact.set_xlim(10, 1000)
        ax_nact.set_ylim(10, 1000)

        rmse = np.sqrt(np.sum(pce_nacts - lhs_nacts)**2.)/n_tot
        mae  = skm.mean_absolute_error(lhs_nacts, pce_nacts)
        r2   = skm.r2_score(lhs_nacts[mask], pce_nacts[mask]) 
        rel_err = 100*(pce_nacts - lhs_nacts)/lhs_nacts
        rel_err = np.ma.masked_greater(rel_err, 10)

        mre = np.mean(rel_err[mask])
        mre_std = np.std(rel_err[mask])

        print """
diag CDNC eval
----------------
          RMSE: {rmse}
           MAE: {mae}
           R^2: {r2}           
  mean rel err: {mre:2.2f}% ({mre_std:2.2f}%)
""".format(n_terms=run_dict['pce'].nterms, rmse=rmse, mae=mae, r2=r2,
           mre=mre, mre_std=mre_std)

        ax_nact.text(0.05, 0.775, stat_label.format(rmse=rmse, mae=mae, r2=r2,
                                                  mre=mre, mre_std=mre_std), 
                     transform=ax_nact.transAxes)

    else: print "...skipping"

    ###########################################################
    ## 3) one-one plot for analytical vs each PCE, but act frac
    print "Figure 3 - One-to-one plots for act fracs from experiment"
    if PLOT_3:

        lhs_afs = lhs_nacts/Ns
        pce_afs = pce_nacts/Ns

        ss = np.linspace(0, 1, 100)
        ax_af.set_xlim(0, 1)
        ax_af.set_ylim(0, 1)
        ax_af.set_xlabel("Analytical")
        ax_af.set_ylabel("Emulator")
        ax_af.set_title("Computed Act. Fraction", loc='left')

        ax_af.scatter(lhs_afs, pce_afs, 
                      c=lhs_results, marker='.', s=12, label=run_name,
                      edgecolor='none', cmap=plt.get_cmap("OrRd"))
        ax_af.plot(ss, ss, color='grey', lw=3)#, zorder=100)
        ax_af.plot(ss, ss*.5, color='k', lw=1, alpha=0.8)#, zorder=100)
        ax_af.plot(ss, ss*2., color='k', lw=1, alpha=0.8)#, zorder=100)

    else: print "...skipping"

    ##########################################################
    ## 3) binned plot of spread as a fcn of log(Smax)
    print "Figure 5 - binned spread - log(Smax)"
    if PLOT_5:
        plt.rc('text', usetex=False)

        log_smax_bins = np.linspace(res_min, res_max, 12+1)
        bin_names = ["%2.1f - %2.1f\n%2d" % (l, r, i) for i, (l, r) \
                     in enumerate(zip(log_smax_bins[:-1], log_smax_bins[1:])) ]

        bin_numbers = np.digitize(lhs_results, log_smax_bins[1:-1])

        bin_df = pd.DataFrame(data={ "bin_numbers": bin_numbers, 
                                     "ana_lhs": lhs_results, 
                                     "pce_lhs": pce_lhs_results, 
                                     "ana_nacts": lhs_nacts, 
                                     "pce_nacts": pce_nacts,
                                     "ana_afs": lhs_afs,
                                     "pce_afs": pce_afs,
                                     "rel_lhs_diff": (pce_lhs_results - lhs_results)/lhs_results,
                                     "rel_nacts_diff": (pce_nacts - lhs_nacts)/lhs_nacts,
                                     "rel_afs_diff": (pce_afs - lhs_afs)/lhs_afs })

        Smaxes = z_func(bin_df.ana_lhs)*100.
        rel_errs = bin_df.rel_lhs_diff*100.
        mask = (Smaxes > 0.01) & (Vs > 0.2) & (Vs < 10.0)
        print "   Number of vals in rel err sample:", len(rel_errs[mask])
        ax_oo_bin.scatter(Smaxes[mask], rel_errs[mask],
                          edgecolor='None', color='k', alpha=0.5, marker='.')
        ax_oo_bin.semilogx()
        ax_oo_bin.set_ylim(-50, 50)
        ax_oo_bin.set_xlim(1e-2, 5.)
        #ax_oo_bin.set_xticklabels(bin_names, rotation=90)
        #ax_oo_bin.set_xlabel(fn_fix("bin_number"))
        ax_oo_bin.set_xlabel("Smax ($\%$)")
        ax_oo_bin.set_ylabel(fn_fix("Relative Error in Smax ($\%$)"))

        sns.barplot(bin_numbers, ci=None, palette="OrRd", ax=ax_bins)
        ax_bins.set_xlabel(fn_fix("bin_number"))

        plt.rc('text', usetex=True)
    else: print "...skipping"

    ##########################################################
    ## 4) binned plot of spread in nact as a fcn of log(Smax)
    print "Figure 5 - binned spread - Nact"
    if PLOT_5:
        plt.rc('text', usetex=False)

        Smaxes = z_func(bin_df.ana_lhs)*100.
        rel_errs = bin_df.rel_nacts_diff*100.
        mask = (Smaxes > 0.01)  & (Vs > 0.2) & (Vs < 10.0)
        ax_nact_bin.scatter(Smaxes[mask], rel_errs[mask], 
                            edgecolor='None', color='k', alpha=0.5, marker='.')
        ax_nact_bin.semilogx()
        #ax_nact_bin.set_xticklabels(bin_names, rotation=90)
        ax_nact_bin.set_ylim(-150, 150)
        ax_nact_bin.set_xlim(1e-2, 5.)
        #ax_nact_bin.set_xlabel(fn_fix("bin_number"))
        ax_nact_bin.set_xlabel("Smax ($\%$)")
        ax_nact_bin.set_ylabel(fn_fix("Relative Error in CDNC ($\%$)"))

        plt.rc('text', usetex=True)
    else: print "...skipping"

    ##########################################################
    ## 4) binned plot of spread in nact as a fcn of log(Smax)
    print "Figure 5 - binned spread - act fracs"
    if PLOT_5:
        plt.rc('text', usetex=False)

        Smaxes = z_func(bin_df.ana_lhs)*100.
        rel_errs = bin_df.rel_afs_diff*100.
        mask = (Smaxes > 0.01)  & (Vs > 0.2) & (Vs < 10.0)
        ax_af_bin.scatter(Smaxes[mask], rel_errs[mask],
                          edgecolor='none', color='k', alpha=0.5, marker='.')
        ax_af_bin.semilogx()
        #ax_af_bin.set_xticklabels(bin_names, rotation=90)
        ax_af_bin.set_ylim(-150, 150)
        ax_af_bin.set_xlim(1e-2, 5.)
        #ax_af_bin.set_xlabel(fn_fix("bin_number"))
        ax_af_bin.set_xlabel("Smax ($\%$)")
        ax_af_bin.set_ylabel(fn_fix("Relative Error in act. frac. ($\%$)"))

        plt.rc('text', usetex=True)
    else: print "...skipping"

    ###########################################################
    if PLOT_1 or PLOT_2 or PLOT_3 or PLOT_5:
        sns.despine()

        #plt.draw()
        plt.tight_layout()
        plt.pause(0.1)
        #fig.tight_layout()
        #break

        c = raw_input("save/continue?")
        if c == "y":
            plt.savefig(plot_dir+exp_name+"_"+run_name[-1]+".pdf", 
                        transparent=True, bbox_inches='tight')

###########################################################
## 3) RMSE in ln(response)
print "Figure 4 - Error"
if PLOT_4:
    fig = plt.figure(num="error_change", figsize=(5,4))
    ax_err = fig.add_subplot(111)

    err_df = pd.DataFrame(data={'nterms': n_terms[::-1], 
                                'RMSE': RMSEs[::-1]})
    print err_df
    err_df.plot(x='nterms', y="RMSE", ax=ax_err)

    ax_err.scatter(n_terms, RMSEs)
    plt.pause(0.1)

else: print "...skipping"

###########################################################
## 6) Matrix of scatterplots to evaluate LHS distribution
print "Figure 6 - LHS scatterplots"
if PLOT_6:
    plt.rc('text', usetex=False)

    ## Not quite working at the moment due to bug in seaborn
    #fig = plt.figure(num="LHS scatterplot")

    var_names = [v[0] for v in exp_dict['variables']]
    design_fix = design[:, ::50].copy()
    labels = np.random.randint(1, 4, len(design_fix[0, :]))
    for i, v in enumerate(var_names):
        if v.startswith("log"): 
            var_names[i] = v[2:]
            design_fix[i, :] = z_func(design_fix[i, :])

    design_df = pd.DataFrame(design_fix.T, columns=var_names)
    design_df['label'] = labels

    #g = sns.PairGrid(design_df, hue='label', diag_sharey=False, size=1.5)
    #g.map_lower(plt.scatter)
    #g.map_diag(sns.kdeplot, lw=3)

    from itertools import combinations
    var_combos = np.array(list(combinations(var_names, 2)))
    inds = np.random.choice(range(len(var_combos)), 4)
    var_combos = var_combos[inds]

    fig, axes = plt.subplots(2, 2, num="LHS scatterplot")
    for i, ax in enumerate(axes.ravel()):
        var_x, var_y = var_combos[i]
        design_df.plot(var_x, var_y, kind='scatter', ax=ax)
    plt.tight_layout()
    plt.pause(0.1)

    plt.rc('text', usetex=True)
else: print "...skipping"
'''
