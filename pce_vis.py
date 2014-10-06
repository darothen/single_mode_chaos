import os, pickle
import numpy as np
import pandas as pd

from statsmodels.distributions import ECDF

import model_run
from dist_tools import design_lhs_exp, map_transfer_fcns

import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import sklearn.metrics as skm

sns.set(style="ticks")
sns.set_palette("Set2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fn_fix = lambda s: s.replace("_", "\_")

stat_label = "RMSE: {rmse:1.2f}\n" + \
             " MAE: {mae:1.2f}\n"  + \
             " R$^2$: {r2:1.2f}\n" + \
             " MRE: {mre:2.2f}$\%$ ({mre_std:2.2f}$\%$)"

exp_name = "tri_modal_ols"

plot_dir = "plots/"

res_min = -14
res_max = -2
#res_min = 0
#res_max = 3000
use_log = False

PLOT_1 = True ## PDFs/CDFS
PLOT_2 = True ## Smax one-to-one plot
PLOT_3 = True ## Nact/act-frac one-to-one plot
PLOT_4 = False ## Error vs # terms in poly
PLOT_5 = True ## Binned plots
PLOT_6 = True ## LHS Scatter plots

###########################################################

## Unload the configured experiment
exp_dict = pickle.load(open("%s_exp.dict" % exp_name, 'r'))
results_dict = pickle.load(open("%s_results.dict" % exp_name, 'r'))

## Clean up - remove expansion order 0 (it's obviously never good)
if "expansion_order_0" in results_dict: del results_dict['expansion_order_0']
if "expansion_order_5" in results_dict: del results_dict['expansion_order_5']

n_runs = len(results_dict.keys())
for run_name, folder in results_dict.iteritems():
    pce = pickle.load(open(os.path.join("save", folder, "pce.p"), 'rb'))
    results_dict[run_name] = { 
        'foldername': folder, 
        'pce': pce,
    }

dataset     = np.load("%s_LHS_sample.npz" % exp_name)
design      = dataset['design']
Ns = np.sum(np.exp(design[:3,:]), axis=0).shape
Vs = np.exp(design[4, :])
z_design    = dataset['z_design']
lhs_results = dataset['results']
lhs_nacts   = dataset['Nacts']
print exp_name
print "lhs", lhs_results.min(), lhs_results.max()
print "\n\n\n"

lhs_ecdf = ECDF(lhs_results.ravel())
pct_levs = np.linspace(1, 99, 99)
lhs_percentiles = [np.percentile(lhs_results.ravel(), x) for x in pct_levs]

n_terms = [] 
RMSEs   = []

###############################################################
## Loop over all the experiments
run_names = sorted(results_dict.keys())
for i, run_name in enumerate(run_names):
    #if run_name != "expansion_order_4": continue
    run_dict = results_dict[run_name]
        
    n_terms.append(run_dict['pce'].nterms)

    pce_lhs_results = run_dict['pce'](z_design)
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
            ax_oo.plot(ss, ss, color='grey', lw=3)#, zorder=100)
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
        pce_nacts = np.array([fn_nact(z, np.exp(smax)) for z, smax in zipped])

        ss = np.linspace(0, 30000, 100)
        ax_nact.set_xlim(0, 30000)
        ax_nact.set_ylim(0, 30000)
        ax_nact.set_xlabel("Analytical")
        ax_nact.set_ylabel("Emulator")
        ax_nact.set_title("Computed CDNC", loc='left')

        ax_nact.scatter(lhs_nacts, pce_nacts, 
                        c=lhs_results, marker='.', s=12, label=run_name,
                        edgecolor='none', cmap=plt.get_cmap("OrRd"))
        ax_nact.plot(ss, ss, color='grey', lw=3)#, zorder=100)

        rmse = np.sqrt(np.sum(pce_nacts - lhs_nacts)**2.)/n_tot
        mae  = skm.mean_absolute_error(lhs_nacts, pce_nacts)
        r2   = skm.r2_score(lhs_nacts, pce_nacts) 
        rel_err = 100*(pce_nacts - lhs_nacts)/lhs_nacts
        rel_err = np.ma.masked_greater(rel_err, 10)

        mre = np.mean(rel_err)
        mre_std = np.std(rel_err)

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

        '''
        sns.violinplot(bin_df.rel_lhs_diff*100., bin_numbers, #names=bin_names,
                        color="OrRd",
                       ax=ax_oo_bin)
        '''
        Smaxes = np.exp(bin_df.ana_lhs)*100.
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

        '''
        sns.violinplot(bin_df.rel_nacts_diff*100., bin_df.bin_numbers, #names=bin_names, 
                       color="OrRd",
                       ax=ax_nact_bin)
        '''
        Smaxes = np.exp(bin_df.ana_lhs)*100.
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
        '''
        sns.violinplot(bin_df.rel_afs_diff*100., bin_df.bin_numbers, #names=bin_names, 
                       color="OrRd",
                       ax=ax_af_bin)
        '''
        Smaxes = np.exp(bin_df.ana_lhs)*100.
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
        if v.startswith("ln"): 
            var_names[i] = v[2:]
            design_fix[i, :] = np.exp(design_fix[i, :])

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
