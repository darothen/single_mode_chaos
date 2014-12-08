"""
Hand-crafted wrapper of PCM results derived using DAKOTA. 

Version: 11/18/2014

    1) The `__main__` portion of this script will process the PCE derivation output
       and create a binary file with the coefficient tables. It'll have to be ported
       with the parameterization script.

    2) Only the listed runs and orders will be saved and made available.

    3) Modified to only extract the Smax expansions

"""

import h5py
import numpy as np
import os
import pandas as pd
import pickle

from parcel_model.activation import _unpack_aerosols, lognormal_activation

from scipy.special import orthogonal as op

from numba import jit

DEBUG = False
BLOCK = True # Block from the 'main' program being run
TEST  = True

STORE_FN = "pcm_param.h5"

RUNS = {
    "SM_OLS": [2, 3, 4, ],
}

RUNS_SAVE = {}

if BLOCK:
    ## Load data into memory
    store = h5py.File(STORE_FN, "r")
    for pcm, levels in RUNS.iteritems():
        save_dict = {}
        for level in levels:
            level = "expansion_order_"+str(level)
            ## Access the necessary PCM parameters from the stored data
            group = store["%s/%s" % (pcm, level)]
            coeffs = np.array(group['coeffs'][:])
            orders = np.array(group['orders'][:])
            max_orders = np.array(group['max_orders'][:])
            save_dict[level] = {'coeffs':coeffs, 'max_orders':max_orders, 'orders':orders}
        RUNS_SAVE[pcm] = save_dict
    store.close()


def uni_to_uni(x, ai, bi, af=-1., bf=1.):    
    """ Transform a uniform random variable to one with another set
    of lower/upper bounds.

    Parameters
    ----------
    x : float
        coordinate in original uniform variable space
    ai, bi : float
        lower and upper bounds of the original uniform distribution
    af, bf : float
        lower and upper bounds of the destination uniform distribution

    Returns
    -------
    float, random variable in new uniform distribution
    """
    if x < ai: return af
    if x > bi: return bf
    return ((bf-af)/(bi-ai))*(x - ai) + af

def project(x):
    """
    Takes a vector with the parameters used in the PCE and projects them to
    the vector space spanned by the orthogonal polynomial basis.

    Parameters
    ----------
    x : [lnN, mu, sigma, kappa, lnV, T, P]

    Returns
    -------
    z : x projected into the vector space spanned by Legendre polynomials

    .. note:: see file docstring for more information on the vector space
    """

    bnds = [[1.0, 4.0],
            [-3.0, -1.0],
            [1.2, 3.0],
            [0.0, 1.2],
            [-2.0, 1.0],
            [240., 310.], 
            [50000., 105000.],
            [0.1, 1.0],]

    y = []
    for i, (xi, (ai, bi)) in enumerate(zip(x, bnds)):
        y.append(uni_to_uni(xi, ai, bi))
    return np.array(y)

@jit("f8(i8[:,:],f8[:],f8[:,:])")
def poly_eval(orders, coeffs, zin):
    nterms = orders.shape[0]
    nvars = orders.shape[1]

    evaluation = 0.
    for row in xrange(nterms):
        row_eval = 1.
        for col in xrange(nvars):
            #print (row+1, col+1, orders[row, col], zin[col, orders[row, col]])
            row_eval *= zin[col, orders[row, col]]
        #print (row+1, row_eval)
        evaluation += row_eval*coeffs[row]

    return evaluation

def param(V, T, P, aerosols=[], accom=0.1,
          mus=[], sigmas=[], Ns=[], kappas=[],
          pcm=RUNS.keys()[0], level="expansion_order_4"):

    if aerosols:
        d = _unpack_aerosols(aerosols)
        mu    = d['mus'][0]
        sigma = d['sigmas'][0]
        kappa = d['kappas'][0]
        N     = d['Ns'][0]
    else:
        ## Assert that the aerosol was already decomposed into component vars
        assert mus
        mu = mus[0]
        assert sigmas
        sigma = sigmas[0]
        assert Ns
        N = Ns[0]
        assert kappas
        kappa = kappas[0]

    logN  = np.log10(N)
    logmu = np.log10(mu)
    logV  = np.log10(V)

    ## Project into orthogonal polynomial space
    z = project([logN, logmu, sigma, kappa, logV, T, P, accom])

    coeffs = RUNS_SAVE[pcm][level]['coeffs'] 
    orders = RUNS_SAVE[pcm][level]['orders'] 
    max_orders = RUNS_SAVE[pcm][level]['max_orders'] 

    nterms = len(coeffs)
    nvars = len(z)

    ## Step 1) Evaluate the necessary polynomial orders
    poly_evals = []
    for zi, order in zip(z, max_orders):
        poly_evals.append([op.eval_legendre(n, zi) for n in xrange(order+1)])
    poly_evals = np.array(poly_evals)

    evaluation = poly_eval(orders, coeffs, poly_evals)

    ## Step 3) This gives us log10(Smax), so now we return Smax
    Smax = 10.**(evaluation)

    n_act, act_frac = lognormal_activation(Smax, mu*1e-6, sigma, N, kappa, T=T)

    return Smax, n_act, act_frac

if __name__ == "__main__":

    if TEST:
        import parcel_model as pm
        V = 0.2
        T = 293.
        P = 85000.

        aerosols = [pm.AerosolSpecies('test', pm.Lognorm(mu=0.004, sigma=2.0, N=1.),
                                      kappa=0.507, bins=200), ]

        #for V in np.logspace(-1, 1, 10):
        for V in [0.2]: 

            x = [np.log10(850.), np.log10(0.05), 2.0, 0.507, np.log10(V), T, P]
            print "V", V
            print "PCM", param(V, T, P, aerosols, accom=0.1) 
            print "ARG", pm.arg2000(V, T, P, aerosols, accom=0.1)
            print "MBN", pm.mbn2014(V, T, P, aerosols, accom=0.1)
            print "--"*30

    if BLOCK:
        print "Locked out"
        import sys; sys.exit()

    store = h5py.File(STORE_FN, "w")

    for run, orders in RUNS.iteritems():
        print run

        run_group = store.create_group(run)

        fn_results = run+"_results.dict"
        fn_exp = run+"_exp.dict"

        results = pickle.load(open(fn_results, 'r'))
        exp = pickle.load(open(fn_exp, 'r'))

        for level, folder in results.iteritems():
            print "   ", level
            level_val = int(level.split("_")[-1])

            path_to_output = os.path.join("save", folder)
            pce = pickle.load(open(os.path.join(path_to_output, "pce.p"), 'r'))['Smax']

            level_group = run_group.create_group(level)
            coeffs = level_group.create_dataset("coeffs", pce.A.shape, dtype=pce.A.dtype)
            coeffs[:] = pce.A[:]

            orders = level_group.create_dataset("orders", pce.orders.shape,
                                                dtype=pce.orders.dtype)
            orders[:] = pce.orders[:]

            max_orders = level_group.create_dataset("max_orders", pce.max_orders.shape, 
                                                    dtype=pce.max_orders.dtype)
            max_orders[:] = np.max(pce.orders, axis=0)

            stacked = np.hstack([coeffs[:][:, np.newaxis], orders[:]])
            np.savetxt("%s_%d.ascii" % (run, level_val), stacked)

    store.close()
