"""
Hand-crafted wrapper of PCM results derived using DAKOTA. 

Version: 9/2/2014

## The set of variables about which to perform the PCE
## symbol, name, [prior, *parameters], meta(dictionary)
    ## ~10 - ~10000
    ['lnN', 'lnN', 1, ['uniform', 2.3, 9.3], 2.3, {}],
    ['mu', 'mu', 2, ['uniform', 0., 0.8], 0., {}],
    ['sigma', 'sigma', 3, ['uniform', 1.2, 3.0], 1.2, {}],
    ['kappa', 'kappa', 4, ['uniform', 0.0, 1.2], 0., {}],
    ['lnV', 'lnV', 5, ['uniform', -1.61, 2.5], -1.61, {}],
    ['T', 'T', 6, ['uniform', 240., 300.], 200., {}],
    ['P', 'P', 7, ['uniform', 50000., 105000.], 50000., {}],

## All modes are uniform

Notes:

    1) The `__main__` portion of this script will process the PCE derivation output
       and create a binary file with the coefficient tables. It'll have to be ported
       with the parameterization script.

    2) Only the listed runs and orders will be saved and made available.
"""

import h5py
import numpy as np
import os
import pandas as pd
import pickle

from parcel_model.activation import activate_lognormal_mode, _unpack_aerosols

from scipy.special import orthogonal as op

DEBUG = False
BLOCK = False # Block from the 'main' program being run

STORE_FN = "pcm_param.h5"

RUNS = {
    "pcm_ols_parcel": [2, 3, 4, 5],
    "pcm_lars_parcel": [2, 3, 4, 5],
    "pcm_lasso_parcel": [2, 3, 4, 5],
}

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

    bnds = [[2.3, 9.3],
            [0., 0.8],
            [1.2, 3.0],
            [0.0, 1.2],
            [-1.61, 2.5],
            [240., 300.], 
            [50000., 105000.]]

    y = []
    for i, (xi, (ai, bi)) in enumerate(zip(x, bnds)):
        y.append(uni_to_uni(xi, ai, bi))
    return np.array(y)

def param(V, T, P, aerosols=[],
          mus=[], sigmas=[], Ns=[], kappas=[],
          pcm=RUNS.keys()[0], level="expansion_order_3"):

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

    lnN = np.log(N)
    lnV = np.log(V)

    ## Project into orthogonal polynomial space
    z = project([lnN, mu, sigma, kappa, lnV, T, P])

    ## Access the necessary PCM parameters from the stored data
    store = h5py.File(STORE_FN, "r")
    group = store["%s/%s" % (pcm, level)]
    coeffs = group['coeffs'][:]
    orders = group['orders'][:]
    max_orders = group['max_orders'][:]
    store.close()

    nterms = len(coeffs)
    nvars = len(z)

    if DEBUG: 
        print "max_orders", max_orders
        print "nterms", nterms
        print "nvars", nvars

    ## Step 1) Evaluate the necessary polynomial orders
    poly_evals = []
    for zi, order in zip(z, max_orders):
        poly_evals.append([op.eval_legendre(n, zi) for n in xrange(order+1)])
    if DEBUG: print poly_evals

    ## Step 2) Loop through arrays and evaluate polynomial
    evaluation = 0.
    for i in xrange(nterms):
        row = orders[i]
        if DEBUG: print i, coeffs[i], orders[i]

        row_eval = 1.
        for n in xrange(nvars):
            order = row[n]
            p_zi = poly_evals[n][order]
            if DEBUG: print "   ", n, p_zi
            row_eval *= p_zi
        if DEBUG: print "RUNNING EVAL -> ", evaluation,
        evaluation += coeffs[i]*row_eval
        if DEBUG: print evaluation

    ## Step 3) This gives us ln(Smax), so now we return Smax
    Smax = np.exp(evaluation)

    n_act, act_frac = activate_lognormal_mode(Smax, mu*1e-6, sigma, N, kappa, T=T)

    return Smax, n_act, act_frac

if __name__ == "__main__":

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

            path_to_output = os.path.join("save", folder)
            pce = pickle.load(open(os.path.join(path_to_output, "pce.p"), 'r'))

            level_group = run_group.create_group(level)
            coeffs = level_group.create_dataset("coeffs", pce.A.shape, dtype=pce.A.dtype)
            coeffs[:] = pce.A[:]

            orders = level_group.create_dataset("orders", pce.orders.shape,
                                                dtype=pce.orders.dtype)
            orders[:] = pce.orders[:]

            max_orders = level_group.create_dataset("max_orders", pce.max_orders.shape, 
                                                    dtype=pce.max_orders.dtype)
            max_orders[:] = np.max(pce.orders, axis=0)

    store.close()
