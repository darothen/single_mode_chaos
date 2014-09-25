"""
Hand-crafted wrapper of PCM results derived using DAKOTA. 

Version: 9/24/2014

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

    3) By default, will attempt to read in a saved binary file with the coefficient tables
       but will fail gracefully if not available.
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

RUNS_SAVE = {}

## Load data into memory
try:
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
except IOError:
    print "Could not open coefficient file %s, check to see it is available." % STORE_FN

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

@jit("f8(i8[:,:],f8[:],f8[:,:])")
def poly_eval(orders, coeffs, zin):
    """
    Interior, polynomial-evaluation loop of the parameterization.

    Parameters
    ----------
    orders : ndarray of integers
        the nterms x nvars matrix corresponding to the order of each
        legendre basis comprising each monomial in the polynomial
    coeffs : ndarray of floats
        an nterm length vector of the chaos coefficients
    zin : ndarray of floats 
        evaluations of each order orthogonal polynomial for the given
        number of inputs, with shape nvars x max_order, where max_order
        is the highest order term in the expansion.

    Returns
    -------
    evaluation : float
        evaluation of the polynomial

    """

    nterms = orders.shape[0]
    nvars = orders.shape[1]

    evaluation = 0.
    for row in xrange(nterms):
        row_eval = 1.
        for col in xrange(nvars):
            row_eval *= zin[col, orders[row, col]]
        evaluation += row_eval*coeffs[row]

    return evaluation

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

    ## Retrieve the specified chaos expansion polynomial details
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

    ## Step 2) Evaluate the polynomial 
    evaluation = poly_eval(orders, coeffs, poly_evals)

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
