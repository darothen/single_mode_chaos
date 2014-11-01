"""
Collection of utilities for basic statistical distribution
transformations.

"""

import numpy as np
from scipy.special import erfinv

from pyDOE import lhs
from scipy.stats.distributions import norm, uniform

def design_lhs_exp(variables, maps, offsets=None, samples=int(1e4)):
    """ Design an LHS experiment """

    design = lhs(len(variables), samples=samples, criterion='m', iterations=100)
    z_design = np.zeros_like(design)

    print "Computing LHS design..."
    for i, v in enumerate(variables):
        dist, a, b = v[3]

        if v[0].startswith("ln"):
            ## 9/4/2014
            ## This is an experimental correction to re-project the
            ## logarithmic variables into their normal coordinate 
            ## system. It should only effect the sampling, and hopefully
            ## improve it by forcing it to even things out over the 
            ## actually range we care about
            a = np.exp(a)
            b = np.exp(b)
            offsets[i] = np.exp(offsets[i])

        elif v[0].startswith("log"):
            ## 10/26/2014
            ## In accordance with above, but for log10 vars
            a = 10.**a
            b = 10.**b
            offsets[i] = 10.**offsets[i]

        if offsets:
            ## These corrections with "offsets" re-center the interval
            ## so that the left endpoint is 0. I found that if arbitrary
            ## lower/upper limits were used, sometimes the PPF routines
            ## would really mess up in inverting the CDF.
            a, b = a-offsets[i], b-offsets[i]
        if dist == 'uniform':
            design[:, i] = uniform(a, b).ppf(design[:, i])
        elif dist == 'normal':
            design[:, i] = norm(a, b).ppf(design[:, i])
        elif dist == 'loguniform':
            design[:, i] = loguni_ppf(design[:, i], a, b)
        else:
            raise ValueError("no dist defined for %s" % dist)

        if offsets:
            ## Project back in to the correct limits
            design[:, i] += offsets[i]
            a, b = a+offsets[i], b+offsets[i]

        if v[0].startswith("ln"):
            ## 9/4/2014
            ## Second half of correction
            a = np.log(a)
            b = np.log(b)
            design[:, i] = np.log(design[:, i])
        elif v[0].startswith("log"):
            ## 10/26/2014
            a = np.log10(a)
            b = np.log10(b)
            design[:, i] = np.log10(design[:, i])

        z_design[:, i] = maps[i](design[:, i], a, b)
    design = design.T # in x-coords
    z_design = z_design.T

    return design, z_design

def map_transfer_fcns(dists, pce_directive, verbose=False):
    maps = []
    for dist in dists:
        if verbose: print dist,
        if "wiener" in 'pce_directive':
            if dist == 'uniform':
                fcn = uni_to_norm
            elif dist == 'normal':
                fcn = normal_to_standard
            else:
                raise ValueError("wiener + %s not implemented" % dist)
        else:
            if dist == 'uniform':
                fcn = uni_to_uni
            elif dist == "loguniform":
                fcn = loguni_to_uni
            elif dist == "normal":
                fcn = normal_to_standard
            else:
                raise ValueError("askey + %s not implemented" % dist)
        if verbose: print fcn.func_name
        maps.append(fcn)

    return maps
            
################
## UNIFORM -> ??
################

def uni_to_norm(x, a, b):
    """ Transform a uniform random variable to a standard (normal) 
    random variable.

    Parameters
    ----------
    x : float
        coordinate in uniform variable space
    a, b : float
        lower and upper bounds of the uniform distribution

    Returns
    -------
    float, random variable in SRV space
    """
    return np.sqrt(2)*erfinv(2.*(x-a)/(b-a) - 1.0)

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

##################
## NORMAL -> ??
##################
def normal_to_standard(x, mu, sigma):
    """ Transform a normal random variable to a standard random variable.

    Parameters
    ----------
    x : float
        coordinate in original normal variable space
    mu, sigma : float
        mean and std dev of original normal variable

    Returns
    -------
    float, random variable transformed to standard random variable
    """
    return (x - mu)/sigma

##################
## LOGUNIFORM -> ??
##################
def loguni_to_uni(x, alpha, beta, a=-1., b=1.):
    """ Transform a logunifrom random variable to a uniform one.

    Parameters
    ----------
    x : float
        coordinate in original uniform variable space
    alpha, beta : float
        lower and upper bounds of the original loguniform distribution
    a, b : float
        lower and upper bounds of the destination uniform distribution

    Returns
    -------
    float, random variable in new uniform distribution
    """
    C = (a*np.log(beta) - b*np.log(alpha))/(b - a)
    D = (np.log(beta/alpha))/(b - a)
    return (1./D)*(np.log(x) - C)

def loguni_pdf(x, alpha, beta):
    """ PDF function for Loguniform distribution

    Parameters
    ----------
    x : float
        coordinate in loguniform variable space
    alpha, beta : float
        lower/upper bounds of loguniform span

    Returns
    -------
    float, PDF of the specified loguniform distribution evaluated at point

    """
    return 1./x/np.log(beta/alpha)

def loguni_cdf(x, alpha, beta):
    """ CDF function for Loguniform distribution

    Parameters
    ----------
    x : float
        coordinate in loguniform variable space
    alpha, beta : float
        lower/upper bounds of loguniform span

    Returns
    -------
    float, CDF of the specified loguniform distribution, evaluated
        from left-to-right over the real number line

    """
    return np.log(x/alpha)/np.log(beta/alpha)

def loguni_ppf(q, alpha, beta):
    """ PPF function for Loguniform distribution

    Parameters
    ----------
    x : float
        coordinate in loguniform variable space
    alpha, beta : float
        lower/upper bounds of loguniform span

    Returns
    -------
    float, PPF of the specified loguniform distribution,

    """
    return np.exp(q*np.log(beta/alpha) + np.log(alpha))
