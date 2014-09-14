#!/usr/bin/env python

"""
Collection of utilities for processing DAKOTA output. Will attempt to convert output
polynomial/coefficient matrix into a string representation suitable for direct 
evaluation

NOTES: 
    1) Agnostic about response function, except written to assume only
       one is present.
"""
import os, pickle, re, sys

import numpy as np
from operator import itemgetter
import pandas as pd
from sympy import Symbol, Rational, Poly, S

from dakota_poly import PolyChaosExp

def poly_around(A, decimals=0):
    """
    Evenly round a polynomial to the given number of decimals.
    """
    B = A.copy()
    for key in A.keys():
        pass

## ORTHOGONAL POLYNOMIALS:
z = Symbol('z')
half = Rational(1, 2)
eighth = Rational(1, 8)
sixteenth = Rational(1, 16)
one_twenty_eigth = Rational(1, 128)


def p(n, z=z): 
    """Legendre polynomial of order `n`"""
    if n == 0:
        poly = Rational(1)
    elif n == 1:
        poly = z
    elif n == 2:
        poly = 0.5*(3*(z**2) - 1)
    elif n == 3:
        poly = 0.5*(5*(z**3) - 3*z)
    elif n == 4:
        poly = 0.125*(35*(z**4) - 30*(z**2) + 3)
    elif n == 5:
        poly = 0.125*(63*(z**5) - 70*(z**3) + 15*z)
    #elif n == 6:
    #    poly = 0.0625*(231*(z**6) - 315*(z**4) + 105*(z**2) - 5)
    #elif n == 7:
    #    poly = 0.0625*(429*(z**7) - 693*(z**5) + 315*(z**3) - 35*z)
    else:
        raise ValueError("n(=%d) too large" % n)
    return poly

def he(n, z=z):
    """Hermite polynomial of order `n`"""
    if n == 0:
        poly = Rational(1)
    elif n == 1:
        poly = z
    elif n == 2:
        poly = z**2 - Rational(1)
    elif n == 3:
        poly = z**3 - 3*z
    elif n == 4:
        poly = z**4 - 6*(z**2) + Rational(3)
    elif n == 5:
        poly = z**5 - 10*(z**3) + 15*z
    elif n == 6:
        poly = z**6 - 15*(z**4) + 45*(z**2) - Rational(15)
    elif n == 7:
        poly = z**7 - 21*(z**5) + 105*(z**3) - 105*z
    elif n == 8:
        poly = z**8 - 28*(z**6) + 210*(z**4) - 420*(z**2) + 105
    else:
        raise ValueError("n(=%d) too large" % n)
    return poly

if __name__ == "__main__":

    directory = sys.argv[1]

    ## READ CONFIGURATION FILE
    setup_dict = pickle.load(open(os.path.join(directory, 'config.p'), 'rb'))
    variables = setup_dict['variables']
    variables = sorted(variables, key=itemgetter(2))
    u = [v[0] for v in variables]
    u_symbols = map(Symbol, u)

    ## READ DAKOTA OUTPUT FILE

    with open(os.path.join(directory, "model_pce.out"), 'r') as f:
        lines = f.readlines()

    ## Seek in output for final polynomial report
    i_beg, i_end = None, None
    for i, line in enumerate(lines):
        if line.startswith("Coefficients of Polynomial Chaos Expansion"):
            i_beg = i+1
        if i_beg and line.startswith("-------------------------"):
            i_end = i
            break

    #print "Coefficients are between %d and %d" % (i_beg, i_end)
    coeff_lines = lines[i_beg:i_end]
    del(coeff_lines[1]) # remove the break row between header and data

    ## Save to scratch file
    coeff_fn = os.path.join(directory, "coeffs_all")
    with open(coeff_fn, "w") as f:
        for line in coeff_lines: 
            #print line
            f.write(line)

    ## Convert to standard CSV
    coeffs = pd.read_table(coeff_fn, delim_whitespace=True)
    coeffs = coeffs[coeffs['coefficient'] != 0.]
    coeffs.to_csv(coeff_fn)
    #print coeffs.head(10)
    #print coeffs

    ## Build internal PCE representation / object
    pce = PolyChaosExp(coeffs, len(u_symbols))
    pickle.dump(pce, open(os.path.join(directory, "pce.p"), "wb"))

    #print pce

##############################
