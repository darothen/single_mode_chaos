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

def extract_pce(responses, lines):
    print "Extracting chaos expansion coefficients..."
    pces = {}
    i_parsed = 0
    for response in responses:
        print "   ", response
        ## Seek in output for final polynomial report
        i_beg, i_end = None, None
        for i, line in enumerate(lines):
            if i < i_parsed: continue
            line = line.strip()
            if line.startswith("Coefficients of Polynomial Chaos Expansion for %s" \
                               % response):
                i_beg = i+1
                print "    FOUND START", i_beg
            elif i_beg:
                if line.startswith("-------------------------") or \
                   line.startswith("Coefficients"):
                    i_end = i
                    print "    FOUND END", i_end
                    break
        i_parsed = i

        print "    Coefficients are between %d and %d" % (i_beg, i_end)
        coeff_lines = lines[i_beg:i_end]

        del(coeff_lines[1]) # remove the break row between header and data

        ## Save to scratch file
        coeff_fn = os.path.join(directory, "%s_coeffs_all" % response)
        with open(coeff_fn, "w") as f:
            for line in coeff_lines: 
                #print line
                f.write(line)

        ## Convert to standard CSV
        coeffs = pd.read_table(coeff_fn, delim_whitespace=True)
        coeffs = coeffs[coeffs['coefficient'] != 0.]
        coeffs.to_csv(coeff_fn)
        #print coeffs

        ## Build internal PCE representation / object
        pces[response] = PolyChaosExp(coeffs, len(u_symbols))

    ## Save to disk
    fn = os.path.join(directory, "pce.p")
    print "writing to %s" % fn
    pickle.dump(pces, open(fn, "wb"))

def extract_sobol(responses, lines):
    print "Extracting Sobol indices..."

    for i, line in enumerate(lines):
        if line.startswith("Global sensitivity indices for each response function:"):
            break
    print "Sobol indices begin at ", i
    lines = lines[i+1:]

    sobols = {}
    for response in responses:
        print "   ", response
        ## 2) Seek forward in output to the matrix of Sobol' indices
        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith("%s Sobol" % response):
                i_beg = i+2
            else:
                if (("Sobol" in line) or 
                   (line.startswith("-------------------"))):
                    i_end = i
                    break
        lines_subset = lines[i_beg:i_end]
        print response, "->", i_beg, i_end

        ## 3) Main and Total interactions
        proc_list = []
        for j, line in enumerate(lines_subset):
            if "Interaction" in line: 
                j += 1
                break
            main, total, term = line.strip().split()
            main, total = float(main), float(total)
            proc_list.append((main, total, term))
        mains, totals, terms = zip(*proc_list)

        ## 4) Multi-term interactions
        proc_list = []
        for k in xrange(j, len(lines_subset)):
            line = lines_subset[k]
            if not line.strip(): break # Capture the blank line ending the block

            bits = line.strip().split()
            interaction = float(bits[0])
            interact_terms = " ".join(sorted(bits[1:]))
            proc_list.append((interaction, interact_terms))
        all_interact, all_interact_terms = zip(*proc_list)

        ## 5) Collect all the main terms together
        mains = mains + all_interact
        all_terms = terms + all_interact_terms
        n_terms = [len(term.split()) for term in all_terms]

        ## 6) Create DataFrame for saving
        main_series = pd.Series(mains, index=all_terms)
        total_series = pd.Series(totals, index=terms)
        nterms_series = pd.Series(n_terms, index=all_terms)

        sobol_df = pd.DataFrame({'Main': main_series, 
                                 'Total': total_series,
                                 'n_terms': nterms_series})
        sobols[response] = sobol_df

        ## Chop off processed lines
        lines = lines[i_end:]

    fn = os.path.join(directory, "sobol.p")
    print "writing to %s" % fn
    pickle.dump(sobols, open(fn, "wb"))

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
    responses = setup_dict['responses']

    ## READ DAKOTA OUTPUT FILE
    with open(os.path.join(directory, "model_pce.out"), 'r') as f:
        lines = f.readlines()

    extract_pce(responses, lines)
    extract_sobol(responses, lines)

##############################
