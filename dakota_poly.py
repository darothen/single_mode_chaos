"""
Simple representation of polynomial output by DAKOTA PCE routines.

"""

import os, pickle, re, sys

import numpy as np
import pandas as pd
from scipy.special import orthogonal as op

class PolyChaosExp():

    def __init__(self, coeffs, nvars, debug=False):
        """

        Parameters
        ----------
        coeffs : pandas.DataFrame 
            the coefficient table output at the end of a PCE run by DAKOTA,
            saved as a pandas DataFrame

        nvars : int
            the number of variables used in the PCE parameter space

        Attributes

        A - coefficients
        orders - matrix of (ncoeffs x nvars) with the order of ecah poly
        polys - type of poly

        """
        self.nvars = nvars
        self.debug = debug

        ## Sanitize coefficients
        coeffs = coeffs[coeffs['coefficient'] != 0.]

        self.A = np.asarray(coeffs['coefficient'])
        self.nterms = len(self.A)

        ## Construct polynomial order matrix:
        ##    rows - each term in the polynomial
        ##    cols - the degree of each polynomial mixin
        self.orders = []
        for i, (ind, row) in enumerate(coeffs.iterrows()):
            if self.debug: print i, ind, row
            row_orders = []
            first = i == 0

            if first: ## We want to save the orthogonal poly! 
                self.polys = []
            for i in xrange(1, nvars+1):
                poly, order = self._separate(row["u%d" % i])
                if first: self.polys.append(poly)

                row_orders.append(order)
            self.orders.append(row_orders)
        self.orders = np.array(self.orders)

        ## Save the maximum polynomial to add in orthogonal eval
        self.max_orders = np.max(self.orders, axis=0)

    def evaluate(self, z):
        """
        Evaluate a chaos expansion at the given parameter tuple z
        """

        ## Sanity checks:
        #  1) transform z into an array
        z = np.asarray(z)
        #  2) make sure it has the right amount of terms
        assert len(z) == self.nvars

        # Step 1): Evaluate the necessary polynomial orders
        poly_evals = []
        for zi, poly, order in zip(z, self.polys, self.max_orders):
            fn = self._poly_fcn_map(poly)
            poly_evals.append([fn(n, zi) for n in xrange(order+1)])

        # Step 2): Loop through arrays and evaluate polynomial
        evaluation = 0.
        for i in xrange(self.nterms):
            if self.debug: print i, self.A[i], self.orders[i]
            row = self.orders[i]
            
            row_eval = 1.
            for n in xrange(self.nvars):
                order = row[n]
                p_zi = poly_evals[n][order]
                if self.debug: print "   ", n, p_zi
                row_eval *= p_zi
            if self.debug: print "RUNNING EVAL -> ", evaluation,
            evaluation += self.A[i]*row_eval
            if self.debug: print evaluation

        return evaluation

    ###############

    def _separate(self, s):
        PATTERN = r"([a-z]+)([0-9]+)"   
        matches = re.match(PATTERN, s, re.IGNORECASE)
        if matches:
            items = matches.groups()
            poly = items[0]
            order = int(items[1])
            return poly, order
        else:
            raise ValueError("Could not understand string, '%s'" % s)


    def _poly_fcn_map(self, poly):
        mapping = {
            'P': op.eval_legendre,
            'He': op.eval_hermitenorm,
            #'He': eval_hermiteprob,
        }
        return mapping[poly]

    ###############

    def __call__(self, z):
        return self.evaluate(z)


