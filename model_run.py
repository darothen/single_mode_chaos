#!/usr/bin/env python
import sys, shutil

import numpy as np 


def arg(lnN, mu, sigma, kappa, lnV, T, P, fn_toggle=None):
    import numpy as np
    import parcel_model as pm
        
    N = np.exp(lnN)
    V = np.exp(lnV)
    #mu = np.exp(lnmu)
    
    if not fn_toggle:
    # Parameterization
    #    Smax, _, _ = pm.arg2000(V, T, P, 
    #                            mus=[mu, ], sigmas=[sigma, ], Ns=[N, ], 
    #                            kappas=[kappa, ])
    #    return np.log(Smax)
    ## Parcel Model
        aerosol_modes = [
            pm.AerosolSpecies('aer', 
                               pm.Lognorm(mu=mu, sigma=sigma, N=N),
                               kappa=kappa, bins=200),
        ]
        try:
            model = pm.ParcelModel(aerosol_modes, V, T, -0.0, P, accom=0.1)
            par_out, aer_out = model.run(t_end=1800., dt=0.01, max_steps=2000,
                                         solver='cvode', output='dataframes',
                                         terminate=True)
                                         #solver_args={'time_limit': 10.})
            Smax = par_out['S'].max()
        except pm.ParcelModelError:
            Smax = -1.0
            
        if Smax < 0: 
            #return -9999.0
            Smax, _, _ = pm.arg2000(V, T, P, 
                                    mus=np.array([mu, ]), 
                                    sigmas=np.array([sigma, ]), 
                                    Ns=np.array([N, ]), 
                                    kappas=np.array([kappa, ]), accom=0.1)
            print ">>>> model integration FAILED", V, T, P, mu, sigma, kappa, N, "|", Smax 

            return np.log(Smax)
        else:
            return np.log(Smax)            
    else:
        mu *= 1e-6
        N_act, act_frac = pm.activate_lognormal_mode(fn_toggle, 
                                                     mu, sigma, N, kappa, T=T)
        return N_act
    


if __name__ == "__main__":

    fin, fout = sys.argv[1:]
    proc_id = int(fin.split(".")[-1])

    ## Read the input file
    with open(fin, 'rb') as f:
        lines = f.readlines()
        lnN = float(lines[1].strip().split()[0])
        mu = float(lines[2].strip().split()[0])
        sigma = float(lines[3].strip().split()[0])
        kappa = float(lines[4].strip().split()[0])
        lnV = float(lines[5].strip().split()[0])
        T = float(lines[6].strip().split()[0])
        P = float(lines[7].strip().split()[0])

    shutil.copy2(fin, "bak/"+fin)

    #############################

    results = arg(lnN, mu, sigma, kappa, lnV, T, P)
        
    #############################

    if not isinstance(results, list):
        results = [results, ]

    ## Write the output file
    with open(fout, "wb") as f:
        for r in results:
            f.write("%e f\n" % r)