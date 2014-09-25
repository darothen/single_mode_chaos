import pandas as pd
import pickle

exp_name ="pcm_lasso_parcel"

def extract_sobol(filepath):
    ## 1) Open script for reading
    with open(filepath, "r") as f:
        lines = f.readlines()

    ## 2) Seek forward in output to the matrix of Sobol' indices
    for i, line in enumerate(lines):
        if "Sobol" in line:
            break
    i = i+2 # skip over the blank line and first header
    
    ## 3) Main and Total interactions
    proc_list = []
    for j in xrange(i, len(lines)):
        line = lines[j]
        if "Interaction" in line: 
            j += 1
            break
        main, total, term = line.strip().split()
        main, total = float(main), float(total)
        proc_list.append((main, total, term))
    mains, totals, terms = zip(*proc_list)

    ## 4) Multi-term interactions
    proc_list = []
    for k in xrange(j, len(lines)):
        line = lines[k]
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

    return sobol_df

###########################################################

## Unload the configured experiment
exp_dict = pickle.load(open("%s_exp.dict" % exp_name, 'r'))
results_dict = pickle.load(open("%s_results.dict" % exp_name, 'r'))

print "Reading %s" % exp_name

all_sobols = {}
for key, save_dir in results_dict.iteritems():
    print "   " + key
    out_dir = "save/" + save_dir

    sobol_df = extract_sobol(out_dir + "/model_pce.out")
    sobol_df.to_csv(out_dir+"sobol.df", 
                    index_label='index')
    all_sobols[key] = sobol_df

## Save
with open("%s_sobol.dict" % exp_name, "w") as f:
    pickle.dump(all_sobols, f)

print "... done"