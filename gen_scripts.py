#!/usr/bin/env python
import pickle
import sys

from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader('templates'))

if __name__ == "__main__": 

    setup_dict = pickle.load(open('config.p', 'rb'))

    ## Model run script
    print "Writing model_run.py with",
    template_fn = "model_run_script.template" if setup_dict['parallel'] else \
                  "model_run_linked.template"  
    print template_fn, "..."
    template = env.get_template(template_fn)
    template_items = {
        "imports": setup_dict['imports'], 
        "function_name": setup_dict['function_name'],
        "function": setup_dict['function'],
        "variables": setup_dict['variables'],
        "variable_symbols": [v[0] for v in setup_dict['variables']],
    }

    output = template.render(**template_items)

    with open("model_run.py", "wb") as f:
        f.write(output)

    ## DAKOTA control script
    print "Writing model_pce.in with",
    uniform_vars, loguniform_vars, normal_vars = [], [], []
    for v in setup_dict['variables']:
        dist_type = v[3][0]
        if dist_type == 'uniform':
            uniform_vars.append(v)
        elif dist_type == 'loguniform':
            loguniform_vars.append(v)
        elif dist_type == 'normal':
            normal_vars.append(v)
        else:
            raise ValueError("Distribution %s not known" % dist_type)

    template_fn = "model_pce_parallel.template" if setup_dict['parallel'] else \
                  "model_pce.template"  
    print template_fn, "..."
    template = env.get_template(template_fn)

    template_items = {
        "graphics": setup_dict['graphics'],
        "pce_directive": setup_dict['pce_directive'],
        "uniform_vars": uniform_vars,
        "loguniform_vars": loguniform_vars,
        "normal_vars": normal_vars,
        "responses": setup_dict["responses"],
    }

    output = template.render(**template_items)

    with open("model_pce.in", "wb") as f:
        f.write(output)
        