# Crude Python results

import numpy as np
import pandas as pd



def run_crude(date, positivity, effect, pw):
    
    # Construct file name and read it in. 
    filename = "{0}-{1}-{2}".format(positivity, effect, pw)
    infile = "data/population/{0}/{1}.csv".format(date, filename)
    full = pd.read_csv(infile)

    # Intercept for PS mod4-08
    full['intercept'] = 1 # Intercept value

    # Number of simulations to loop through
    sim_start = 1
    sim_end = 20
    sim_vec = np.arange(sim_start, sim_end+1)

    # Parameters to estimate
    params = ['mu_crude_0', 'mu_crude_1', 'delta_crude']

    # Create DataFrame to insert simulation values into
    final = pd.DataFrame(index = range(0,(sim_end + 1 - sim_start)*len(params)), columns = ["Sim", "Param", "Coef", "Variance", "LCL", "UCL"])
    final['Param'] = params*(sim_end + 1 - sim_start)
    final['Sim'] = np.repeat(sim_vec, len(params))

    for i in range(sim_start, sim_end+1):
        d = full.loc[full["sim"]==i]
        tab = pd.crosstab(d['a'], d['y'], margins=True)
        
        # Risks and variances
        mu1_crude = tab.loc[1,1]/tab.loc[1,'All']
        mu1_var = mu1_crude*(1-mu1_crude)/tab.loc[1, 'All']
        mu0_crude = tab.loc[0,1]/tab.loc[0,'All']
        mu0_var = mu0_crude*(1-mu0_crude)/tab.loc[0, 'All']

        # Risk ratio and variance
        lnrr_crude = np.log(mu1_crude/mu0_crude)
        lnrr_var_crude = 1/tab.loc[1,1]-1/tab.loc[1, 'All']+1/tab.loc[0,1]-1/tab.loc[0,'All']

        ests = np.array([mu1_crude, mu0_crude, lnrr_crude])
        vars = np.array([mu1_var, mu0_var, lnrr_var_crude])

        # CIs
        lcl = ests - 1.96*np.sqrt(vars)
        ucl = ests + 1.96*np.sqrt(vars)

        final.loc[final['Sim'] == i, 'Coef'] = ests
        final.loc[final['Sim'] == i, 'Variance'] = vars
        final.loc[final['Sim'] == i, 'LCL'] = lcl
        final.loc[final['Sim'] == i, 'UCL'] = ucl

    outfile = "data/simulations/{0}/crude-{1}-{2}-{3}.csv".format(date, filename, str(np.min(sim_vec)).zfill(5), str(np.max(sim_vec)).zfill(5))
    final.to_csv(outfile, index = False)

newdate = "2024-04-29"
run_crude(date = newdate, positivity = "full", effect = "homogeneous", pw = "75")
run_crude(date = newdate, positivity = "full", effect = "homogeneous", pw = "50")
run_crude(date = newdate, positivity = "full", effect = "homogeneous", pw = "25")

run_crude(date = newdate, positivity = "full", effect = "none", pw = "75")
run_crude(date = newdate, positivity = "full", effect = "none", pw = "50")
run_crude(date = newdate, positivity = "full", effect = "none", pw = "25")

run_crude(date = newdate, positivity = "partial", effect = "homogeneous", pw = "75")
run_crude(date = newdate, positivity = "partial", effect = "homogeneous", pw = "50")
run_crude(date = newdate, positivity = "partial", effect = "homogeneous", pw = "25")

run_crude(date = newdate, positivity = "partial", effect = "none", pw = "75")
run_crude(date = newdate, positivity = "partial", effect = "none", pw = "50")
run_crude(date = newdate, positivity = "partial", effect = "none", pw = "25")