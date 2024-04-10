import numpy as np
import pandas as pd

import delicatessen
from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_regression, ee_rogan_gladen
from delicatessen.utilities import inverse_logit

# Define m-estimation function 

def psi(theta):
    # Dividing parameters into corresponding parts and labels from slides
    alpha = theta[0:2]              # Logistic model coefficients
    mu_ht_0, mu_ht_1 = theta[2], theta[3]   # HT Causal risks
    mu_hajek_0, mu_hajek_1 = theta[4], theta[5]   # HT Causal risks
    mu_smr_0, mu_smr_1 = theta[6], theta[7]
    delta_ht, delta_hajek, delta_smr = theta[8], theta[9], theta[10]    # Causal contrasts

    # Logistic regression model for propensity score
    ee_logit = ee_regression(theta=alpha,       # Regression model
                            y=a,               # ... for exposure
                            X=W,               # ... given confounders
                            model='logistic')  # ... logistic model

    # Transforming logistic model coefficients into causal parameters
    pscore = inverse_logit(np.dot(W, alpha))    # Propensity score
    ipt = d['a']/pscore + (1-d['a'])/(1-pscore)  # Corresponding IPT weights
    smr = d['a']*1 + (1-d['a'])*pscore/(1-pscore)

    # HT IPT Estimating function for causal risk under a=1
    ee_ht_r1 = d['a']*d['y']*ipt - mu_ht_1 # Weighted conditional mean
    # HT IPT Estimating function for causal risk under a=0
    ee_ht_r0 = (1-d['a'])*d['y']*ipt - mu_ht_0 # Weighted conditional mean

    # Hajek IPT Estimating function for causal risk under a=1
    ee_hajek_r1 = d['a']*ipt*(d['y']-mu_hajek_1) # Weighted conditional mean
    # Hajek IPT Estimating function for causal risk under a=0
    ee_hajek_r0 = (1-d['a'])*ipt*(d['y']-mu_hajek_0)  # Weighted conditional mean

    # SMR IPT Estimating function for causal risk under a=1
    ee_smr_r1 = d['a']*smr*(d['y']-mu_smr_1) # Weighted conditional mean
    # Hajek IPT Estimating function for causal risk under a=0
    ee_smr_r0 = (1-d['a'])*smr*(d['y']-mu_smr_0)  # Weighted conditional mean


    # Estimating function for HT causal risk difference
    ee_ht_rr = np.ones(d.shape[0])*((np.log(mu_ht_1) - np.log(mu_ht_0)) - delta_ht)
    # Estimating function for Hajek causal risk difference
    ee_hajek_rr = np.ones(d.shape[0])*((np.log(mu_hajek_1) - np.log(mu_hajek_0)) - delta_hajek)
    # Estimating function for SMR causal risk difference
    ee_smr_rr = np.ones(d.shape[0])*((np.log(mu_smr_1) - np.log(mu_smr_0)) - delta_smr)

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_logit,   # EF of logistic model
                    ee_ht_r0,      # HT EF of causal risk a=0
                    ee_ht_r1,      # HT EF of causal risk a=1
                    ee_hajek_r0,   # Hajek EF of causal risk a=0
                    ee_hajek_r1,   # Hajek EF of causal risk a=1
                    ee_smr_r0,     # SMR EF of causal risk a=0
                    ee_smr_r1,     # SMR EF of causal risk a=1
                    ee_ht_rr,      # HT EF of causal contrast
                    ee_hajek_rr,   # Hajek EF of causal contrast
                    ee_smr_rr ])   # SMR EF of causal contrast


# Read in population 
def run_mestimation(date, positivity, effect, pw, iter):
    filename = "{0}-{1}-{2}".format(positivity, effect, pw)
    infile = "data/population/{0}/{1}.csv".format(date, filename)
    full = pd.read_csv(infile)
    
    full['intercept'] = 1 # Intercept value

    # Number of simulations to loop through
    sim_start = 1
    sim_end = 200
    sim_vec = np.arange(sim_start, sim_end+1)

    # Parameters to estimate
    params = ['alpha_0', 'alpha_1', 'mu_ht_0', 'mu_ht_1', 'mu_hajek_0', 'mu_hajek_1', 'mu_smr_0', 'mu_smr_1', 'delta_ht', 'delta_hajek', 'delta_smr']

    # Create DataFrame to insert simulation values into
    final = pd.DataFrame(index = range(0,(sim_end + 1 - sim_start)*len(params)), columns = ["Sim", "Param", "Coef", "Variance", "LCL", "UCL"])
    final['Param'] = params*(sim_end + 1 - sim_start)
    final['Sim'] = np.repeat(sim_vec, len(params))

    for i in range(sim_start, sim_end+1):
        d = full.loc[full["sim"]==i]
        a = np.asarray(d['a'])
        W = np.asarray(d[['intercept','w']])
        
        try:        
            estr = MEstimator(psi, init=[-1, -1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, np.log(0.8), np.log(0.8), np.log(0.8)])
            estr.estimate(maxiter = iter)

            final.loc[final['Sim'] == i, 'Coef'] = estr.theta
            final.loc[final['Sim'] == i, 'Variance'] = np.diag(estr.variance)
            ci = estr.confidence_intervals()
            final.loc[final['Sim'] == i, 'LCL'] = ci[:, 0]
            final.loc[final['Sim'] == i, 'UCL'] = ci[:, 1]
        except: 
            final.loc[final['Sim'] == i, 'Coef'] = [np.nan, ]*len(params)
            final.loc[final['Sim'] == i, 'Variance'] = [np.nan, ]*len(params)
            final.loc[final['Sim'] == i, 'LCL'] = [np.nan, ]*len(params)
            final.loc[final['Sim'] == i, 'UCL'] = [np.nan, ]*len(params)

    outfile = "data/simulations/{0}/{1}-{2}-{3}.csv".format(date, filename, str(np.min(sim_vec)).zfill(5), str(np.max(sim_vec)).zfill(5))
    final.to_csv(outfile, index = False)


delta_none = 0
delta_diff = np.log(0.8)

## Look into partial pos scenarios starting values/ print thetas
run_mestimation(date = "2024-04-08", positivity = "full", effect = "homogeneous", pw = "75", iter = 5000)
run_mestimation(date = "2024-04-08", positivity = "full", effect = "homogeneous", pw = "50", iter = 5000)
run_mestimation(date = "2024-04-08", positivity = "full", effect = "homogeneous", pw = "25", iter = 5000)

run_mestimation(date = "2024-04-08", positivity = "full", effect = "none", pw = "75", iter = 5000)
run_mestimation(date = "2024-04-08", positivity = "full", effect = "none", pw = "50", iter = 5000)
run_mestimation(date = "2024-04-08", positivity = "full", effect = "none", pw = "25", iter = 5000)

run_mestimation(date = "2024-04-08", positivity = "partial", effect = "homogeneous", pw = "75", iter = 75000)
run_mestimation(date = "2024-04-08", positivity = "partial", effect = "homogeneous", pw = "50", iter = 75000)
run_mestimation(date = "2024-04-08", positivity = "partial", effect = "homogeneous", pw = "25", iter = 75000)

run_mestimation(date = "2024-04-08", positivity = "partial", effect = "none", pw = "75", iter = 75000)
run_mestimation(date = "2024-04-08", positivity = "partial", effect = "none", pw = "50", iter = 75000)
run_mestimation(date = "2024-04-08", positivity = "partial", effect = "none", pw = "25", iter = 75000) 