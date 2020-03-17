#!/usr/bin/env python

# split into four populations, asymmetric migration in each epoch
# n(para): 32

import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import random
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
import pandas as pd

#infile="c1234_dadi.data"
#pop_ids=["pop0","pop1","pop2","pop3"]
#projections=[int(69),int(15),int(67),int(40)]
#mu=float(0.02)
#gtime=float(0.005)

infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]]
projections=[int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9])]
#params=[float(sys.argv[6]),float(sys.argv[7]),float(sys.argv[8]),float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11])]
params=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

# mutation rate per sequenced portion of genome per generation: for A.millepora, 0.02
mu=float(sys.argv[10])
# generation time, in thousand years: 0.005  (5 years)
gtime=float(sys.argv[11]) 


#reorder data to match order of population divergences
df = pd.read_table(infile)
cols = df.columns.tolist()
pop_order = [3,1,0,2]
col_order = [0,1,2] + [3+x for x in pop_order] + [7] + [8+x for x in pop_order] + [12,13]
cols = [cols[i] for i in col_order]
df = df[cols]
df2 = df.rename(columns={"pop0.1": "pop0", "pop1.1": "pop1", "pop2.1": "pop2", "pop3.1": "pop3"})
df2.to_csv(infile+"_temp", sep = "\t", index = False)
pop_ids = [pop_ids[i] for i in pop_order]
projections = [projections[i] for i in pop_order]

dd = Misc.make_data_dict(infile+"_temp")
# set Polarized=False below for folded AFS analysis
data = Spectrum.from_data_dict(dd, pop_ids, projections, polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# split into unequal pop sizes with asymmetrical migration

def sc1234(params,ns):
    nu123_1,nu4_1,nu13_2,nu2_2,nu4_2,nu1_3,nu2_3,nu3_3,nu4_3,T1,T2,T3,m4_123,m123_4,m4_13,m2_13,m13_4,m2_4_e2,m13_2,m4_2_e2,m4_1,m2_1,m3_1,m1_4,m2_4_e3,m3_4,m1_2,m4_2_e3,m3_2,m1_3,m4_3,m2_3 = params
	
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0]+ns[1]+ns[2]+ns[3])
    fs = moments.Spectrum(sts)
    
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1]+ns[2]+ns[3])
    #row: sink population, column: source population
    mig1 = np.array([[0, m123_4], [m4_123, 0]])
    fs.integrate([nu4_1, nu123_1], T1, m = mig1)
    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2]+ns[3])
    mig2 = np.array([[0, m2_4_e2, m13_4], [m4_2_e2, 0, m13_2], [m4_13, m2_13, 0]])
    fs.integrate([nu13_2, nu2_2, nu4_2], T2, m = mig2)
    
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])
    mig3 = np.array([[0, m2_4_e3, m1_4, m3_4], [m4_2_e3, 0, m1_2, m3_2], [m4_1, m2_1, 0, m3_1], [m4_3, m2_3, m1_3, 0]])
    fs.integrate([nu1_3, nu2_3, nu3_3, nu4_3], T3, m = mig3)
    
    return fs
 
func=sc1234
upper_bound = [100,100,100,100,100,100,100,100,100,100,100,100,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

poptg = moments.Inference.optimize_log(params, data, func, lower_bound=lower_bound, upper_bound=upper_bound, verbose=False, maxiter=30)

# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="sc1234_"+ind+".png", pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta
all_boot=moments.Misc.bootstrap(dd,pop_ids,projections)
uncert=moments.Godambe.GIM_uncert(func,all_boot,poptg,data)

# printing parameters and their SDs
print "RESULT","sc1234",ind,len(params),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],poptg,theta,uncert
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3, pop_ids=pop_ids)
plt.savefig("sc1234_"+ind+"_"+sys.argv[1]+'.pdf')

