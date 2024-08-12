from scipy.stats import qmc
import pdb
import numpy as np
import os
import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool
import post_processing as proc
import results_mf6
from results_mf6 import save_results_for_saltelli, plot_color_mesh_saltelli, get_pars, get_outputs, paleo_volumes
from SALib.sample import sobol as sbsample
from SALib.analyze import sobol
import pairs_plot
import paper_plots
import pickle
import supporting_plots as sp
import itertools

def model_run(inputs):
    # function to run model with fixed geometry and variable hydraulic gradient
    ensemble_name, name, real = inputs
    pars = ModelParameters(name, H=20, D=20, K_aquifer=real[0], 
                           K_aquitard=real[1]*10, alpha_L=real[2], anis_aquifer=1, anis_aquitard=10,
                           z0=55, L=30000, x_b=4000,  
                           h_modern=-1, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=real[3], dz=real[4], dt=2500)    
    sim = mf6.build_steady_model(pars)
    succ = True
    try:
        mf6.run_model(sim)
    except:
        succ= False
    
    if succ: 
        try:
            save_results_for_saltelli(name, ensemble_name)
        except: 
            pass
        
    mf6.remove_model_files(name)  

def plot_raleighs(name, n):
    for i in range(n):
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_results.npy")
        conc = results[0][0]
        qx = results[0][2]
        qz = results[0][3]
        pars = load_parameters(name)
        pass

def main():
    # Main function for running sensitivity analysis.
    # The difference between this and other main function is that aquifer geometry is constant
    # An added sensitivity is to the predevelopment hydraulic gradient
    name = "dis"
    Ka = [1 ,100]
    KbV = [ 1e-4, 1e-3]
    alpha = [0.1, 10]
    dx = [25, 5, 1]
    dz = [1, 0.5]
    values = np.array(list(itertools.product(*[Ka, KbV, alpha, dx, dz])))
    inputs = []
    for i in range(48):
        inputs.append(["dis", f"dis{i}", values[i]])

    model_run(inputs[2])  
    #p=Pool(processes=8)
    # p.map(model_run, inputs)
#     plot_color_mesh_saltelli("dis", 48)

if __name__ == "__main__":
    main()
    # plot_raleighs("dis", 1)