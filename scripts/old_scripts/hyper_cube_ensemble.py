from scipy.stats import qmc
import pdb
import numpy as np
import os
import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool
import post_processing as proc
from mf6 import get_results, get_results_for_hc
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.tri as tri
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from pyeee import ee

# produce samples 
def prod_samples(name, n, d, ranges, scales):
    sampler = qmc.LatinHypercube(d=d)
    sample = sampler.random(n=n)
    
    # for each column
    sample_par_space = np.zeros_like(sample)
    for i, (rang, scale, column) in enumerate(zip(ranges, scales, sample.T)):
        if scale == "log":
            func = lambda x : rang[0]*10**(np.log10(rang[1]/rang[0])*x)
            sample_par_space[:, i] = [func(x) for x in column]
        else:
            func = lambda x : rang[0] +(rang[1]-rang[0])*x
            sample_par_space[:, i] = [func(x) for x in column]
    	   
    # save 
    if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/"):
        os.makedirs(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/")
    
    np.savetxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", sample_par_space, 
    		delimiter=",")
    
    return sample_par_space
    
# create and run model
def create_run_model(name, real):

    pars = ModelParameters(name, H=real[0], D=real[1], K_aquifer=real[2], K_aquitard=real[3], alpha_L=real[4], anis_aquifer=real[5], z0=np.tan(real[6])*real[7]+real[0]+real[1], L=real[7], h_modern=real[8], porosity=real[9], outputs="last")    
    sim = mf6.build_steady_model(pars)
    try:
        mf6.run_model(sim)
        return True
    except:
        return False

# find result
def find_metrics(name, ensemble_name):
	
    pars = load_parameters(name)
    conc, qz = get_results_for_hc(name)
    fresh_list = proc.fresh_volume([conc[-2],conc[-1]], pars)
    mix_list = proc.mixing([conc[-2],conc[-1]], pars)
    if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/"):
        os.makedirs(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/")
        
    np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}.csv", fresh_list+mix_list[0], delimiter=",")
    np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}_qz.csv", qz[:, 0, :], delimiter=",")
   
    return 0.0


# parallel func
def run_hc_real(inputs):
    ensemble_name, name, real = inputs
    succ = create_run_model(name, real)
    
    if succ: 
        try:
            _ = find_metrics(name, ensemble_name)
        except: 
            pass 

    mf6.remove_model_files(name)
    

    	
def run_hc_real_single(inputs):
    ensemble_name, name, real = inputs
    succ = create_run_model(name, real)
    
    if succ: 
        _ = find_metrics(name, ensemble_name)
    	
    mf6.remove_model_files(name)
    
def extract_results(ensemble_names):

    for name in ensemble_names:
        results = []
        for i in range(50):# len(os.listdir(f"/home/superuser/sloping_confined_aquifer/results/{name}/"))):
            try:
                output = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}.csv", delimiter=",")
                results.append(output)
            except: 
                results.append([np.nan, np.nan, np.nan, np.nan,])
        np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", np.array(results), delimiter=",")
        
        
def extract_fluxes(ensemble_names):

    
    for name in ensemble_names:
        results = []
        for i in range(50):
            pars = load_parameters(f"{name}{i}", backup=True)
            qz = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_qz.csv", delimiter=",")
            flux=[]
            #pdb.set_trace()
            for cell in pars.top_cells:
                if cell[2] >= pars.x_b/pars.dx:
                    flux.append(qz[cell[0]][cell[2]])
            results.append(flux)
        results_np = np.zeros((len(results), len(max(results, key=len))))
        for i, result in enumerate(results):
            results_np[i][:len(result)]=result
        np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all_fluxes.csv", results_np, delimiter=",")
    # pdb.set_trace()


def pairs_plot(results, scales, pars, units, lims=None, title=None):
    f, axs = plt.subplots(len(units), len(units), figsize=(6,6))
    #axs[0, 1].axis('off')
    
    for row in range(len(units)):
        for col in range(len(units)):
            if row>col:
                # pdb.set_trace()                
                            
                if scales[row] == "log":
                    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                axs[row, col].scatter(results[:, col].T, results[:, row].T, c=results[:, -1].T, marker="o", alpha=0.7)   
                             
            elif row == col:
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")             
                axs[row, col].scatter(results[:, col].T, results[:, -1].T, c="grey", marker="o", alpha=0.7)
            else:
                axs[row, col].axis('off')

            
            if col == 0 and row == len(units)-1:
                axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
            elif col == 0:
                if row == 0:
                    axs[row, col].set_ylabel(r"$\frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
                else:# [{units[row]}]")
                    axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                plt.setp(axs[row, col].get_xticklabels(), visible=False)
            elif row == len(units)-1:
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
                axs[row, col].set_yticklabels([])
            else:
                axs[row, col].set_xticklabels([])
                axs[row, col].set_yticklabels([])
    
            axs[row, col].set_box_aspect(1)  

    #cb = plt.colorbar(scatter, ax=axs[0, 1], cmap="viridis", extend="both", location="right")
    #cb.ax.set_title(f"{title}")
    #cb.ax.yaxis.set_ticks_position("left")
    #plt.tight_layout()       
   
    plt.show()
    
    
def threeDscatter(name, results, scales):
    
    f = plt.figure()
    ax = f.add_subplot(projection='3d')
    ax.scatter(np.log10(results[:, 0]), np.log10(results[:, 1]), results[:, 2], c=results[:, 3], cmap="viridis", vmin=0, vmax=1)
    # ax.scatter(np.log10(results[:, 0]), np.log10(results[:, 1]), results[:, 2], cmap="viridis", vmin=0, vmax=1)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.set_xlabel(r"$K_{aquifer}$ [log(m/day)]")
    ax.set_ylabel(r"$K_{aquitard}$ [log(m/day)]")
    ax.set_zlabel(r"$h_{post}$ [masl]")
    # ax.set_box_aspect([1,1,1])
    plt.show()
    plt.savefig(f"/home/superuser/sloping_confined_aquifer/results/{name}/3D.png")

def fit_polynomial_regression(x1, x2, y, scales, lims):

    poly = PolynomialFeatures(degree=2)
    if scales[0]=="log":
        x1 = np.log10(x1)
        x1_mesh = np.log10(np.logspace(np.log10(lims[0][0])/np.log10(10),np.log10(lims[0][1]/np.log10(10))))
    else: 
        x1_mesh = np.linspace(lims[0][0],lims[0][1])
    if scales[1]=="log":
        x2 = np.log10(x2)
        x2_mesh = np.log10(np.logspace(np.log10(lims[1][0])/np.log10(10),np.log10(lims[1][1]/np.log10(10))))
    else:
        x2_mesh = np.linspace(lims[1][0],lims[1][1])
    	
    X1, X2 = np.meshgrid(x1_mesh, x2_mesh)
    
    x = np.dstack([x1, x2])[0]
    x_ = poly.fit_transform(x)
    X = np.dstack([X1.flatten(),X2.flatten()])[0]
    X_ = poly.fit_transform(X)
    
    clf = LinearRegression()
    clf.fit(x_, y)
    
    Y = clf.predict(X_)
    if scales[0]=="log":
        X1 = np.power(10, X1)
    if scales[1]=="log":
    	X2 = np.power(10, X2)
        
    return X1, X2, Y.reshape(50, -1)

def main():
    name = "paper_test_one"
    # samples = prod_samples(name, 50, 10, 
    #     [[10, 50], # aquifer thickness
    #      [10, 20], # aquitard thickness
    #      [1, 1000], # K aquifer
    #      [0.001, 0.0001], # K aquitard
    #      [0.1, 10], # Dispersivity
    #      [1, 100], # Anisotropy
    #      [1e-4, 2e-3], # aquifer slope
    #      [10000, 20000], # aquifer length
    #      [2,-10],
    #      [0.2,0.8]], # head change
    #     ["log", "linear", "log", "log", "linear", "linear", "linear", "linear", "linear", "linear"])
    samples = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/paper_test_one.csv", delimiter=",")
    #inputs = [[name, f"{name}{i}", sample] for i, sample in enumerate(samples)]

    #p=Pool(processes=4)
    #p.map(run_hc_real, inputs)
    #for i in range(len(inputs)):
    #run_hc_real_single(inputs[0])
    extract_results([name])

    #scales = ["linear", "linear", "log", "log", "linear", "linear", "linear", "linear", "linear", "linear"]
    scales = ["linear", "log", "log", "linear", "linear", "linear", "linear"]
    #pars = ["H", "D", r"$K_a$", "$K_b$", r"$\alpha$", "anis", r"$\beta$", "L", r"$h_{post}$", "n"]
    pars = ["H", r"$K_a$", "$C_b^z$", r"$\alpha$", r"$\beta$", "L", r"$h_{post}$"]
    #units = ["m", "m", "m/day", "m/day", "-", "-", "radians", "m", "m", "-"]
    units = ["m", "m/day", "m/day", "-", "radians", "m", "m"]
    results = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", delimiter=",")   
    inputs = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", delimiter=",")
    exhaustion_metric=(results[:, 0]-results[:, 1])/results[:, 0]
    inputs_cleaned = inputs[~np.isnan(exhaustion_metric)]
    exhaustion_cleaned = exhaustion_metric[~np.isnan(exhaustion_metric)]
    inputs_cleaned_post = np.copy(inputs_cleaned[:, :-1])
    inputs_cleaned_post[:, 3] = (inputs_cleaned_post[:, 3]/inputs_cleaned_post[:, 5])/inputs_cleaned_post[:, 1]
    inputs_cleaned_post = np.delete(inputs_cleaned_post, 5, 1)
    inputs_cleaned_post = np.delete(inputs_cleaned_post, 1, 1)
    pairs_plot(np.concatenate((inputs_cleaned_post, np.atleast_2d(exhaustion_cleaned).T), axis=1), scales, pars, units)
    
    # results = collect_results(name, ensemble_pars=False)
    # results = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", delimiter=",")
    # pairs_plot(name, results, ["log", "log", "linear"], [r"$K_{aquifer}$", r"$K_{aquitard}$", r"$h_{post}$"], ["$m/day$", "$m/day$", "$masl$"], lims=[[1, 100],[0.001, 1],[0, -10]])
    # threeDscatter(name, results, ["log", "log", "linear"])
    
    
if __name__ == "__main__":
    #find_metrics("all0", "all")
    main()   	
    
