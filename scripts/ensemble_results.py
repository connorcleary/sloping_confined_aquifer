from multiprocessing import Pool
import results_mf6 as mf6r
from pars import ModelParameters, load_parameters
import numpy as np
import matplotlib.pyplot as plt
import os
from mf6 import get_results

def plot_set(inputs):
    
    set, set_name=inputs
    f, axs = plt.subplots(6, len(set)+1, sharex=True)
    mid =  int(np.ceil(0.5*(1+len(set))))

    f.suptitle(set_name)
    
    d = {}

    i = 0
    for r in set:
        if i >= mid-1:
           i=i+1
        d[f"axs3_{i}twin"] = axs[3][i].twinx()
        d[f"axs4_{i}twin"] = axs[4][i].twinx()
        d[f"axs5_{i}twin"]= axs[5][i].twinx()
    
    conc, _, qx, qz, times = get_results("a0")
    conc0p, conc0m = mf6r.plot_results("a0", conc, qx, qz, times, plot=False)
    fresh0, toe0, mix0, sgd0, toe0twin, mix0twin, sgd0twin = mf6r.plot_metrics("a0", times, plot=False)
    base = [conc0p, conc0m, fresh0, toe0, mix0, sgd0, toe0twin, mix0twin, sgd0twin]

    i = 0
    for r in set:
        if i >= mid-1:
           i=i+1
        conc, _, qx, qz, times = get_results(r)
        # mf6r.metrics(r, conc, qz)
        axs[0][i], axs[1][i] = mf6r.plot_results(r, conc, qx, qz, times, plot=False)
        axs[2][i], axs[3][i], axs[4][i], axs[5][i], d[f"axs3_{i}twin"], d[f"axs5_{i}twin"],d[f"axs5_{i}twin"]  = mf6r.plot_metrics(r, times, plot=False)
        i=i+1
     
        axs[0][mid], axs[1][mid], axs[2][mid], axs[3][mid], axs[4][mid], axs[5][mid],d[f"axs3_{mid}twin"], d[f"axs5_{mid}twin"],d[f"axs5_{mid}twin"] = base

    if not os.path.exists(f"/home/ccleary/shared/results/initial_cases"):
        os.makedirs(f"/home/ccleary/shared/results/initial_cases")
    f.savefig(f"/home/ccleary/shared/results/initial_cases/{set_name}", dpi=300) 
        


def plot_ensemble(sets, set_names, name):

    # get everything for base case
    conc, _, qx, qz, times = get_results("a0")
    mf6r.metrics("a0", conc, qz)
    conc0p, conc0m = mf6r.plot_results("a0", conc, qx, qz, times, plot=False)
    fresh0, toe0, mix0, sgd0, toe0twin, mix0twin, sgd0twin = mf6r.plot_metrics("a0", times, plot=False)
    base = [conc0p, conc0m, fresh0, toe0, mix0, sgd0, toe0twin, mix0twin, sgd0twin]
    # plot_set(sets[0], set_names[0], base)
    inputs = [[s, n] for s, n in zip(sets, set_names)]
    p = Pool(processes=16)
    p.map(plot_set, inputs)


def main():
    sets = [[f"a{i}", f"a{i+1}"] for i in range(1, 13, 2)]
    print(sets)
    set_names = ["K_aquifer", "K_aquitard", "anis_aquifer", "anis_aquitard", "h_modern", "aquitard_D"]
    plot_ensemble(sets, set_names, "inital_cases")

if __name__=="__main__":
        main()
