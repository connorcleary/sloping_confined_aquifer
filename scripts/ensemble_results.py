from multiprocessing import Pool
import results_mf6 as mf6r
from pars import ModelParameters, load_parameters
import numpy as np
import matplotlib.pyplot as plt
import os
from mf6 import get_results
import pandas as pd

def plot_metrics(name, times, axs2, axs3, axs4, axs5, axs3twin, axs4twin, axs5twin, y_labels="none", legend=False):

    df = pd.read_csv(f"/home/ccleary/shared/results/{name}/metrics.csv")
    
    # plot fresh
    axs2.plot(times[0::5], df.fresh_list)
    axs3twin.plot(times[0::5], df.toe_list)
    axs3twin.plot(times[0::5], df.tip_list)
    axs4.plot(times[0::5], df.volume_list)
    axs4twin.plot(times[0::5], df.centroid_list, linestyle="--")
    axs5.plot(times[0::5], df.fresh_flux_list)
    axs5.plot(times[0::5], df.sal_flux_list)
    axs5.set_xlabel("time from present (days)")
    axs5twin.plot(times[0::5], df.fresh_centroid_list, linestyle="--")
    axs5twin.plot(times[0::5], df.sal_centroid_list, linestyle="--")


    if y_labels == "left":
        axs2.set_ylabel("fresh volume (m^3)")
        axs4.set_ylabel("volume (m^3)")
        axs5.set_ylabel("flux (m^3)")
        axs2.legend(["fresh volume"])
        axs4.legend(["mixing zone"])
        axs5.legend(["fresh sgd", "saline sgd"])
        
    elif y_labels == "right":
        axs3twin.set_ylabel("distance offshore (m)")
        axs3twin.legend(["toe", "tip"])
        axs4twin.set_ylabel("distance offshore (m)")
        axs4twin.legend(["centoid"])
        axs5twin.set_ylabel("distance offshore (m^3)")
        axs5twin.legend(["fresh centroid", "saline centroid"])
        
    return axs2, axs3, axs4, axs5, axs3twin, axs4twin, axs5twin


def plot_results(name, conc, qx, qz, times, axs0, axs1, y_labels="none", legends=False):

    pars = load_parameters(name)
    i_predev = np.argmin(np.abs(times))
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc_predev = np.flipud(conc[i_predev][:, 0, :])
    qx_predev = np.flipud(qx[i_predev][:, 0, :])
    qz_predev = np.flipud(qz[i_predev][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc_predev[i, j] == np.float32(1.e30):
                conc_predev[i,j] = np.nan

    conccm = axs0.pcolormesh(x, y, conc_predev, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx_predev = qx_predev / np.sqrt(qx_predev**2+qz_predev**2)
    qz_predev = qz_predev / np.sqrt(qx_predev**2+qz_predev**2)
    axs0.quiver(x[::80], y[::50], qx_predev[::50, ::80],-1*qz_predev[::50, ::80], 
                     color="white", width=0.002)
    axs0.set_aspect(50)
    #ax.axhline(pars.sea_levels[num*2], c='k', alpha=0.5, zorder = 3, linestyle=':', label=r"sea_level", linewidth=3)

    conc_post = np.flipud(conc[-1][:, 0, :])
    qx_post = np.flipud(qx[-1][:, 0, :])
    qz_post = np.flipud(qz[-1][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc_post[i, j] == np.float32(1.e30):
                conc_post[i,j] = np.nan

    conccm = axs1.pcolormesh(x, y, conc_post, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx_post = qx_post / np.sqrt(qx_post**2+qz_post**2)
    qz_post = qz_post / np.sqrt(qx_post**2+qz_post**2)
    axs1.quiver(x[::80], y[::50], qx_post[::50, ::80],-1*qz_post[::50, ::80], 
                     color="white", width=0.002)
    axs1.set_aspect(50)
    #ax.axhline(pars.sea_level

    if y_labels=="left":
        axs0.set_ylabel("Distance above present sea level (m)")
        
    axs1.set_xlabel("Distance offshore (m)")

    return axs0, axs1
 
def plot_set(inputs):
    
    set, set_name=inputs
    f, axs = plt.subplots(6, len(set)+1, sharey="row", figsize=(20,10))
    mid =  int(np.ceil(0.5*(1+len(set))))-1

    f.suptitle(set_name)
    
    d = {}

    for i in range(len(set)+1):
        d[f"axs3_{i}twin"] = axs[3][i].twinx()
        d[f"axs4_{i}twin"] = axs[4][i].twinx()
        d[f"axs5_{i}twin"] = axs[5][i].twinx()
    
    for i in range(len(set)):
        d["axs3_0twin"].get_shared_y_axes().join(d["axs3_0twin"], d[f"axs3_{i+1}twin"])


    conc, _, qx, qz, times = get_results("a0")
    axs[0][mid], axs[1][mid] = plot_results("a0", conc, qx, qz, times, axs[0][mid], axs[1][mid])
    axs[2][mid], axs[3][mid], axs[4][mid], axs[5][mid], d[f"axs3_{mid}twin"], d[f"axs5_{mid}twin"],d[f"axs5_{mid}twin"] \
        = plot_metrics("a0", times, axs[2][mid], axs[3][mid], axs[4][mid], axs[5][mid], d[f"axs3_{mid}twin"], d[f"axs4_{mid}twin"],d[f"axs5_{mid}twin"])

    i = 0
    y_labels= "left"
    for r in set:
        if i >= mid:
           i=i+1
        conc, _, qx, qz, times = get_results(r)
        # mf6r.metrics(r, conc, qz)
        axs[0][i], axs[1][i] = plot_results(r, conc, qx, qz, times, axs[0][i], axs[1][i], y_labels=y_labels)
        axs[2][i], axs[3][i], axs[4][i], axs[5][i], d[f"axs3_{i}twin"], d[f"axs4_{i}twin"],d[f"axs5_{i}twin"]  \
            = plot_metrics(r, times, axs[2][i], axs[3][i], axs[4][i], axs[5][i], d[f"axs3_{i}twin"], d[f"axs4_{i}twin"],d[f"axs5_{i}twin"], y_labels=y_labels)
        i=i+1
        if i < len(set)-2:
            y_labels="none"
        else:
            y_labels="right"
    
    f.tight_layout()

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
    # plot_set([sets[0], set_names[0]])
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
