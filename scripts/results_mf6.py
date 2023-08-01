from pars import ModelParameters, load_parameters
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation
import post_processing as proc
import pandas as pd
import flopy.utils.binaryfile as bf
from matplotlib.colors import LinearSegmentedColormap
from mf6 import get_results
from matplotlib import animation
import pdb

def metrics(name, conc, qz):
    pars = load_parameters(name)
    if not os.path.exists(f"/home/superuser/results/{name}"):
        os.makedirs(f"/home/superuser/results/{name}")
    fresh_list = proc.fresh_volume(conc, pars)
    toe_list, tip_list = proc.interface(conc, pars)
    volume_list, centroid_list = proc.mixing(conc, pars)
    fresh_flux_list, fresh_centroid_list, sal_flux_list, sal_centroid_list = proc.sgd(conc, qz, pars)
    dict = {"fresh_list": fresh_list, "toe_list": toe_list, "tip_list": tip_list,
            "volume_list": volume_list, "centroid_list": centroid_list, 
            "fresh_flux_list": fresh_flux_list, "fresh_centroid_list": fresh_centroid_list,
            "sal_flux_list": sal_flux_list, "sal_centroid_list": sal_centroid_list}
    df = pd.DataFrame(dict)
    df.to_csv(f"/home/superuser/results/{name}/metrics.csv")

def plot_metrics(name, times, plot=True):
    if not os.path.exists(f"/home/superuser/results/{name}"):
        os.makedirs(f"/home/superuser/results/{name}")

    pars = load_parameters(name)
    df = pd.read_csv(f"/home/superuser/results/{name}/metrics.csv")
    f, axs = plt.subplots(4, 1, figsize=(9,9), sharex=True)
    
    # plot fresh
    axs[0].plot(times, df.fresh_list)
    axs[0].set_ylabel("fresh volume (m^3)")
    axs[0].legend(["fresh volume"])
    axs1twin = axs[1].twinx()
    axs1twin.plot(times, df.toe_list)
    axs1twin.set_ylabel("distance offshore (m)")
    axs1twin.plot(times, df.tip_list)
    axs1twin.legend(["toe", "tip"])
    axs[2].plot(times, df.volume_list)
    axs[2].set_ylabel("volume (m^3)")
    axs[2].legend(["mixing zone"])
    axs2twin = axs[2].twinx()
    axs2twin.set_ylabel("distance offshore (m)")
    axs2twin.plot(times, df.centroid_list, linestyle="--")
    axs2twin.legend(["centoid"])
    axs[2].legend(["mixing zone"])
    axs[3].plot(times, df.fresh_flux_list)
    axs[3].plot(times, df.sal_flux_list)
    axs[3].set_xlabel("time from present (days)")
    axs[3].set_ylabel("flux (m^3)")
    axs[3].legend(["fresh sgd", "saline sgd"])
    axs[3].set_yscale("log")
    axs3twin = axs[3].twinx() 
    axs3twin.set_ylabel("distance offshore (m^3)")
    axs3twin.plot(times, df.fresh_centroid_list, linestyle="--")
    axs3twin.plot(times, df.sal_centroid_list, linestyle="--")
    axs3twin.legend(["fresh centroid", "saline centroid"])
    if plot:
        plt.savefig(f"/home/superuser/results/{name}/metrics", dpi=300)
    else:
        return axs[0], axs[1], axs[2], axs[3], axs1twin, axs2twin, axs3twin


def plot_results(name, conc, qx, qz, times, plot=True):
    if not os.path.exists(f"/home/superuser/results/{name}"):
        os.makedirs(f"/home/superuser/results/{name}")
    
    pars = load_parameters(name)
    f,axs = plt.subplots(2,1, figsize=(6,6))

    i_predev = np.argmin(np.abs(times))
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc_predev = np.flipud(conc[i_predev][:, 0, :])
    qx_predev = np.flipud(qx[i_predev][:, 0, :])
    qz_predev = np.flipud(qz[i_predev][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_predev[i, j]) == 1.e30:
                conc_predev[i,j] = np.nan

    conccm = axs[0].pcolormesh(x, y, conc_predev, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx_predev = qx_predev / np.sqrt(qx_predev**2+qz_predev**2)
    qz_predev = qz_predev / np.sqrt(qx_predev**2+qz_predev**2)
    axs[0].quiver(x[::40], y[::5], qx_predev[::5, ::40],-1*qz_predev[::5, ::40], 
                     color="white", width=0.002)
    axs[0].set_aspect(100)
    #ax.axhline(pars.sea_levels[num*2], c='k', alpha=0.5, zorder = 3, linestyle=':', label=r"sea_level", linewidth=3)

    conc_post = np.flipud(conc[-1][:, 0, :])
    qx_post = np.flipud(qx[-1][:, 0, :])
    qz_post = np.flipud(qz[-1][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_post[i, j]) == 1.e30:
                conc_post[i,j] = np.nan

    conccm = axs[1].pcolormesh(x, y, conc_post, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx_post = qx_post / np.sqrt(qx_post**2+qz_post**2)
    qz_post = qz_post / np.sqrt(qx_post**2+qz_post**2)
    axs[1].quiver(x[::40], y[::5], qx_post[::5, ::40],-1*qz_post[::5, ::40], 
                     color="white", width=0.002)
    axs[1].set_aspect(100)
    #ax.axhline(pars.sea_level

    axs[0].set_title(f"predevelopment salinity")
    axs[1].set_title(f"postdevelopment salinity")
    axs[0].set_ylabel("Distance above present sea level (m)")
    axs[1].set_xlabel("Distance offshore (km)")
    
    if plot:
        plt.savefig(f"/home/superuser/results/{name}/concentrations", dpi=300)
    else:
        return axs[0], axs[1]


def all_results(name):
    conc, _, qx, qz, times = get_results(name)
    metrics(name, conc, qz)
    plot_metrics(name, times)
    plot_results(name, conc, qx, qz, times)

def animate_func(num, ax, pars, concentration, head):
    ax.clear()
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc = np.flipud(concentration[num][:, 0, :])
    h = np.flipud(head[num][:, 0, :])
	
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc[i, j] == np.float32(1.e30) or h[i,j] <= np.float32(-5.e2):
                if not conc[i, j] == 35.0:
                    conc[i,j] = np.nan

    conccm = ax.pcolormesh(x, y, conc, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #headc = ax.contour(x, y, h, cmap="viridis", vmax=pars.Lz-pars.z0, vmin=0)
    #ax.clabel(headc, headc.levels, inline=True, fontsize=10)
    ax.set_aspect(100)
    return ax

def plot_gif(name):
    pars = load_parameters(name)
    concentration,head, times= get_results(name)
    qx_plot = []
    qz_plot = []
    head_plot = []
    concentration_plot = []
    N=0
    
    for i, (c, h) in enumerate(zip(concentration, head)):
        if i % 10 == 0:
            print("plotting")
            concentration_plot.append(c)
            head_plot.append(h)
            # qx_plot.append(x)
            # qz_plot.append(z)
            # head_plot.append(h)
            N += 1
	
    fig = plt.figure(figsize =(14, 9))
    ax = plt.axes()
    ani = animation.FuncAnimation(fig, animate_func, fargs = (ax, pars, concentration_plot, head_plot), interval=100,   
                                   frames=N)
    f = f"animate_func_{name}.gif"
    writergif = animation.PillowWriter(fps=10)#N/6)
    ani.save(f, writer=writergif)
    #plt.show()
def plot_gif(name):
    pass


def save_results_for_saltelli(name, ensemble_name):
    pars = load_parameters(name)
    conc, head, qx, qz, times = get_results(name)
    volume_list, centroid_list = proc.mixing(conc, pars)
    fresh_list = proc.fresh_volume(conc, pars)
    toe_list, tip_list = proc.interface(conc, pars)
    metrics =np.stack((volume_list, centroid_list, fresh_list, toe_list, tip_list))
    if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/"):
        os.makedirs(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/")
     
    pre_dev = [conc[1], head[1], qx[1], qz[1]]
    post_dev = [conc[-1], head[-1], qx[-1], qz[-1]]
    results = [pre_dev, post_dev]
    
    np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}_metrics.csv", metrics, delimiter=",")
    np.save(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}_results.npy", results)

def plot_color_mesh_saltelli(name, n):
    for i, in range(1):
        pars = load_parameters(name)
        x = np.linspace(-pars.x_b, pars.L, pars.ncol)
        y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_results.npy")
        f, axs = plt.subplots(2, 2)
        axs[0, 0].pcolormesh(results[0, 0][:, 0, :], cmap="viridis", vmax=35.7, vmin=0)
        axs[0, 1].pcolormesh(results[1, 0][:, 0, :], cmap="viridis", vmax=35.7, vmin=0)
        axs[0, 0].quiver(results[0, 2][:, 0, :], results[0, 3][:, 0, :])
        axs[0, 1].quiver(results[1, 2][:, 0, :], results[1, 3][:, 0, :])
        axs[0, 0].set_box_aspect(0.1)
        axs[0, 1].set_box_aspect(0.1)
        axs[1, 0].set_box_aspect(0.1)
        axs[1, 1].set_box_aspect(0.1)
        axs[1, 0].contour(results[0, 1][:, 0, :])
        axs[1, 1].contour(results[0, 1][:, 0, :])
        plt.savefig(f"/home/superuser/sloping_confined_aquifer/results/{name}/plots/{i}.png")

if __name__=="__main__":
    name = input("name:")
    all_results(name)
