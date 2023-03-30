from pars import ModelParameters, load_parameters
import numpy as np
import sloping_confined_aquifer as sca
import matplotlib.pyplot as plt
import plot_helpers as plth
import os
from matplotlib import animation
import post_processing as proc
import pandas as pd
import flopy.utils.binaryfile as bf
from matplotlib.colors import LinearSegmentedColormap
from mf6 import get_results

def metrics(name, conc, qz):

    if not os.path.exists(f"/home/ccleary/shared/results/{name}"):
        os.makedirs(f"/home/ccleary/shared/results/{name}")
    pars = load_parameters(name)
    fresh_list = proc.fresh_volume(conc, pars)
    toe_list, tip_list = proc.interface(conc, pars)
    volume_list, centroid_list = proc.mixing(conc, pars)
    fresh_flux_list, fresh_centroid_list, sal_flux_list, sal_centroid_list = proc.sgd(conc, qz, pars)
    dict = {"fresh_list": fresh_list, "toe_list": toe_list, "tip_list": tip_list,
            "volume_list": volume_list, "centroid_list": centroid_list, 
            "fresh_flux_list": fresh_flux_list, "fresh_centroid_list": fresh_centroid_list,
            "sal_flux_list": sal_flux_list, "sal_centroid_list": sal_centroid_list}
    df = pd.DataFrame(dict)
    df.to_csv(f"/home/ccleary/shared/results/{name}/metrics.csv")

def plot_metrics(name, times, plot=True):
    if not os.path.exists(f"/home/ccleary/shared/results/{name}"):
        os.makedirs(f"/home/ccleary/shared/results/{name}")

    pars = load_parameters(name)
    df = pd.read_csv(f"/home/ccleary/shared/results/{name}/metrics.csv")
    f, axs = plt.subplots(4, 1, figsize=(9,9), sharex=True)
    
    # plot fresh
    axs[0].plot(times[0::5], df.fresh_list)
    axs[0].set_ylabel("fresh volume (m^3)")
    axs[0].legend(["fresh volume"])
    axs1twin = axs[1].twinx()
    axs1twin.plot(times[0::5], df.toe_list)
    axs1twin.set_ylabel("distance offshore (m)")
    axs1twin.plot(times[0::5], df.tip_list)
    axs1twin.legend(["toe", "tip"])
    axs[2].plot(times[0::5], df.volume_list)
    axs[2].set_ylabel("volume (m^3)")
    axs[2].legend(["mixing zone"])
    axs2twin = axs[2].twinx()
    axs2twin.set_ylabel("distance offshore (m)")
    axs2twin.plot(times[0::5], df.centroid_list, linestyle="--")
    axs2twin.legend(["centoid"])
    axs[2].legend(["mixing zone"])
    axs[3].plot(times[0::5], df.fresh_flux_list)
    axs[3].plot(times[0::5], df.sal_flux_list)
    axs[3].set_xlabel("time from present (days)")
    axs[3].set_ylabel("flux (m^3)")
    axs[3].legend(["fresh sgd", "saline sgd"])
    axs3twin = axs[3].twinx() 
    axs3twin.set_ylabel("distance offshore (m^3)")
    axs3twin.plot(times[0::5], df.fresh_centroid_list, linestyle="--")
    axs3twin.plot(times[0::5], df.sal_centroid_list, linestyle="--")
    axs3twin.legend(["fresh centroid", "saline centroid"])
    if plot:
        plt.savefig(f"/home/ccleary/shared/results/{name}/metrics", dpi=300)
    else:
        return axs[0], axs[1], axs[2], axs[3], axs1twin, axs2twin, axs3twin


def plot_results(name, conc, qx, qz, times, plot=True):
    if not os.path.exists(f"/home/ccleary/shared/results/{name}"):
        os.makedirs(f"/home/ccleary/shared/results/{name}")
    
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
            if conc_predev[i, j] == np.float32(1.e30):
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
            if conc_post[i, j] == np.float32(1.e30):
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
        plt.savefig(f"/home/ccleary/shared/results/{name}/concentrations", dpi=300)
    else:
        return axs[0], axs[1]


def all_results(name):
    conc, _, qx, qz, times = get_results(name)
    metrics(name, conc, qz)
    plot_metrics(name, times)
    plot_results(name, conc, qx, qz, times)


if __name__=="__main__":
    name = input("name:")
    all_results(name)
