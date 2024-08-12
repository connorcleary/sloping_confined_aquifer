from pars import ModelParameters, load_parameters
import numpy as np
import scripts.old_scripts.sloping_confined_aquifer as sca
import matplotlib.pyplot as plt
import scripts.old_scripts.plot_helpers as plth
import os
from matplotlib import animation
import post_processing as proc
import pandas as pd
import flopy.utils.binaryfile as bf
from matplotlib.colors import LinearSegmentedColormap

def metrics(name):
    pars = load_parameters(name)
    conc, _, _, _, qz = sca.load_results(name)
    fresh_list = proc.fresh_volume(conc, pars)
    toe_list, tip_list = proc.interface(conc, pars)
    volume_list, centroid_list = proc.mixing(conc, pars)
    fresh_flux_list, fresh_centroid_list, sal_flux_list, sal_centroid_list = proc.sgd(conc, qz, pars)
    dict = {"fresh_list": fresh_list, "toe_list": toe_list, "tip_list": tip_list,
            "volume_list": volume_list, "centroid_list": centroid_list, 
            "fresh_flux_list": fresh_flux_list, "fresh_centroid_list": fresh_centroid_list,
            "sal_flux_list": sal_flux_list, "sal_centroid_list": sal_centroid_list}
    df = pd.DataFrame(dict)
    df.to_csv(f"results/{name}/metrics.csv")

def plot_metrics(name):
    pars = load_parameters(name)
    df = pd.read_csv(f"results/{name}/metrics.csv")
    f, axs = plt.subplots(4, 1, figsize=(9,9), sharex=True)
    
    # plot fresh
    axs[0].plot(pars.times, df.fresh_list)
    axs[0].set_ylabel("fresh volume (m^3)")
    axs[0].legend(["fresh volume"])
    axs1twin = axs[1].twinx()
    axs1twin.plot(pars.times, df.toe_list)
    axs1twin.set_ylabel("distance offshore (m)")
    axs1twin.plot(pars.times, df.tip_list)
    axs1twin.legend(["toe", "tip"])
    axs[2].plot(pars.times, df.volume_list)
    axs[2].set_ylabel("volume (m^3)")
    axs[2].legend(["mixing zone"])
    axs2twin = axs[2].twinx()
    axs2twin.set_ylabel("distance offshore (m)")
    axs2twin.plot(pars.times, df.centroid_list, linestyle="--")
    axs2twin.legend(["centoid"])
    axs[2].legend(["mixing zone"])
    axs[3].plot(pars.times, df.fresh_flux_list)
    axs[3].plot(pars.times, df.sal_flux_list)
    axs[3].set_xlabel("time from present (days)")
    axs[3].set_ylabel("flux (m^3)")
    axs[3].legend(["fresh sgd", "saline sgd"])
    axs3twin = axs[3].twinx() 
    axs3twin.set_ylabel("distance offshore (m^3)")
    axs3twin.plot(pars.times, df.fresh_centroid_list, linestyle="--")
    axs3twin.plot(pars.times, df.sal_centroid_list, linestyle="--")
    axs3twin.legend(["fresh centroid", "saline centroid"])
    plt.savefig(f"results/{name}/metrics", dpi=300)

def plot_results(name):
    pars = load_parameters(name)
    f,axs = plt.subplots(2,1, figsize=(6,6))
    conc, head,qx,_,qz = sca.load_results(name)
    i_predev = np.where(np.array((pars.times), dtype=int) == 0)[0][0] -1
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

    plt.savefig(f"results/{name}/concentrations", dpi=300)


def animate_func(num, ax, pars, concentration, qx, qz, head):
    ax.clear()
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc = np.flipud(concentration[num][:, 0, :])
    qx = np.flipud(qx[num][:, 0, :])
    qz = np.flipud(qz[num][:, 0, :])
    head = np.flipud(head[num][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc[i, j] == np.float32(1.e30) or head[i,j] <= np.float32(-5.e2):
                conc[i,j] = np.nan

    conccm = ax.pcolormesh(x, y, conc, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx = qx / np.sqrt(qx**2+qz**2)
    qz = qz / np.sqrt(qx**2+qz**2)
    ax.quiver(x[::40], y[::5], qx[::5, ::40],-1*qz[::5, ::40], 
                     color="white", width=0.002)
    # hc = ax.contour(x, y, head, cmap="viridis", vmax=1, vmin=-1, levels = np.linspace(-10, 1, 10))
    # ax.clabel(hc, hc.levels, inline=True, fontsize=10)
    # conccb = plt.colorbar(conccm, shrink=1, ax=ax)
    # conccb.ax.set_title('Salinity (kg/m^3)', fontsize = 'small')

    ax.set_aspect(100)
    #ax.axhline(pars.sea_levels[num*2], c='k', alpha=0.5, zorder = 3, linestyle=':', label=r"sea_level", linewidth=3)
    if pars.times[num*5] <= 0:
        time = "BP"
    else:
        time = "AP"
    
    ax.set_title(f"salinity at {int(np.abs(pars.times[num*5]/365))} years {time}")
    ax.set_ylabel("Distance above present sea level (m)")
    ax.set_xlabel("Distance offshore (km)")
    return ax
    

def plot_gif(name):
    pars = load_parameters(name)
    concentration,head,qx,_,qz = sca.load_results(name)
    qx_plot = []
    qz_plot = []
    head_plot = []
    concentration_plot = []
    N= 0
    for i, (c, x, z, h) in enumerate(zip(concentration, qx, qz, head)):
        if i % 1 == 0:
            print("plotting")
            concentration_plot.append(c)
            qx_plot.append(x)
            qz_plot.append(z)
            head_plot.append(h)
            N += 1
    

    fig = plt.figure(figsize =(14, 9))
    ax = plt.axes()
    ani = animation.FuncAnimation(fig, animate_func, fargs = (ax, pars, concentration_plot, qx_plot, qz_plot, head_plot), interval=100,   
                                   frames=N)
    f = f"animate_func_{name}.gif"
    writergif = animation.PillowWriter(fps=10)#N/6)
    ani.save(f, writer=writergif)
    #plt.show()


def plot_last_step(name):

    pars = load_parameters(name)
    name = pars.name
    model_ws = f"./model_files/{name}"

    # open binary files
    ucnobj = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws, f'{name}.hds'))
 
    # get head and concentration data
    concentration = ucnobj.get_alldata()[-1]
    head = headobj.get_alldata()[-1]

    qx = np.zeros_like(concentration)
    qz = np.zeros_like(concentration)
    

    qx = cbbobj.get_data(text="flow right face")[-1]
    qz = cbbobj.get_data(text="flow lower face")[-1]

    f, axs = plt.subplots(3,1)
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc = np.flipud(concentration[:, 0, :])
    qx = np.flipud(qx[:, 0, :])
    qz = np.flipud(qz[:, 0, :])
    head = np.flipud(head[:, 0, :])
    courant = np.zeros_like(conc)

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc[i, j] == np.float32(1.e30) or head[i,j] <= np.float32(-5.e2):
                conc[i,j] = np.nan
                courant[i, j] = "0"
            elif abs(qx[i, j]*pars.dt_modern) >= pars.dx:
                if abs(qz[i, j]*pars.dt_modern) >= pars.dz:
                    courant[i,j] = "4"
                else:
                    courant[i,j] = "2"
            elif abs(qz[i, j]*pars.dt_modern) >= pars.dz:
                courant[i,j] = "3"
            else:
                courant[i,j] = "1"



    conccm = axs[0].pcolormesh(x, y, conc, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    axs[0].quiver(x[::40], y[::5], qx[::5, ::40],-1*qz[::5, ::40], 
                     color="white", width=0.002)
    hc = axs[1].contourf(x, y, head, cmap="viridis", levels = np.linspace(-1, 1, 11))
    plt.colorbar(hc, shrink=1, ax=axs[1])
    #axs[1].clabel(hc, hc.levels, inline=True, fontsize=10)
    # conccb = plt.colorbar(conccm, shrink=1, ax=ax)
    # conccb.ax.set_title('Salinity (kg/m^3)', fontsize = 'small')

    cmap = LinearSegmentedColormap.from_list("courant", [(1.0, 1.0, 1.0), (0.0, 0.93, 0.0), (1.0, 1.0, 0.0), (1.0, 0.65, 0.0), (1.0, 0.0, 0.0)])
    axs[2].pcolor(x, y, courant, cmap=cmap)

    axs[0].set_aspect(100)
    axs[1].set_aspect(100)
    axs[2].set_aspect(100)

    plt.show()


    


