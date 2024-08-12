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
from get_cells import get_cells
import matplotlib.colors as colors 


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
    pars = load_parameters(name, backup=True)
    concentration,head, times= get_results(name, backup=True)
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


def save_results_for_saltelli(name, ensemble_name):
    pars = load_parameters(name)
    conc, head, qx, qz, times = get_results(name)
    volume_list, centroid_list = proc.mixing(conc, pars)
    fresh_list = proc.fresh_volume(conc, pars)
    toe_list, tip_list = proc.interface(conc, pars)
    flux_list = proc.abstracted_flux(conc, qx, pars)
    dhh, dhv, x_sal, z_sal = proc.time_constants(conc, head, pars)

    metrics =np.stack((volume_list, centroid_list, fresh_list, toe_list, tip_list, flux_list, dhh, dhv, x_sal, z_sal))
    if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/"):
        os.makedirs(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/")
     
    pre_dev = [conc[1], head[1], qx[1], qz[1]]
    post_dev = [conc[-1], head[-1], qx[-1], qz[-1]]
    results = [pre_dev, post_dev]
    
    np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}_metrics.csv", metrics, delimiter=",")
    np.save(f"/home/superuser/sloping_confined_aquifer/results/{ensemble_name}/{name}_results.npy", results)

def plot_color_mesh_saltelli(name, n):
    for m in range(n):
        try:
            pars = load_parameters(f"{name}{m}", backup=False)
            x = np.linspace(-pars.x_b, pars.L, pars.ncol)/1000
            y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
            results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{m}_results.npy")
            f, axs = plt.subplots(2, 2, figsize=(8.3, 3), layout="constrained", sharex=True, sharey=True)
            plt.rcParams.update({'font.size': 9})
            for i in range(pars.nlay):
                for j in range(pars.ncol):
                    if abs(results[0, 0][i, 0, j]) == 1.e30:
                        results[0, 0][i, 0, j] = np.nan
                        results[0, 1][i, 0, j] = np.nan
                    if abs(results[1, 0][i, 0, j]) == 1.e30:
                        results[1, 0][i, 0, j] = np.nan
                        results[1, 1][i, 0, j] = np.nan
            hf_pre = proc.hf(results[0, 0][:, 0, :], results[0, 1][:, 0, :], pars)
            hf_post = proc.hf(results[1, 0][:, 0, :], results[1, 1][:, 0, :], pars)
            qx_pre = results[0, 2][::-5, 0, ::50] / np.sqrt(results[0, 2][::-5, 0, ::50]**2+np.abs(results[0, 3][::-5, 0, ::50])**2)
            qz_pre = results[0, 3][::-5, 0, ::50] / np.sqrt(results[0, 2][::-5, 0, ::50]**2+np.abs(results[0, 3][::-5, 0, ::50])**2)
            qx_post = results[1, 2][::-5, 0, ::50] / np.sqrt(results[1, 2][::-5, 0, ::50]**2+np.abs(results[1, 3][::-5, 0, ::50])**2)
            qz_post = results[1, 3][::-5, 0, ::50] / np.sqrt(results[1, 2][::-5, 0, ::50]**2+np.abs(results[1, 3][::-5, 0, ::50])**2)
            axs[0, 0].pcolormesh(x, y, results[0, 0][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
            cm = axs[0, 1].pcolormesh(x, y, results[1, 0][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
            axs[0, 0].quiver(x[::50], y[::5], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
            axs[0, 1].quiver(x[::50], y[::5], qx_post/2, qz_post/2, color="white", width=0.002, scale=20)
            axs[0, 0].set_box_aspect(0.3)
            axs[0, 1].set_box_aspect(0.3)
            axs[1, 0].set_box_aspect(0.3)
            axs[1, 1].set_box_aspect(0.3)
            head_min = np.nanmin(hf_post)
            head_max = np.nanmax(hf_pre)
            cp0 = axs[1, 0].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
            cp = axs[1, 1].contour(x, y, hf_post[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
            axs[0,0].set_ylabel("Depth [masl]")
            axs[1,0].set_ylabel("Depth [masl]")
            axs[1,0].set_xlabel("Distance offshore [km]")
            axs[1,1].set_xlabel("Distance offshore [km]")
            cb = plt.colorbar(cm, shrink=1)
            cb.ax.set_title('C [PSU]')
            axs[1,0].clabel(cp0, inline=True, fontsize=7)
            axs[1,1].clabel(cp, inline=True, fontsize=7)
            axs[0, 0].set_title("Predevelopment")
            axs[0, 1].set_title("Postdevelopment")

            axs[0,0].annotate("Concentration \nand flux", xy=(0, 0.5), xytext=(-axs[0,0].yaxis.labelpad - 5, 0),
                xycoords=axs[0,0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center', ma='center')
            
            axs[1,0].annotate("Head", xy=(0, 0.5), xytext=(-axs[1,0].yaxis.labelpad - 30, 0),
                xycoords=axs[1,0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center', ma='center')
            
            if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/results/{name}/plots/"):
                os.makedirs(f"/home/superuser/sloping_confined_aquifer/results/{name}/plots/")

            plt.savefig(f"/home/superuser/sloping_confined_aquifer/results/{name}/plots/{m}.png", dpi=600)
            
        except:
            pass
        # plot_gif(f"{name}{m}")

def extract_pre_post_volume(name, n):
    total_og = []
    total_change = []
    for m in range(n):
        pars = load_parameters(f"{name}{m}", backup=False)
        if name==("paper3"):
            pars.dx=25
            pars.ncol=1360
        if pars.bottom_cells == None:
            aquitard_cells, aquifer_cells, _, bottom_cells, _, _, _ = get_cells(pars)
            pars.aquitard_cells= aquitard_cells
            pars.aquifer_cells = aquifer_cells
            pars.bottom_cells = bottom_cells
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{m}_results.npy")
        concs = [results[0][0], results[1][0]]
        fresh_list = proc.fresh_volume_total(concs, pars)
        total_change.append(fresh_list[0]-fresh_list[1])
        total_og.append(fresh_list[0])


    return total_og, total_change

def paleo_volumes(name, n):
    paleo_volumes = []
    
    for m in range(n):
        pars = load_parameters(f"{name}{m}", backup=True)
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{m}_results.npy")
        concs = [results[0][0], results[1][0]]
        qxs = [results[0][2], results[1][2]]
        x_paleo = 0
        if name==("paper3"):
            pars.dx=25
            pars.ncol=1360
        if pars.bottom_cells == None:
            aquitard_cells, aquifer_cells, _, bottom_cells, _, _, _ = get_cells(pars)
            pars.aquitard_cells= aquitard_cells
            pars.aquifer_cells = aquifer_cells
            pars.bottom_cells = bottom_cells
        

        for cell in pars.bottom_cells:
            if qxs[0][cell[0], cell[1], cell[2]] < 0:
                x_paleo = cell[2]
                break
        if x_paleo == 0:
            paleo_volumes.append([0, 0])
        else:
            pv = [0, 0]
            for i in range(2):               
                for cell in pars.aquifer_cells:
                    if cell[2]>=x_paleo and concs[i][cell[0], cell[1], cell[2]] <= 0.35:
                        pv[i] += 1
                pv[i] = pv[i]*pars.dx*pars.dz

            paleo_volumes.append(pv)

    return paleo_volumes

def toe_position(name, n):
    toe = []
    for m in range(n):
        pars = load_parameters(f"{name}{m}", backup=True)
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{m}_results.npy")
        concs = [results[0][0], results[1][0]]
        toe.append(proc.interface(concs, pars)[0])
    return toe

def get_pre_concs(name, n):
    concs = []
    for m in range(n):
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{m}_results.npy")
        concs.append(results[0, 0][:, 0, :])
    return concs 

def get_pars(name, n):
    inputs = []
    for i in range(n):
        pars = load_parameters(f"{name}{i}", backup=True)
        # inputs.append([pars.H, pars.D, pars.K_aquifer, pars.K_aquitard, pars.alpha_L, pars.anis_aquifer, pars.beta, pars.L, pars.h_modern, pars.porosity])
        # inputs.append([pars.K_aquifer, pars.K_aquitard, pars.alpha_L, pars.anis_aquifer, pars.h_modern, pars.porosity, pars.gradient])
        inputs.append([pars.K_aquifer, pars.K_aquitard/pars.anis_aquitard, pars.alpha_L, pars.anis_aquifer, pars.h_modern, pars.porosity])
    return np.array(inputs)

def get_outputs(name, n, inputs):
    metrics = []
    for i in range(n):
        metrics.append(np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_metrics.csv", delimiter=","))
    
    fractional_volume_change = []
    mixing = []
    total_volume_change = []
    original_volume = []
    total_flux = []
    predev_expected_flux = []
    for i in range(n):
        mixing.append((max([metrics[i][0][-1]-metrics[i][0][1],0]))/(metrics[i][2][1]-metrics[i][2][-1]))
        fractional_volume_change.append((metrics[i][2][1]-metrics[i][2][-1])/metrics[i][2][1])
        total_volume_change.append((metrics[i][2][1]-metrics[i][2][-1])*inputs[i][-1])
        original_volume.append((metrics[i][2][1])*inputs[i][-1])
        total_flux.append(sum(metrics[i][5][2:])*500) # 500 is the multiplying this by the number of days between observations
        if total_flux[-1] > 0: 
            pass
        predev_expected_flux.append(metrics[i][5][1]*36500) # multiplying by the number of dates 

    return fractional_volume_change, total_volume_change, mixing, original_volume, total_flux, predev_expected_flux

def get_outputs_for_timescales(name, n, inputs):
    dhh = [] # horizontal head change
    dl0 = [] # distance to salt in aquifer
    dd0 = [] # depth from aquifer to saline groundwater in aquitard
    dl1 = [] # for end of post development phase
    dd1 = [] # for end of post development phase
    dhv = [] # vertical head change
    T = []
    T_test = []
    for i in range(n):
        metrics = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_metrics.csv", delimiter=",")
        results = np.load(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_results.npy")
        pars = load_parameters(f"{name}{i}", backup=True)
        shoreline = pars.bottom_cells[int(pars.x_b/pars.dx)]
        toe_cell = pars.bottom_cells[int((pars.x_b+metrics[3,1])/pars.dx)]
        # find head at toe 
        dhh.append(results[1, 1][shoreline[0], shoreline[1], shoreline[2]]
                   -results[1, 1][toe_cell[0], toe_cell[1], toe_cell[2]])
        dl0.append(metrics[3,1])
        dl1.append(metrics[3,-1])
        dd0.append(pars.D - ((pars.nlay - np.argmax(np.fliplr(np.array([(36>results[0, 0][:, 0, int(pars.x_b/pars.dx)+1]) & (results[0, 0][:, 0, int(pars.x_b/pars.dx)+1]>0.35)]))))*pars.dz-(pars.Lz-pars.z0)))
        if dd0[-1] == -35.0:
            dd0[-1] = pars.D
        dd1.append(pars.D - ((pars.nlay - np.argmax(np.fliplr(np.array([(36>results[1, 0][:, 0, int(pars.x_b/pars.dx)+1]) & (results[1, 0][:, 0, int(pars.x_b/pars.dx)+1]>0.35)]))))*pars.dz-(pars.Lz-pars.z0)))
        dhv.append(results[1,1][int(shoreline[0]-pars.H/pars.dz), shoreline[1], shoreline[2]]
                   - results[1,1][int(((pars.nlay - np.argmax(np.fliplr(np.array([(36>results[0, 0][:, 0, int(pars.x_b/pars.dx)+1]) & (results[0, 0][:, 0, int(pars.x_b/pars.dx)+1]>0.35)]))))*pars.dz-(pars.Lz-pars.z0))/pars.dz), shoreline[1], shoreline[2]])
        # definition based on fraction salinized
        # if (metrics[2][1]-metrics[2][-1])/metrics[2][1] < 1:
        #     T.append(100/((metrics[2][1]-metrics[2][-1])/metrics[2][1]))
        # else: 
        #     T.append(500/365*np.argmax([metrics[2][1:]==0]))
        
        # definition based on distances

        # if salinization has occured
        if dl1[-1] <= 0 or dd1[-1] <= 0:
            if np.argmax([metrics[-1][1:]<=0]) != 0 and np.argmax([metrics[-2][1:]<=0]) != 0:
                T.append(np.min(np.abs([500/365*np.argmax([metrics[-1][1:]<=0]), 500/365*np.argmax([metrics[-2][1:]<=0])])))
            elif np.argmax([metrics[-1][1:]<=0]) !=0:
                T.append(500/365*np.argmax([metrics[-1][1:]<=0])) 
            elif np.argmax([metrics[-2][1:]<=0]) != 0:
                T.append(500/365*np.argmax([metrics[-2][1:]<=0]))
            else:
                T.append(np.min([100/(np.abs((dl0[-1]-dl1[-1]))/dl0[-1]), 100/((dd0[-1]-dd1[-1])/dd0[-1])]))
        # if salinization has not occured
        else: 
            T.append(np.min([100/(np.abs((dl0[-1]-dl1[-1]))/dl0[-1]), 100/((dd0[-1]-dd1[-1])/dd0[-1])]))
        
        if T[-1] == np.inf:
            T[-1] = 1e3
        
    return {"dhh": np.array(dhh), "dl0": np.array(dl0), "dl1": np.array(dl1),"dd0": np.array(dd0), "dd1": np.array(dd1),  "dhv": np.array(dhv), "T_sal": np.array(T)}

if __name__=="__main__":
    name = input("name:")
    all_results(name)
