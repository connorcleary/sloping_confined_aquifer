import numpy as np
import os
import mf6
from pars import ModelParameters, load_parameters
import results_mf6
from results_mf6 import save_results_for_saltelli, plot_color_mesh_saltelli, get_pars, get_outputs, paleo_volumes
import matplotlib.pyplot as plt
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
from results_smaller import load_smaller
import post_processing as proc
from matplotlib import cm
from matplotlib.lines import Line2D

def raleigh(conc, qx, pars, pars2):

    qxs = [qx[cell[0]][cell[1]][cell[2]] for cell in pars2.top_cells[3000:]]
    concs = [conc[cell[0]][cell[1]][cell[2]] for cell in pars2.top_cells[3000:]]
    return [np.abs((0.025)*(conci/35)*pars.alpha_L*pars.K_aquifer/pars.anis_aquifer/(1*pars.diff+np.abs(pars.alpha_L*qxi)/pars.porosity)) for qxi, conci in zip(qxs, concs)]

def back_dispersion(qz, pars, pars2):
    qzs = [qz[int(cell[0]-40)][cell[1]][cell[2]] for cell in pars2.top_cells[3000:]]
    return [(pars.alpha_L*pars.alpha_anisV*qzi)/pars.porosity/pars.diff for qzi in qzs]

# def raleigh_vert(qz, pars):
#     rv = []
#     for cell in pars.top_cells:
#         rv.append(25*pars.K_aquifer/pars.anis_aquifer/1000/qz[cell[0]][cell[1]][cell[2]]) 
#     return rv

# def raleigh_hori(qx, pars):
#     rh = []
#     for cell in pars.top_cells:
#         rh.append(25*pars.K_aquifer*pars.alpha_L/pars.anis_aquifer/1000/(pars.diff+qx[cell[0]][cell[1]][cell[2]]/pars.porosity*pars.alpha_L)) 
#     return rh

# def back_dispersion(qz, pars):
#     bd = []
#     for cell in pars.top_cells:
#         bd.append(pars.alpha_L*pars.alpha_anisV*qz[int(cell[0]-pars.D/pars.dz)][cell[1]][cell[2]]) 


# def plot_dimensionless(qx, qz, name, pars, x):
#     f, axs = plt.subplots(2, 1, figsize=(4, 4))
#     rv = raleigh_vert(qz[1], pars)
#     rh = raleigh_hori(qx[1], pars)
#     # axs[0].plot(x, rv)
#     axs[0].plot(x, rh)
#     axs[0].set_yscale("log")
#     plt.savefig(f"/home/superuser/sloping_confined_aquifer/results/base30/{Ka}_{KbV}_{h_modern}_post.png", dpi=600)


def plot_predev_bases():
    # ensemble_name="base10"
    # Kas = [1, 10]
    # KbVs = [1e-5, 1e-4, 1e-3]
    # alpha = 0.1
    # anis = 100
    # n = 0.5
    # h_modern = 1
    x_step = 1400
    y_step = 10
    
    # for Ka in Kas:
    #     for KbV in KbVs:
    #         try:
    #             f, axs = plt.subplots(1, 2, figsize=(4, 6), layout="constrained", sharex=True, sharey=True)
    #             name=f"{Ka}_{KbV}_{h_modern}"
    ensemble_name="base30"            
    Kas = [100, 10, 1]
    KbVs = [1e-3, 1e-4, 1e-5]
    alphas = [1, 10, 100]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    
    for Ka in Kas:
        for KbV in KbVs:
            f, axs = plt.subplots(1, 2, figsize=(6, 4), layout="constrained", sharex=True, sharey=True)
            name=f"10km_{Ka}_{KbV}"
            results = load_smaller(name)
            conc = results["conc"]
            head = results["head"]
            qx = results["qx"]
            qz = results["qz"]
            pars = load_parameters(f"{name}", backup=True)
            x = np.linspace(-pars.x_b, pars.L, pars.ncol)/1000
            y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
            plt.rcParams.update({'font.size': 9})
            for i in range(pars.nlay):
                for j in range(pars.ncol):
                    if abs(conc[1][i, 0, j]) == 1.e30:
                        conc[1][i, 0, j] = np.nan
                        head[1][i, 0, j] = np.nan

            hf_pre = proc.hf(conc[1][:, 0, :], head[1][:, 0, :], pars)
            qx_pre = qx[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
            qz_pre = qz[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
            axs[0].pcolormesh(x, y, conc[1][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
            axs[0].quiver(x[::x_step], y[::y_step], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
            axs[0].set_box_aspect(0.3)
            
            head_min = np.nanmin(hf_pre)
            head_max = np.nanmax(hf_pre)
            cp0 = axs[1].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
            axs[1].set_box_aspect(0.3)
            axs[0].set_ylabel("Depth [masl]")
            axs[1].clabel(cp0, inline=True, fontsize=7)
            plt.savefig(f"/home/superuser/sloping_confined_aquifer/results/base10/{name}_pre.png", dpi=600)

            # plot_dimensionless(qx[1], qz[1], name, pars)
            
        # except:
            # print(f"model {{Ka}_{KbV}_{h_modern}} not found")
    
def plot_results_1():   
    x_step = 1400
    y_step = 10
    
    # for Ka in Kas:
    #     for KbV in KbVs:
    #         try:
    #             f, axs = plt.subplots(1, 2, figsize=(4, 6), layout="constrained", sharex=True, sharey=True)
    #             name=f"{Ka}_{KbV}_{h_modern}"
    ensemble_name="base30"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    # alphas = [1, 10, 100]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    f,axs = plt.subplots(3, 2, figsize=(8.5, 5), layout="constrained", sharex=True, sharey=True)
    i = 0 
    for Ka in Kas:
        KbV = KbVs[i]  
        name=f"10km_{Ka}_{KbV}"
        results = load_smaller(name)
        conc = results["conc"]
        head = results["head"]
        qx = results["qx"]
        qz = results["qz"]
        pars = load_parameters(f"{name}", backup=True)
        x = np.linspace(-pars.x_b, pars.L, pars.ncol)/1000
        y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
        plt.rcParams.update({'font.size': 9})
        for k in range(pars.nlay):
            for j in range(pars.ncol):
                if abs(conc[1][k, 0, j]) == 1.e30:
                    conc[1][k, 0, j] = np.nan
                    head[1][k, 0, j] = np.nan

        hf_pre = proc.hf(conc[1][:, 0, :], head[1][:, 0, :], pars)
        qx_pre = qx[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
        qz_pre = qz[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
        pcm = axs[i][0].pcolormesh(x, y, conc[1][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
        axs[i][0].quiver(x[::x_step], y[::y_step], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
        axs[i][0].set_box_aspect(0.3)
        
        head_min = np.nanmin(hf_pre)
        head_max = np.nanmax(hf_pre)
        cp0 = axs[i][1].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
        axs[i][1].set_box_aspect(0.3)
        axs[i][0].set_ylabel("Depth [masl]")
        axs[i][1].clabel(cp0, inline=True, fontsize=7)
        
        i = i+1

    cb = f.colorbar(pcm, ax=axs[:, 0], shrink=0.6, location='left', ticks=[0, 0.35, 3.5, 35])
    cb.set_label("C [PSU]")
    cb.ax.set_yticklabels(["0", "0.35", "3.5", "35"])
    axs[2][0].set_xlabel("Distance offshore [km]")
    axs[2][1].set_xlabel("Distance offshore [km]")
    axs[0][0].title.set_text('Concentration and Flux')
    axs[0][1].title.set_text('Head')
    plt.show()    

def plot_results_2():
    ensemble_name="base30"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    alphas = [1, 10, 100]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
            K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
            z0=43.333, L=10000, x_b=3000,  
            h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 

    pars2 = mf6.build_steady_model(pars2, True)
    i = 0 
    for Ka in Kas: 
        for KbV in KbVs:
            i = i + 1
            name = rf"10km_{Ka}_{KbV}"
            results = load_smaller(name)
            raleighi = raleigh(results["conc"][1], results["qx"][1], load_parameters(name, backup = True), pars2)
            raleighs.append(raleighi)
            epsilon = back_dispersion(results["qz"][1], load_parameters(name, backup = True), pars2)
            epsilons.append(epsilon)
            names.append(rf"Case {i}, $K_a$={Ka}, $K_b^V$={KbV}")

    f, axs = plt.subplots(1, 2, figsize=(8.45, 4))
    epsilons = np.array(epsilons)
    raleighs = np.array(raleighs)
    eps_ave = np.mean(epsilons.reshape(4, 300, -1), axis=2)
    rals_ave = np.mean(raleighs.reshape(9, 100, -1), axis=2)
    axs[0].plot(np.linspace(0, 10, 10000), epsilons.T)
    lines = axs[1].plot(np.linspace(0, 10, 100), rals_ave.T)
    axs[1].legend(iter(lines), (names), loc='center left', bbox_to_anchor=(1, 0.5))
    # axs[0].set_yscale("log")
    axs[1].set_yscale("log")
    axs[1].set_xlabel("Distance offshore [km]")
    axs[0].set_xlabel("Distance offshore [km]")
    axs[1].set_ylabel("Ra' [-]")
    axs[0].set_ylabel(r"$\epsilon_D$ [-]")
    axs[0].set_ylim(bottom=0)
    f.tight_layout()

    
    plt.show()

def plot_results_3():
    x_step = 25
    y_step = 4
    
    # for Ka in Kas:
    #     for KbV in KbVs:
    #         try:
    #             f, axs = plt.subplots(1, 2, figsize=(4, 6), layout="constrained", sharex=True, sharey=True)
    #             name=f"{Ka}_{KbV}_{h_modern}"
    ensemble_name="base30"            
    Kas = [10, 100]
    KbVs = [1e-4]
    x_ranges=[[5,5.25],[3,3.25]]
    # alphas = [1, 10, 100]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    f,axs = plt.subplots(1, 2, figsize=(8.5, 4), layout="constrained", sharey=True)
    i = 0 
    for Ka in Kas:
        KbV = KbVs[0]  
        name=f"10km_{Ka}_{KbV}"
        results = load_smaller(name)
        conc = results["conc"]
        qx = results["qx"]
        qz = results["qz"]
        pars = load_parameters(f"{name}", backup=True)
        x = np.linspace(-pars.x_b, pars.L, pars.ncol)/1000
        y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
        plt.rcParams.update({'font.size': 9})
        for k in range(pars.nlay):
            for j in range(pars.ncol):
                if abs(conc[1][k, 0, j]) == 1.e30:
                    conc[1][k, 0, j] = np.nan

        qx_pre = qx[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
        qz_pre = qz[1][::-y_step, 0, ::x_step] / np.sqrt(qx[1][::-y_step, 0, ::x_step]**2+np.abs(qz[1][::-y_step, 0, ::x_step])**2)
        pcm = axs[i].pcolormesh(x, y, conc[1][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
        axs[i].quiver(x[::x_step], y[::y_step], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
        axs[i].set_box_aspect(1)
        axs[i].set_xlim(x_ranges[i][0], x_ranges[i][1])
        axs[i].set_ylim(-40, -20)
        axs[0].set_ylabel("Depth [masl]")
        axs[i].set_xlabel("Distance offshore [km]")
        i = i+1
        axs[0].set_title("Case 5: $K_a$=10")
        axs[1].set_title("Case 8: $K_a$=100")
    
    cb = f.colorbar(pcm, ax=axs[1], shrink=0.6, location='right', ticks=[0, 0.35, 3.5, 35])
    cb.set_label("C [PSU]")
    cb.ax.set_yticklabels(["0", "0.35", "3.5", "35"])
    plt.show()

def plot_results_4():
    ensemble_name="base30"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    alphas = [1, 10, 100]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
            K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
            z0=43.333, L=10000, x_b=3000,  
            h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 

    pars2 = mf6.build_steady_model(pars2, True)
    i = 0
    toes = []
    tips = []
    total_volumes = []
    mixing_volumes = []
    paleo_volumes = []
    f, axs = plt.subplots(2, 2, sharex=True, figsize=(8.5,6))

    for Ka in Kas:
        toe = []
        tip = []
        tv = []
        mv = []
        pv = []
        for KbV in KbVs:
            name = rf"10km_{Ka}_{KbV}"
            results = load_smaller(name)
            to, ti = proc.interface([results["conc"][1]], pars2)
            toe.append(to[0])
            tip.append(ti[0])
            tv.append(proc.fresh_volume_total([results["conc"][1]], pars2)[0])
            mv.append(proc.mixing([results["conc"][1]], pars2)[0])
            pv.append(proc.paleo_volume([results["conc"][1]], [results["qx"][1]], pars2)[0])
        toes.append(toe)
        tips.append(tip)
        total_volumes.append(tv)
        mixing_volumes.append(mv)
        paleo_volumes.append(pv)

    cmap = cm.get_cmap("Greens")
    axs[0, 0].plot(Kas, np.array(total_volumes).T[0, :], c=cmap(0.25), marker="x")
    axs[0, 0].plot(Kas, np.array(total_volumes).T[1, :], c=cmap(0.5), marker="x")
    axs[0, 0].plot(Kas, np.array(total_volumes).T[2, :], c=cmap(0.75), marker="x")
    axs[0, 0].set_xscale("log")
    axs[0, 0].set_ylabel(r"$V_{ofg}^{pre}$ [$m^3$]")
    

    axs[0, 1].plot(Kas, np.array(toes).T[0, :], c=cmap(0.25), marker="x")
    axs[0,1].plot(Kas, np.array(toes).T[1, :], c=cmap(0.5), marker="x")
    axs[0,1].plot(Kas, np.array(toes).T[2, :], c=cmap(0.75), marker="x")
    axs[0, 1].plot(Kas, np.array(tips).T[0, :], c=cmap(0.25), ls="--", marker=".")
    axs[0,1].plot(Kas, np.array(tips).T[1, :], c=cmap(0.5), ls="--", marker=".")
    axs[0,1].plot(Kas, np.array(tips).T[2, :], c=cmap(0.75), ls="--", marker=".")
    axs[0,1].set_xscale("log")
    axs[0, 1].set_ylabel(r"$x^{pre}_{toe}, x^{pre}_{tip}$ [m]")

    axs[1, 0].plot(Kas, np.array(paleo_volumes).T[0][0, :], c=cmap(0.25), marker="x")
    axs[1, 0].plot(Kas, np.array(paleo_volumes).T[0][1, :], c=cmap(0.5), marker="x")
    axs[1, 0].plot(Kas, np.array(paleo_volumes).T[0][2, :], c=cmap(0.75), marker="x")
    axs[1, 0].set_xscale("log")
    axs[1, 0].set_xlabel(r"$K_a$ [m/day]")
    axs[1, 0].set_ylabel(r"$V_{paleo}^{pre}$ [$m^3$]")


    axs[1, 1].plot(Kas, np.array(mixing_volumes).T[0, 0, :], c=cmap(0.25), marker="x")
    axs[1, 1].plot(Kas, np.array(mixing_volumes).T[0, 1, :], c=cmap(0.5), marker="x")
    axs[1, 1].plot(Kas, np.array(mixing_volumes).T[0, 2, :], c=cmap(0.75), marker="x")
    axs[1, 1].set_xscale("log")
    axs[1, 1].set_xlabel(r"$K_a$ [m/day]")
    axs[1, 1].set_ylabel(r"$V_{mix}^{pre}$ [$m^3$]")

    line1 = Line2D([0], [0], label="$1e^{-5}$", color=cmap(0.25), marker="x")
    line2 = Line2D([0], [0], label="$1e^{-4}$", color=cmap(0.5), marker="x")
    line3 = Line2D([0], [0], label="$1e^{-3}$", color=cmap(0.75), marker="x")

    style1 = Line2D([0], [0], label="toes", color=cmap(0.5), marker="x")
    style2 = Line2D([0], [0], label="$tips$", color=cmap(0.5), ls="--", marker=".")

    axs[1,1].legend(title=r"$K_b^V$ [m/day]",  handles=[line1, line2, line3], loc='center left', bbox_to_anchor=(1, 0.5))
    axs[0,1].legend(title=r"metric", handles=[style1, style2])
    f.tight_layout()

    plt.show()

def plot_results_5():
    x_step = 1400
    y_step = 10
    ensemble_name="base10"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    heads = [-5, -1, 0]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
            K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
            z0=43.333, L=10000, x_b=3000,  
            h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 

    pars2 = mf6.build_steady_model(pars2, True)
    f,axs = plt.subplots(3, 2, figsize=(8.5, 5), layout="constrained", sharex=True, sharey=True)

    for i, (Ka, KbV, head) in enumerate(zip(Kas, KbVs, heads)):
        name = rf"{Ka}_{KbV}_{head}"
        results = load_smaller(name)
        conc = results["conc"]
        head = results["head"]
        qx = results["qx"]
        qz = results["qz"]
        # pars = load_parameters(f"{name}", backup=True)
        x = np.linspace(-pars2.x_b, pars2.L, pars2.ncol)/1000
        y = np.linspace(-pars2.z0, pars2.Lz-pars2.z0, pars2.nlay)
        plt.rcParams.update({'font.size': 9})
        for k in range(pars2.nlay):
            for j in range(pars2.ncol):
                if abs(conc[-1][k, 0, j]) == 1.e30:
                    conc[-1][k, 0, j] = np.nan
                    head[-1][k, 0, j] = np.nan

        hf_pre = proc.hf(conc[-1][:, 0, :], head[-1][:, 0, :], pars2)
        qx_pre = qx[-1][::-y_step, 0, ::x_step] / np.sqrt(qx[-1][::-y_step, 0, ::x_step]**2+np.abs(qz[-1][::-y_step, 0, ::x_step])**2)
        qz_pre = qz[-1][::-y_step, 0, ::x_step] / np.sqrt(qx[-1][::-y_step, 0, ::x_step]**2+np.abs(qz[-1][::-y_step, 0, ::x_step])**2)
        pcm = axs[i][0].pcolormesh(x, y, conc[-1][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
        axs[i][0].quiver(x[::x_step], y[::y_step], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
        axs[i][0].set_box_aspect(0.3)
        
        head_min = np.nanmin(hf_pre)
        head_max = np.nanmax(hf_pre)
        cp0 = axs[i][1].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
        axs[i][1].set_box_aspect(0.3)
        axs[i][0].set_ylabel("Depth [masl]")
        axs[i][1].clabel(cp0, inline=True, fontsize=7)
        
        i = i+1

    cb = f.colorbar(pcm, ax=axs[:, 0], shrink=0.6, location='left', ticks=[0, 0.35, 3.5, 35])
    cb.set_label("C [PSU]")
    cb.ax.set_yticklabels(["0", "0.35", "3.5", "35"])
    axs[2][0].set_xlabel("Distance offshore [km]")
    axs[2][1].set_xlabel("Distance offshore [km]")
    axs[0][0].title.set_text('Concentration and Flux')
    axs[0][1].title.set_text('Head')
    plt.show()      

def plot_results_6():

    ensemble_name="base10"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    heads = [0, -1, -2, -5]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
            K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
            z0=43.333, L=10000, x_b=3000,  
            h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 
    
    pars2 = mf6.build_steady_model(pars2, True)

    volume_salinized = []
    fraction_salinized = []
    paleo_salinized = []
    mixing_changed = []
    extent_movement = []
    distance_onshore =[]

    for Ka in Kas:
        for KbV in KbVs:
            old_name = f"10km_{Ka}_{KbV}"
            old_results = load_smaller(old_name)

            old_volume = proc.fresh_volume_total([old_results["conc"][1]], pars2)[0]
            new_volume = proc.fresh_volume_total([old_results["conc"][-1]], pars2)[0]

            volume_salinized.append(old_volume-new_volume)
            fraction_salinized.append((old_volume-new_volume)/old_volume)

            old_paleo, x_paleo = proc.paleo_volume([old_results["conc"][1]], [old_results["qx"][1]], pars2)
            old_paleo = old_paleo[0]
            new_paleo = proc.paleo_volume([old_results["conc"][-1]], [old_results["qx"][-1]], pars2, set_x_paleo=x_paleo)[0]
            
            paleo_salinized.append(old_paleo-new_paleo)

            old_mixing = proc.mixing([old_results["conc"][1]], pars2)[0][0]
            new_mixing = proc.mixing([old_results["conc"][-1]], pars2)[0][0]

            mixing_changed.append(new_mixing-old_mixing)
            
            old_extent = proc.extent(old_results["conc"][1], pars2)
            new_extent = proc.extent(old_results["conc"][-1], pars2)
            distance_onshore.append(new_extent)
            extent_movement.append(old_extent-new_extent)

            for head in heads:
                new_name = f"{Ka}_{KbV}_{head}"
                new_results = load_smaller(new_name)

                new_volume = proc.fresh_volume_total([new_results["conc"][-1]], pars2)[0]
                volume_salinized.append(old_volume-new_volume)
                fraction_salinized.append((old_volume-new_volume)/old_volume)

                new_paleo = proc.paleo_volume([new_results["conc"][-1]], [new_results["qx"][-1]], pars2, set_x_paleo=x_paleo)[0] 
                paleo_salinized.append(old_paleo-new_paleo)

                new_mixing = proc.mixing([new_results["conc"][-1]], pars2)[0][0]
                mixing_changed.append(new_mixing-old_mixing)
                
                new_extent = proc.extent(new_results["conc"][-1], pars2)
                distance_onshore.append(new_extent)
                extent_movement.append(old_extent-new_extent)

    f, axs = plt.subplots(3, 2, figsize = (8.45, 12), sharex=True)
    names = [r"$\Delta V_{ofg}^{modern}$ [$m^3$]", r"Fraction OFG salinized [-]", r"$\Delta V_{paleo}^{modern}$ [$m^3$]", 
             r"$\Delta V_{mix}^{modern}$ [$m^3$]", r"$\Delta x_{interface}^{modern}$ [km]", r"$x_{interface}^{post}$ [km]"]
    styles = ["solid", "dotted", "dashed"]
    markers = [".", "+", "x"]
    cmap = cm.get_cmap("Greens")
    colors = [cmap(0.25), cmap(0.5), cmap(0.75)]
    results = [volume_salinized, fraction_salinized, paleo_salinized, mixing_changed, extent_movement, distance_onshore]
    i = 0 
    
    heads = [1, 0, -1, -2, -5]
    scales = ["log", "linear", "log", "log", "log", "symlog"]
    for row in range(3):
        for col in range(2):
            result = np.array(results[i]).reshape((3, 3, 5))
            for ia, Ka in enumerate(Kas):
                for ib, KbV in enumerate(KbVs):
                    axs[row][col].plot(heads, result[ia,ib,:], ls = styles[ia], marker=markers[ia], c=colors[ib], label=rf"$K_a$={Ka}, $K_b^V$={KbV}")
            axs[row][col].set_ylabel(names[i])
            axs[row][col].set_xlim((1,-5))
            axs[row][col].set_yscale(scales[i])
            i = i+1

    axs[-1][0].set_xlabel(r"$h_{modern}$ [m]")
    axs[-1][1].set_xlabel(r"$h_{modern}$ [m]")
    leg = axs[-1][0].legend(loc='upper left', bbox_to_anchor=(0, -0.10), ncol=3)
    leg.set_draggable(state=True)
    plt.show()

    x_step = 1400
    y_step = 10
    ensemble_name="base10"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    heads = [-1, -1, -1]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
            K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
            z0=43.333, L=10000, x_b=3000,  
            h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 

    pars2 = mf6.build_steady_model(pars2, True)
    f,axs = plt.subplots(3, 2, figsize=(8.5, 5), layout="constrained", sharex=True, sharey=True)

    for i, (Ka, KbV, head) in enumerate(zip(Kas, KbVs, heads)):
        name = rf"{Ka}_{KbV}_{head}"
        results = load_smaller(name)
        conc = results["conc"]
        head = results["head"]
        qx = results["qx"]
        qz = results["qz"]
        # pars = load_parameters(f"{name}", backup=True)
        x = np.linspace(-pars2.x_b, pars2.L, pars2.ncol)/1000
        y = np.linspace(-pars2.z0, pars2.Lz-pars2.z0, pars2.nlay)
        plt.rcParams.update({'font.size': 9})
        for k in range(pars2.nlay):
            for j in range(pars2.ncol):
                if abs(conc[-1][k, 0, j]) == 1.e30:
                    conc[-1][k, 0, j] = np.nan
                    head[-1][k, 0, j] = np.nan

        hf_pre = proc.hf(conc[-1][:, 0, :], head[-1][:, 0, :], pars2)
        qx_pre = qx[-1][::-y_step, 0, ::x_step] / np.sqrt(qx[-1][::-y_step, 0, ::x_step]**2+np.abs(qz[-1][::-y_step, 0, ::x_step])**2)
        qz_pre = qz[-1][::-y_step, 0, ::x_step] / np.sqrt(qx[-1][::-y_step, 0, ::x_step]**2+np.abs(qz[-1][::-y_step, 0, ::x_step])**2)
        pcm = axs[i][0].pcolormesh(x, y, conc[-1][::-1, 0, :], cmap="viridis", norm=colors.SymLogNorm(0.4, vmin=0.35, clip=True))
        axs[i][0].quiver(x[::x_step], y[::y_step], qx_pre/2, qz_pre/2, color="white", width=0.002, scale=20)
        axs[i][0].set_box_aspect(0.3)
        
        head_min = np.nanmin(hf_pre)
        head_max = np.nanmax(hf_pre)
        cp0 = axs[i][1].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
        axs[i][1].set_box_aspect(0.3)
        axs[i][0].set_ylabel("Depth [masl]")
        axs[i][1].clabel(cp0, inline=True, fontsize=7)
        
        i = i+1

    cb = f.colorbar(pcm, ax=axs[:, 0], shrink=0.6, location='left', ticks=[0, 0.35, 3.5, 35])
    cb.set_label("C [PSU]")
    cb.ax.set_yticklabels(["0", "0.35", "3.5", "35"])
    axs[2][0].set_xlabel("Distance offshore [km]")
    axs[2][1].set_xlabel("Distance offshore [km]")
    axs[0][0].title.set_text('Concentration and Flux')
    axs[0][1].title.set_text('Head')
    plt.show()      

def plot_results_7():

    ensemble_name="base10"            
    Kas = [1, 10, 100]
    KbVs = [1e-5, 1e-4, 1e-3]
    heads = [0, -1, -2, -5]
    anis = 100
    n = 0.5
    h_modern = 1
    inputs = []
    epsilons = []
    raleighs = []
    names = []
    
    # pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
    #         K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
    #         z0=43.333, L=10000, x_b=3000,  
    #         h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 
    
    # pars2 = mf6.build_steady_model(pars2, True)
    # Ths = []
    # Tvs = []
    # Ts = []
    # overs = []

    # for Ka in Kas:
    #     for KbV in KbVs:
    #         name = rf"10km_{Ka}_{KbV}"
    #         results = load_smaller(name)
    #         _, extent0 = proc.extent(results["conc"][1], pars2, True)
    #         zextent0 = proc.get_z_extent_cell(results["conc"][1], pars2)

    #         for head in heads:

    #             name = rf"{Ka}_{KbV}_{head}"
    #             results = load_smaller(name)
    #             conc = results["conc"]
    #             head = results["head"]
    #             qx = results["qx"]
    #             z = results["qz"]
                
    #             dl0, dhh, dd0, dhv, T = proc.time_constants(results["conc"], results["head"][0], pars2, extent0, zextent0)
    #             T_h = np.abs(dl0**2*n/(dhh*Ka*365))
    #             T_v = np.abs(dd0**2*n/(dhv*KbV*365))
    #             residual = np.min([T_h, T_v]) - T
    #             over = residual > 0

    #             Ths.append(T_h)
    #             Tvs.append(T_v)
    #             Ts.append(T)
    #             overs.append(over)
    

    results = np.load("results.npy")
    overs = results[0].astype(bool)
    Ths = results[1]
    Tvs = results[2]
    Ts = results[3]


    f, ax = plt.subplots(figsize=(5,4), layout="tight")
    ax.scatter(np.minimum(Ths, Tvs), Ts)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.rcParams.update({'font.size': 9})
    ax.hlines(100, 1, 10000, colors="grey", linestyles="dashed")
    ax.plot(np.linspace(10, 100), np.linspace(10, 100), color="grey", ls="--")
    ax.set_ylim(10, 1000)
    ax.set_aspect('equal')
    ax.set_ylabel(r"$T_{sal}$ [years]")
    ax.set_xlabel(r"$T_{h}$ [years]")
    # vmin = np.min(np.log10(Ts))
    # vmax = np.max(np.log10(Ts))
    f.savefig('figure7.png', dpi=600)
    pass
    # Ths = np.array(Ths)
    # Tvs = np.array(Tvs)
    # Ts = np.array(Ts)

    
    # norm = colors.SymLogNorm(0.5, vmin=min(np.log10(outputs["T_sal"][active[0]])), vmin=min(np.log10(outputs["T_sal"][active[0]]))vmax=1e4,)
    # colors_plot = cm.coolwarm_r(norm(outputs["T_sal"][active[0] & (inputs[:,2]<0.5)]))
    # norm=colors.TwoSlopeNorm((vmin+vmax)/2, vmin=vmin, vmax=vmax)
    sc = ax.scatter(Ths[overs] , Tvs[overs], c=np.log10(Ts[overs]), cmap="coolwarm_r", norm=colors.TwoSlopeNorm(2, vmin=vmin, vmax=vmax), marker="^", s=40, edgecolors="w")
    ax.scatter(Ths[~overs] , Tvs[~overs], c=np.log10(Ts[~overs]), cmap="coolwarm_r", norm=colors.TwoSlopeNorm(2, vmin=vmin, vmax=vmax), marker="v", s=40, edgecolors="w")

    # with anisotropy ranking
    # ax.scatter(T_h[active[0] & (inputs[:,2]<0.5)], T_v[active[0] & (inputs[:,2]<0.5)], c="w", cmap="coolwarm_r", norm=colors.SymLogNorm(0.5, vmin=1e2, vmax=1e4,), marker="^", edgecolors=colors_plot)
    # ax.scatter(T_h[active[0] & (inputs[:,2]>0.5) & (inputs[:,2]<5)] , T_v[active[0] & (inputs[:,2]>0.5) & (inputs[:,2]<5)], c=outputs["T_sal"][active[0] & (inputs[:,2]>0.5) & (inputs[:,2]<5)], cmap="coolwarm_r", norm=colors.SymLogNorm(0.5, vmin=1e2, vmax=1e4,), marker="^")
    # ax.scatter(T_h[active[0] & (inputs[:,2]>5)], T_v[active[0] & (inputs[:,2]>5)], c=outputs["T_sal"][active[0] & (inputs[:,2]>5)], cmap="coolwarm_r", norm=colors.SymLogNorm(0.5, vmin=1e2, vmax=1e4,), marker="^", edgecolors="k")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"$T_h$ [years]")
    ax.set_ylabel(r"$T_v$ [years]")
    cb = plt.colorbar(sc, ax=ax, cmap="coolwarm_r", location="right", pad=0.1, shrink=0.75)
    cb.ax.set_title(r"$T_{sal}$ [years]")
    cb.set_ticks([np.log10(20), np.log10(50), np.log10(100), np.log10(200), np.log10(500)])
    cb.set_ticklabels(["20", "50","100","200", "500"])
    ax.set_box_aspect(1)
    #ax.set_xlim(right=1e6)
    # ax.set_ylim(top=1e5)
    legend_elements = [Line2D([0], [0], marker='^', color="w", markerfacecolor="black", label="overestimate", markeredgecolor="black", markeredgewidth=0.5),
                       Line2D([0], [0], marker='v', color="w", markerfacecolor="black", label="underestimate", markeredgecolor="black", markeredgewidth=0.5)]
    ax.legend(handles=legend_elements, ncol=2)
    ax.plot(np.linspace(1e2,1e4), np.linspace(1e2,1e4), lw=1, ls="--", c="k")
    plt.show()


if __name__ == "__main__":
    # plot_predev_bases()
    # plot_results_1()
   # plot_results_2()
    # plot_results_3()
    # plot_results_4()
    # plot_results_5()
    # plot_results_6()
    plot_results_7()