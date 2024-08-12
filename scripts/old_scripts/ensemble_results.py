from multiprocessing import Pool
import results_mf6 as mf6r
from pars import ModelParameters, load_parameters
import numpy as np
import matplotlib.pyplot as plt
import os
from mf6 import get_results, get_last_results
import pandas as pd
import pdb

def plot_metrics(name, times, axs2, axs3, axs4, axs5, axs3twin, axs4twin, axs5twin, y_labels="none", legend=False):

    df = pd.read_csv(f"/home/superuser/results/{name}/metrics.csv")
    
    # plot fresh
    axs2.plot(times, df.fresh_list)
    axs3twin.plot(times, df.toe_list)
    axs3twin.plot(times, df.tip_list)
    axs4.plot(times, df.volume_list)
    axs4twin.plot(times, df.centroid_list, linestyle="--")
    axs5.plot(times, df.fresh_flux_list)
    axs5.plot(times, df.sal_flux_list)
    axs5.set_xlabel("time from present (days)")
    axs5twin.plot(times, df.fresh_centroid_list, linestyle="--")
    axs5twin.plot(times, df.sal_centroid_list, linestyle="--")


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
    print(times)
    i_predev = list(times).index(max(times[times<=0]))
    # i_predev=0
    print(times[i_predev])
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    conc_predev = np.flipud(conc[i_predev][:, 0, :])
    qx_predev = np.flipud(qx[i_predev][:, 0, :])
    qz_predev = np.flipud(qz[i_predev][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_predev[i, j]) == np.float64(1.e30):
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
            if abs(conc_post[i, j]) == np.float64(1.e30):
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


    conc, _, qx, qz, times = get_results("d0")
    mf6r.metrics("d0", conc, qz)
    axs[0][mid], axs[1][mid] = plot_results("d0", conc, qx, qz, times, axs[0][mid], axs[1][mid])
    axs[2][mid], axs[3][mid], axs[4][mid], axs[5][mid], d[f"axs3_{mid}twin"], d[f"axs5_{mid}twin"],d[f"axs5_{mid}twin"] \
        = plot_metrics("d0", times, axs[2][mid], axs[3][mid], axs[4][mid], axs[5][mid], d[f"axs3_{mid}twin"], d[f"axs4_{mid}twin"],d[f"axs5_{mid}twin"])

    i = 0
    y_labels= "left"
    for r in set:
        if i >= mid:
           i=i+1
        conc, _, qx, qz, times = get_results(r)
        mf6r.metrics(r, conc, qz)
        axs[0][i], axs[1][i] = plot_results(r, conc, qx, qz, times, axs[0][i], axs[1][i], y_labels=y_labels)
        axs[2][i], axs[3][i], axs[4][i], axs[5][i], d[f"axs3_{i}twin"], d[f"axs4_{i}twin"],d[f"axs5_{i}twin"]  \
            = plot_metrics(r, times, axs[2][i], axs[3][i], axs[4][i], axs[5][i], d[f"axs3_{i}twin"], d[f"axs4_{i}twin"],d[f"axs5_{i}twin"], y_labels=y_labels)
        i=i+1
        if i < len(set)-2:
            y_labels="none"
        else:
            y_labels="right"
            
    axs[5][0].set_yscale("log")
    axs[5][1].set_yscale("log")
    axs[5][2].set_yscale("log")
    
    f.tight_layout()

    if not os.path.exists(f"/home/superuser/results/initial_cases_d"):
        os.makedirs(f"/home/superuser/results/initial_cases_d")
    f.savefig(f"/home/superuser/results/initial_cases_d/{set_name}", dpi=300) 

def plot_sgd_results(name, times, conc, qx, qz, axs0, axs1, axs2, axs3, y_labels=None):
    
    pars = load_parameters(name)
    i_predev = 0 # list(times).index(max(times[times<=0]))
    
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    
    conc_predev = np.flipud(conc[i_predev][:, 0, :])
    qx_predev = np.flipud(qx[i_predev][:, 0, :])
    qz_predev = np.flipud(qz[i_predev][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_predev[i, j]) == np.float64(1.e30):
                conc_predev[i,j] = np.nan

    sgds = [qz[i_predev][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    
    c_sgds = [conc[i_predev][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    sgds_fresh = np.ma.masked_where(np.array(c_sgds) < 0.35, np.array(sgds))
    sgds_sal = np.ma.masked_where(np.array(c_sgds) >= 0.35, np.array(sgds))
    axs0.plot(x, sgds_fresh, c=[253/255,231/255,37/255])
    axs0.plot(x, sgds_sal, c=[68/255,1/255,84/255])
    # axs0.set_yscale("log")
    sgds_max = np.max(sgds)
	
    conccm = axs1.pcolormesh(x, y, conc_predev, 
             	                cmap="viridis", vmax=35.7, vmin=0)

    qx_predev = qx_predev / np.sqrt(qx_predev**2+qz_predev**2)
    qz_predev = qz_predev / np.sqrt(qx_predev**2+qz_predev**2)
    axs1.quiver(x[::80], y[::50], qx_predev[::50, ::80],-1*qz_predev[::50, ::80], 
                     color="white", width=0.002)
    axs1.set_aspect(50)


    conc_post = np.flipud(conc[-1][:, 0, :])
    qx_post = np.flipud(qx[-1][:, 0, :])
    qz_post = np.flipud(qz[-1][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_post[i, j]) == np.float64(1.e30):
                conc_post[i,j] = np.nan
                
    sgds = [qz[-1][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    
    c_sgds = [conc[-1][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    sgds_fresh = np.ma.masked_where(np.array(c_sgds) < 0.35, np.array(sgds))
    sgds_sal = np.ma.masked_where(np.array(c_sgds) >= 0.35, np.array(sgds))
    
    axs2.plot(x, sgds_fresh, c=[253/255,231/255,37/255])
    axs2.plot(x, sgds_sal, c=[68/255,1/255,84/255])
    sgds_min = np.min(sgds)
    
    axs0.set_ylim([-2e-5, 3e-5])
    axs2.set_ylim([-2e-5, 3e-5])

    conccm = axs3.pcolormesh(x, y, conc_post, 
             	                cmap="viridis", vmax=35.7, vmin=0)
    #conccb = plt.colorbar(conccm, shrink=1)
    # conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    qx_post = qx_post / np.sqrt(qx_post**2+qz_post**2)
    qz_post = qz_post / np.sqrt(qx_post**2+qz_post**2)
    axs3.quiver(x[::80], y[::50], qx_post[::50, ::80],-1*qz_post[::50, ::80], 
                     color="white", width=0.002)
    axs3.set_aspect(50)
    #ax.axhline(pars.sea_level

    if y_labels=="left":
        axs1.set_ylabel("Distance above present sea level (m)")
        axs3.set_ylabel("Distance above present sea level (m)")
        
    axs3.set_xlabel("Distance offshore (m)")

    return axs0, axs1, axs2, axs3
    
def plot_sgd_metrics(name, times, axs4, axs4twin, y_labels="none", legend=False):

    df = pd.read_csv(f"/home/superuser/results/{name}/metrics.csv")
    
    # plot fresh
    axs4.plot(times, df.fresh_flux_list)
    axs4.plot(times, df.sal_flux_list)
    axs4.set_xlabel("time from present (days)")
    axs4twin.plot(times, df.fresh_centroid_list, linestyle="--")
    axs4twin.plot(times, df.sal_centroid_list, linestyle="--")
    axs4.set_yscale("log")

    if y_labels == "left":
        axs4.legend(["fresh sgd", "saline sgd"])
        
    elif y_labels == "right":
        axs4twin.set_ylabel("distance offshore (m^3)")
        axs4twin.legend(["fresh centroid", "saline centroid"])
        
    return axs4, axs4twin
    
def plot_pair_sgd(set, set_name):

    f, axs = plt.subplots(5, 2, sharey="row", figsize=(15,10))
 
    f.suptitle(set_name)
    
    d = {}

    for i in range(len(set)):
        d[f"axs4_{i}twin"] = axs[4][i].twinx()
    
    # for i in range(len(set)):
    #     d["axs3_0twin"].get_shared_y_axes().join(d["axs3_0twin"], d[f"axs3_{i+1}twin"])


    for i, (r, label)  in enumerate(zip(set, ["left", "right"])):
        conc, _, _, qz, times = get_results(r)
        mf6r.metrics(r, conc, qz)
        pdb.set_trace()
        axs[0][i], axs[1][i], axs[2][i], axs[3][i] = plot_sgd_results(r, times, conc, qx, qz, axs[0][i], axs[1][i], axs[2][i], axs[3][i], y_labels=label)
        # axs[4][i],d[f"axs4_{i}twin"] = plot_sgd_metrics(r, times, axs[4][i], d[f"axs4_{i}twin"], y_labels=label)
    
    f.tight_layout()

    if not os.path.exists(f"/home/superuser/results/sgds"):
        os.makedirs(f"/home/superuser/results/sgds")
    f.savefig(f"/home/superuser/results/sgds/{set_name}", dpi=300) 
    
def plot_single_sgd(name):

    f, axs = plt.subplots(3, 2, figsize=(12,6))
    conc, qz, times = get_results(name)
    _, _, qx_last, _, _ = get_last_results(name)
    mf6r.metrics(name, conc, qz)
    pars = load_parameters(name)
    i_predev = list(times).index(max(times[times<=0]))-1
    times=times/365
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    
    conc_predev = np.flipud(conc[i_predev][:, 0, :])
    qz_predev = np.flipud(qz[i_predev][:, 0, :])
    
    qx_predev = np.flipud(qx_last[0][:, 0, :])
    qx_post = np.flipud(qx_last[-1][:, 0, :])

    
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_predev[i, j]) == np.float64(1.e30):
                conc_predev[i,j] = np.nan
               
    
    axs[1,0].set_xlim([-pars.x_b, pars.L])
    axs[1,1].set_xlim([-pars.x_b, pars.L])

    axs[2,0].set_ylim([1e-4, 1e-1])
    axs[2,1].set_ylim([1e-4, 1e-1]) 

    sgds = [qz[i_predev][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    
    c_sgds = [conc[i_predev][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    sgds_fresh = np.ma.masked_where(np.array(c_sgds) < 0.35, np.array(sgds))
    sgds_sal = np.ma.masked_where(np.array(c_sgds) >= 0.35, np.array(sgds))
    axs[1, 0].plot(x, sgds_fresh/10, c=[253/255,231/255,37/255])
    axs[1, 0].plot(x, sgds_sal/10, c=[68/255,1/255,84/255])
    axs[1,0].hlines(y=0, xmin=-pars.x_b, xmax=pars.L, colors='lightgrey', ls='dotted', zorder=-1)
    # axs0.set_yscale("log")
    sgds_max = np.max(sgds)
	
    conccm = axs[0,0].pcolormesh(x, y, conc_predev, cmap="viridis", vmax=35.7, vmin=0)
    
    qx_predev = qx_predev / np.sqrt(qx_predev**2+qz_predev**2)
    qz_predev = qz_predev / np.sqrt(qx_predev**2+qz_predev**2)
    axs[0,0].quiver(x[::200], y[::15], qx_predev[::15, ::200],-1*qz_predev[::15, ::200], 
                     color="white", width=0.002)
    
    conc_post = np.flipud(conc[-1][:, 0, :])
    qz_post = np.flipud(qz[-1][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if abs(conc_post[i, j]) == np.float64(1.e30):
                conc_post[i,j] = np.nan
                
    sgds = [qz[-1][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    
    c_sgds = [conc[-1][cell[0], cell[1], cell[2]] for cell in pars.top_cells]
    sgds_fresh = np.ma.masked_where(np.array(c_sgds) < 0.35, np.array(sgds))
    sgds_sal = np.ma.masked_where(np.array(c_sgds) >= 0.35, np.array(sgds))
    
    axs[1,1].plot(x, sgds_fresh/10, c=[253/255,231/255,37/255])
    axs[1,1].plot(x, sgds_sal/10, c=[68/255,1/255,84/255])
    axs[1,1].hlines(y=0, xmin=-pars.x_b, xmax=pars.L, colors='lightgrey', ls='dotted', zorder=-1)
    sgds_min = np.min(sgds)
    
    axs[1, 0].set_ylim([-3e-6, 2e-6])
    axs[1, 1].set_ylim([-3e-6, 2e-6])
    
    axs[0,0].set_ylabel("masl")
    axs[1,0].set_ylabel(r"flux [$m/day$]")
    axs[2,0].set_ylabel(r"discharge [$m^3/day$]")
    axs[2,0].set_xlabel("years before 'present'")
    axs[2,1].set_xlabel("years after 'present'")
    axs[1,0].set_xlabel("distance offshore [m]")
    axs[1,1].set_xlabel("distance offshore [m]")
    axs[0,0].set_xlabel("distance offshore [m]")
    axs[0,1].set_xlabel("distance offshore [m]")
    axs[0,0].set_title("Pre-development")
    axs[0,1].set_title("Post-development")
    conccm = axs[0,1].pcolormesh(x, y, conc_post, cmap="viridis", vmax=35.7, vmin=0)
    
    qx_post = qx_post / np.sqrt(qx_post**2+qz_post**2)
    qz_post = qz_post / np.sqrt(qx_post**2+qz_post**2)
    axs[0,1].quiver(x[::200], y[::15], qx_post[::15, ::200],-1*qz_post[::15, ::200], 
                     color="white", width=0.002)
    
    axs[0,0].annotate("Concentration", xy=(0, 0.5), xytext=(-axs[0,0].yaxis.labelpad - 5, 0),
                xycoords=axs[0,0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')
    axs[1,0].annotate("Submarine flux", xy=(0, 0.5), xytext=(-axs[1,0].yaxis.labelpad - 5, 0),
                xycoords=axs[1,0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')
    axs[2,0].annotate("Total submarine \n fresh discharge", xy=(0, 0.5), xytext=(-axs[2,0].yaxis.labelpad - 5, 0),
                xycoords=axs[2,0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')
                
    df = pd.read_csv(f"/home/superuser/results/{name}/metrics.csv")
    axs[2, 0].plot(times[:i_predev], df.fresh_flux_list[:i_predev], c=[68/255,1/255,84/255])
    # axs[2, 0].plot(times[:i_predev], df.sal_flux_list[:i_predev], c=[253/255,231/255,37/255])
    axs[2, 1].plot(times[i_predev:], df.fresh_flux_list[i_predev:], c=[68/255,1/255,84/255])
    # axs[2, 1].plot(times[i_predev:], df.sal_flux_list[i_predev:], c=[253/255,231/255,37/255])
            
    axs[2, 0].spines['right'].set_visible(False)
    axs[2, 1].spines['left'].set_visible(False)
    axs[2, 0].yaxis.tick_left()
    axs[2, 0].tick_params(labelright=False)
    axs[2, 1].yaxis.tick_right()
    
    axs[2, 0].set_yscale("log")
    axs[2, 1].set_yscale("log")
    
    # d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    # kwargs = dict(transform=axs[2, 0].transAxes, color='k', clip_on=False)
    # axs[2, 0].plot((1-d, 1+d), (-d, +d), **kwargs)
    # axs[2, 0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    # kwargs.update(transform=axs[2, 1].transAxes)  # switch to the bottom axes
    # axs[2, 1].plot((-d, +d), (1-d, 1+d), **kwargs)
    # axs[2, 1].plot((-d, +d), (-d, +d), **kwargs)
    
    for row in range(3):
    	for col in range(2):
    		axs[row, col].set_box_aspect(1/4)
    
    f.tight_layout()
    if not os.path.exists(f"/home/superuser/sloping_confined_aquifer/results/sgds"):
        os.makedirs(f"/home/superuser/sloping_confined_aquifer/results/sgds")
    f.savefig(f"/home/superuser/sloping_confined_aquifer/results/sgds/{name}", dpi=300)     	                
    plt.show()

def plot_ensemble(sets, set_names, name):

    # get everything for base case
    # plot_set([sets[0], set_names[0]])
    inputs = [[s, n] for s, n in zip(sets, set_names)]
    p = Pool(processes=8)
    p.map(plot_set, inputs)

def main():
    sets = [[f"c{i}", f"c{i+1}"] for i in range(1, 13, 2)]
    sets = [["d1", "d2"], ["d9", "d10"]]
    print(sets)
    
    set_names = ["K_aquifer", "K_aquitard", "anis_aquifer", "anis_aquitard", "h_modern", "aquitard_D"]
    set_names = ["K_aquifer"]
    
    # plot_ensemble(sets, set_names, "set_d_model_tweaks")
    plot_single_sgd("l0")
     
if __name__=="__main__":
        main()
