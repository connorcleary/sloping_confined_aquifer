from pars import ModelParameters, load_parameters
import numpy as np
import sloping_confined_aquifer as sca
import matplotlib.pyplot as plt
import plot_helpers as plth
import os
from matplotlib import animation
#import post_processing as proc

def plot_results(name, timestep=-1, row=0, return_axs=False, figsize=(12,6),
    cmap="viridis", arrow_c="white", aspect=100, vector_T=5, width=0.002, fmt="%3.2f"):
    """
        Plot the head, concentrations and fluxes

        Inputs:
            name: model to plot
            timestep: timestep to plot, default is -1 (the last one)
            row: row to plot
            return_axs: flag whether to return the axes objects
            figsize: figure dimensions (inches)
            cmap: colormap
            arrow_c: color of arrows
            aspect: vertical exageration 
            vector_T: spaces between vector arrows
            width: arrow width
            fmt: format of contour labels
        Outputs:
            axs: axes objects (optional)
    """

    # load parameters and results
    pars = load_parameters(name)
    concentration, head, qx, qy, qz = sca.load_results(name)

    f, axs = plt.subplots(2, 1, figsize=figsize)

    # set up x and y arrays, to be distance above sea level and distance onshore
    x = np.linspace(-pars.Lx*pars.offshore_proportion, pars.Lx-pars.Lx*pars.offshore_proportion, pars.ncol)
    y = np.linspace(-pars.sea_level, pars.Lz-pars.sea_level, pars.nlay)

    # select relevent slice in time and the alongshore direction, and set values above the water table as nan
    concentration_array = concentration[timestep,:,row,:]
    head_array = head[timestep,:,row,:] - pars.sea_level*np.ones_like(head[timestep,:,row,:])
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if concentration_array[i, j] == np.float32(1.e30):
                head_array[i,j] = np.nan
                concentration_array[i, j] = np.nan

    # plot head colormesh
    headcm = axs[0].pcolormesh(x, y, np.flipud(head_array), 
            	                cmap=cmap, vmax=np.nanmax(head_array[:, 1:]), vmin=np.min([-1, pars.h_b]))

    # plot head contours
    # hc = axs[0].contour(x, y, np.flipud(head_array), colors=arrow_c, 
    #                     levels=np.linspace(np.min([0, pars.h_b]), np.nanmax(head_array[:, 1:]), 15))

    # label contours
    # axs[0].clabel(hc, hc.levels, inline=True, fontsize=10, fmt=fmt)

    # plot concentration colormesh
    conccm = axs[1].pcolormesh(x, y, np.flipud(concentration_array), 
            	                cmap=cmap, vmax=35, vmin=0)

    # plot arrows
    # axs[1].quiver(plth.sample_grid_data(x, vector_T, vector_T), plth.sample_grid_data(y, vector_T, vector_T), 
    #                 plth.sample_grid_data(np.flipud(qx[timestep,:,row,:]), vector_T, vector_T), 
    #                 -plth.sample_grid_data(np.flipud(qz[timestep,:,row,:]), vector_T, vector_T), 
    #                 color=arrow_c, width=width)
    axs[0].set_aspect(aspect)
    axs[1].set_aspect(aspect)
    axs[0].set_title("Head")
    axs[1].set_title("Salinity")
    axs[0].set_ylabel("Height above sealevel (m)")
    axs[1].set_ylabel("Height above sealevel (m)")
    axs[1].set_xlabel("Distance onshore (m)")

    f.suptitle(f"Head and salinity distributions for {name}")

    headcb = plt.colorbar(headcm, shrink=1, ax=axs[0])
    conccb = plt.colorbar(conccm, shrink=1, ax=axs[1])
    headcb.ax.set_title('Head (m)', fontsize = 'small')
    conccb.ax.set_title('Salinity (kg/m^3)', fontsize = 'small')
    
    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)
    plt.savefig(f"{ws}\\head_and_concentration", dpi=300)

    # return axs objects if necessary
    if return_axs: return axs

def animate_func(num, ax, pars, concentration, qx, qz):
    ax.clear()
    x = np.linspace(-pars.Lx*pars.offshore_proportion, pars.Lx-pars.Lx*pars.offshore_proportion, pars.ncol)
    y = np.linspace(-pars.sea_level, pars.Lz-pars.sea_level, pars.nlay)
    conc = np.flipud(concentration[num][:, 0, :])
    qx = np.flipud(qx[num][:, 0, :])
    qz = np.flipud(qz[num][:, 0, :])

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if conc[i, j] == np.float32(1.e30):
                conc[i,j] = np.nan

    # conccm = ax.pcolormesh(x, y, conc, 
    #         	                cmap="viridis", vmax=35, vmin=0)
    conccon = ax.contourf(x, y, conc, cmap="binary", vmax=36, vmin=0, levels=[0, 0.5,  5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 35.75])
    # plot arrows
    # ax.quiver(x[::20], y[::5], 10*qx[::5, ::20],-10*qz[::5, ::20], 
    #                  color="white", width=0.002)
    ax.set_aspect(10)
    return ax
    

def plot_gif(name):
    pars = load_parameters(name)
    concentration,_,qx,_,qz = sca.load_results(name)
    qx_plot = []
    qz_plot = []
    concentration_plot = []
    N= 0
    for i, (c, x, z) in enumerate(zip(concentration, qx, qz)):
        if i % 60 == 0:
            concentration_plot.append(c)
            qx_plot.append(x)
            qz_plot.append(z)
            N += 1
    

    fig = plt.figure(figsize =(14, 9))
    ax = plt.axes()
    ani = animation.FuncAnimation(fig, animate_func, fargs = (ax, pars, concentration_plot, qx_plot, qz_plot), interval=100,   
                                   frames=N)
    f = f"animate_func_{name}.gif"
    writergif = animation.PillowWriter(fps=2*5)#N/6)
    ani.save(f, writer=writergif)
    plt.show()


    


