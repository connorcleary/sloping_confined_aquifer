from scipy.stats import qmc
import pdb
import numpy as np
import os
import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool
import post_processing as proc
from mf6 import get_results, get_results_for_hc
from results_mf6 import save_results_for_saltelli, plot_color_mesh_saltelli, extract_pre_post_volume, paleo_volumes, toe_position
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.tri as tri
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from SALib.sample import sobol as sbsample
from SALib.analyze import sobol
import pairs_plot
import seaborn as sns
from sklearn.preprocessing import PolynomialFeatures
from matplotlib.lines import Line2D
from sklearn.linear_model import LinearRegression
import matplotlib
import matplotlib.colors as colors
from matplotlib import cm
import results_mf6

def fit_polynomial_regression(x1, x2, y, scales, lims):

    poly = PolynomialFeatures(degree=2)
    # poly = PolynomialFeatures(interaction_only=True,include_bias = False)
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
    
    print(clf.score(x_, y))
    return X1, X2, Y.reshape(50, -1)
def scatter_and_heatmap(S, pars, signs):
    f, axs = plt.subplots(1, 2)
    cmap = plt.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, len(pars)))
    markers = []
    for sign in signs:
        if sign == np.sign(1.0):
            markers.append("^")
        else: 
            markers.append("v")

    for i in range(len(colors)):
        axs[0].scatter(S["S1"][i], S["ST"][i], c=colors[i], label=pars[i], marker=markers[i], edgecolors="black", linewidths=0.5)
    axs[0].plot(np.linspace(0, 1 ,5), np.linspace(0, 1, 5), lw=0.5, ls="--", color="grey")
    axs[0].set_box_aspect(1)  
    axs[0].legend()
    axs[0].set_xlabel("S")
    axs[0].set_ylabel("ST")
    axs[0].set_ylim([0, 1])
    axs[0].set_xlim([0, 1])
    i_lower = np.tril_indices(len(pars), -1)
    S["S2"][i_lower] = S["S2"].T[i_lower]
    sns.heatmap(S["S2"],  vmin=0, cmap="pink_r", robust=False, annot=True, square=True, yticklabels=pars,xticklabels=pars, ax=axs[1], cbar=False, fmt=".2f")
    plt.show()

def fit_and_plot_regression(name, pars, inputs, outputs):
    inputs = np.asarray(inputs)
    inputs[:,-2] = -inputs[:,-2] + 2
    X_train, X_test, y_train, y_test = train_test_split(inputs, outputs, test_size=0.25)
    reg = LinearRegression().fit(np.log10(inputs), np.log10(outputs))
    R2 = reg.score(np.log10(inputs), np.log10(outputs))
    print(R2)
    # f, ax = plt.subplots(figsize=(4,4))
    # plt.rcParams.update({'font.size': 9})
    # ax.scatter(np.log10(y_test), reg.predict(np.log10(X_test)))
    # ax.plot(np.linspace(-1.2, 0), np.linspace(-1.2,0), color="grey", ls="--")
    # ax.set_xlabel(r"$log_{10} \frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$  (modelled)")
    # ax.set_ylabel(r"$log_{10} \frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$  (predicted)")
    # plt.tight_layout()
    # plt.show()
    outputs = results_mf6.get_outputs_for_timescales(name, 1024, inputs)["T_sal"]
    
    f, axs = plt.subplots(1, 2, figsize=(8,4), sharey=True, layout="tight")
    plt.rcParams.update({'font.size': 9})
    reg_full = LinearRegression().fit(np.log10(inputs), np.log10(outputs))
    predicted = 10**reg_full.predict(np.log10(inputs))
    residuals = outputs-predicted
    h_shifted = inputs[:,-2]
    h_driving = -(-(inputs[:,-2]-2)-(0.025/1)*35)
    # ax.scatter(h_driving, predicted, c=residuals, vmin=-1, vmax=1, cmap="coolwarm", s=5)
    h_range = np.linspace(min(h_shifted), max(h_shifted))
    h_driving_range = np.linspace(min(h_driving), max(h_driving))
    kbv = 10**((np.log10(0.000001)+np.log10(0.001))/2)
    alpha = 1
    anisa = 10
    anisb = 10
    n = 0.5
    ka = 10**((np.log10(1)+np.log10(100))/2)
    kbvs = [0.001, 0.0001, 0.00001, 0.000001]
    kas = [1, 10, 100]
    cmap = matplotlib.cm.get_cmap('Greens')
    colors = [cmap(i) for i in np.linspace(0, 1, len(kas)+1)]
    for i in range(len(kas)):
        axs[0].plot(h_driving_range, [10**reg_full.predict(np.log10(np.array([kas[i], kbv, alpha, anisa, anisb, h_range[j], n]).reshape(1, -1))) for j in range(len(h_range))], c=colors[i+1], label=str(kas[i]))
    axs[0].set_ylim(0, 1)
    axs[0].set_xlabel(r"$h_f-h_{post}$")
    axs[0].set_ylabel(r"Time to salinization [years]")

    cmap = matplotlib.cm.get_cmap('Oranges')
    colors = [cmap(i) for i in np.linspace(0, 1, len(kbvs)+1)]
    for i in range(len(kbvs)):
        axs[1].plot(h_driving_range, [10**reg_full.predict(np.log10(np.array([ka, kbvs[i], alpha, anisa, anisb, h_range[j], n]).reshape(1, -1))) for j in range(len(h_range))], c=colors[i+1], label=str(kbvs[i]))
    axs[1].set_ylim(0, 1000)
    axs[1].set_xlabel(r"$h_f-h_{post}$")
    #axs[1].set_ylabel(r"Fraction OFG Salinized")
    axs[0].legend(title=r"$K_a$ [m/day]")
    axs[1].legend(title=r"$K_b$ [m/day]")
    axs[0].set_box_aspect(1)
    axs[1].set_box_aspect(1)
    plt.show()
    


def plot_rate_of_intrusion(name, n):
    metrics = []
    rates = []
    for i in range(n):
        metrics.append(np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_metrics.csv", delimiter=","))

    for i in range(n):
        for j in range(2, len(metrics[i][2])):
            if (metrics[i][2][1]-metrics[i][2][j])/metrics[i][2][1] == 1:
                break
        average = (metrics[i][2][1]-metrics[i][2][-1])/(j-1)
        rate = []
        for j in range(2, len(metrics[i][2])):
            rate.append((metrics[i][2][j-1]-metrics[i][2][j])/average)
        rates.append(rate)
    
    return rates

def plot_original_volume_vs_volume_change(original, change, inputs_retrieved, pars,scales):
    # plot the original ofg volume against how much is exhausted, for each parameter.
    # identify which parameters have potential consequences if over estimated.
    inputs_retrieved = np.delete(inputs_retrieved, 4, 1)
    pars = np.delete(pars, 4)
    scales = np.delete(scales, 4)
    f, axs, = plt.subplots(2, 3)
    for i in range(len(inputs_retrieved[0])):
        axs[i//3][i%3].set_yscale("log")
        axs[i//3][i%3].set_xscale("log")
        if scales[i]=="log":
            c = np.log10(inputs_retrieved[:, i])
        else:
            c = inputs_retrieved[:, i]
        sc = axs[i//3][i%3].scatter(original, change, c=c)
        cb = plt.colorbar(sc, ax=axs[i//3][i%3], cmap="viridis", extend="both", location="top")
        cb.ax.set_title(pars[i])
        if i//3 == 1:
            axs[i//3][i%3].set_xlabel("${V^{pre}_{ofg}}$")
        if i%3 == 0:
            axs[i//3][i%3].set_ylabel(r"$\frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
    
    # axs[1, 3].axis('off')
    plt.tight_layout()
    plt.show()

def plot_original_volume_and_total_change(original, change, inputs_retrieved, pars, scales):
    f, axs= plt.subplots(2, 7, figsize=(8, 2))
    for i in range(len(inputs_retrieved[0])):
        axs[0][i].set_yscale("log")
        axs[1][i].set_yscale("log")
        if scales[i]=="log":
            axs[0][i].set_xscale("log")
            axs[1][i].set_xscale("log")

        axs[0][i].scatter(inputs_retrieved[:, i], original, c="grey", s=2)
        axs[1][i].scatter(inputs_retrieved[:, i], change, c="grey", s=2)
        axs[1][i].set_xlabel(pars[i])
        axs[1][i].set_box_aspect(1)
        axs[0][i].set_box_aspect(1)
        axs[0][i].set_xticklabels([])
        if i > 0:
            axs[0][i].set_yticklabels([])
            axs[1][i].set_yticklabels([])


    axs[0][0].set_ylabel(r"$V^{pre}_{ofg}$")
    axs[1][0].set_ylabel(r"$\Delta V^{100}_{ofg}$")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def plot_paleo_volume_and_changes(paleo_volume, inputs_retrieved, pars, scales):
    f, axs= plt.subplots(3, 7, figsize=(8, 4))
    paleo_volume = np.asarray(paleo_volume)
    paleo_volume_fraction = []
    for i in range(len(inputs_retrieved)):
        if paleo_volume[i, 0] == 0:
            paleo_volume_fraction.append(0)
        else: 
            paleo_volume_fraction.append((paleo_volume[i, 0]-paleo_volume[i, 1])/paleo_volume[i, 0])
    for i in range(len(inputs_retrieved[0])):
        axs[0][i].set_yscale("log")
        axs[1][i].set_yscale("log")
        if scales[i]=="log":
            axs[0][i].set_xscale("log")
            axs[1][i].set_xscale("log")
            axs[2][i].set_xscale("log")

        axs[0][i].scatter(inputs_retrieved[:, i], paleo_volume[:, 0], c="grey", s=2)
        axs[1][i].scatter(inputs_retrieved[:, i], paleo_volume[:, 0]-paleo_volume[:, 1], c="grey", s=2)
        axs[2][i].scatter(inputs_retrieved[:, i], paleo_volume_fraction, c="grey", s=2)
        axs[2][i].set_xlabel(pars[i])
        axs[1][i].set_box_aspect(1)
        axs[0][i].set_box_aspect(1)
        axs[2][i].set_box_aspect(1)

        axs[0][i].set_xticklabels([])
        axs[1][i].set_xticklabels([])

        if i > 0:
            axs[0][i].set_yticklabels([])
            axs[1][i].set_yticklabels([])
            axs[2][i].set_yticklabels([])

    axs[0][0].set_ylabel(r"$V^{pre}_{paleo}$")
    axs[1][0].set_ylabel(r"$\Delta V^{100}_{paleo}$")
    axs[2][0].set_ylabel(r"$\frac{\Delta V^{100}_{paleo}}{V^{pre}_{paelo}}$")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def plot_final_toe_position(toe, inputs_retrieved, pars, scales):
    pass

def plot_abstraction_efficiency(fluxs, change, inputs_retrieved, pars, scales, units):
    """
    Inputs: fluxs and change: used to calculate efficiency
        inputs_retrieved: last is the post_development head, which is the independent variabe
    """

    
    f, axs, = plt.subplots(3, 1, figsize=(4, 6), sharex=False, layout="tight")
    for i in range(len(inputs_retrieved[0])-1):
        # sc = axs[i//4][i%4].scatter(inputs_retrieved[:, i], np.asarray(change)*inputs_retrieved[:, -2])
        sc = axs[i].scatter(inputs_retrieved[:, -1], inputs_retrieved[:,i], c=((np.asarray(fluxs))-(np.asarray(change)))/(np.asarray(change)), cmap="cool", s=5)
        axs[i].set_ylabel(f"{pars[i]} [{units[i]}]")
        if scales[i] =="log":
            axs[i].set_yscale("log")
    axs[2].scatter(inputs_retrieved[:, 0], 1-((np.asarray(fluxs))+(np.asarray(change)))/(np.asarray(change)),s=5)
    axs[2].set_xscale("symlog") 
    # axs[0][0].set_ylabel(r"$\frac{V^{100}_{extr}}{\Delta V^{100}_{ofg}}$")
    # axs[1][0].set_ylabel(r"$\frac{V^{100}_{extr}}{\Delta V^{100}_{ofg}}$")

    # cb = plt.colorbar(sc, ax=axs[1][3], cmap="viridis", extend="both", location="right")
    # cb.ax.set_title(pars[-3])
    
    # axs[1, 3].axis('off')
    plt.show()

def get_relationships(inputs, outputs):
    reg = LinearRegression().fit(inputs, outputs)
    coeff = reg.coef_
    return np.sign(coeff)

def plot_total_vs_paleo_abstraction(inputs_retrieved, paleo_volumes, total_volume_change, fraction_volume_change, paleo_volume_change):
    
    f, axs = plt.subplots(2, 1, figsize=(2.5, 5))
    paleo_flag = paleo_volumes[:, 0] > 0
    axs[0].scatter(inputs_retrieved[:, 0], inputs_retrieved[:, 1]/inputs_retrieved[:, 3], c=paleo_flag)
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[1].set_yscale("log")
    axs[1].scatter(inputs_retrieved[:, 4], total_volume_change, c=paleo_flag)
    axs[1].scatter(inputs_retrieved[:, 4][paleo_flag], paleo_volumes[:,0][paleo_flag]-paleo_volumes[:,1][paleo_flag], c="red")
    plt.show()

def Si_combined(problem, original_volume, total_volume_change, fractional_volume_change, pars, inputs_retrieved):
    Si = []
    Signs = []
    Si.append(sobol.analyze(problem, np.log10(np.asarray(original_volume))))
    Signs.append(get_relationships(inputs_retrieved, original_volume))
    Si.append(sobol.analyze(problem, np.log10(np.asarray(total_volume_change))))
    Signs.append(get_relationships(inputs_retrieved, total_volume_change))
    Si.append(sobol.analyze(problem, np.log10(np.asarray(fractional_volume_change))))
    Signs.append(get_relationships(inputs_retrieved, fractional_volume_change))

    f, axs = plt.subplots(2, 2, figsize=(8,7), layout="tight")
    plt.rcParams.update({'font.size': 9})
    for i in range(3):
        cmap = plt.get_cmap('tab10')
        colors = cmap(np.linspace(0, 1, len(pars)))
        markers = []
        for sign in Signs[i]:
            if sign == np.sign(1.0):
                markers.append("^")
            else: 
                markers.append("v")

        for j in range(len(colors)):
            axs[i//2][i%2].scatter(Si[i]["S1"][j], Si[i]["ST"][j], c=colors[j], label=pars[j], marker=markers[j], edgecolors="black", linewidths=0.5)
        
        axs[i//2][i%2].plot(np.linspace(0, 1 ,5), np.linspace(0, 1, 5), lw=0.5, ls="--", color="grey")
        axs[i//2][i%2].set_box_aspect(1)  
        # axs[i//2][i%2].legend()
        axs[i//2][i%2].set_xlabel(r"$S$")
        axs[i//2][i%2].set_ylabel(r"$ST$")
        axs[i//2][i%2].set_ylim([0, 1])
        axs[i//2][i%2].set_xlim([0, 1])
    legend_elements = [Line2D([0], [0], marker='o', color="w", markerfacecolor=colors[i], label=pars[i], markeredgecolor="black", markeredgewidth=0.5) for i in range(len(pars))]
    axs[0, 0].legend(handles=legend_elements, loc="center right", bbox_to_anchor=(0, 0.5), draggable=True)
    i_lower = np.tril_indices(len(pars), -1)
    Si[2]["S2"][i_lower] = Si[2]["S2"].T[i_lower]
    sns.heatmap(Si[2]["S2"],  vmin=0, cmap="pink_r", robust=False, annot=True, square=True, yticklabels=pars,xticklabels=pars, ax=axs[1][1], cbar=False, fmt=".2f")
    axs[0, 0].set_title("Predevelopment OFG volume")
    axs[0, 1].set_title("OFG volume salinized")
    axs[1, 0].set_title("Fraction OFG salinized")
    axs[1, 1].set_title("Interaction effects on fraction salinized")
    plt.show()

def Si_combined_seperate(problem, original_volume, total_volume_change, fractional_volume_change, pars, inputs_retrieved, axis=0, log=True):
    
    Si = []
    Signs = []
    if log == True:
        Si.append(sobol.analyze(problem, np.log10(np.asarray(original_volume))))
        Signs.append(get_relationships(inputs_retrieved, original_volume))
        Si.append(sobol.analyze(problem, np.log10(np.asarray(total_volume_change))))
        Signs.append(get_relationships(inputs_retrieved, total_volume_change))
        Si.append(sobol.analyze(problem, np.log10(np.asarray(fractional_volume_change))))
        Signs.append(get_relationships(inputs_retrieved, fractional_volume_change))
    else:
        Si.append(sobol.analyze(problem, np.asarray(original_volume)))
        Signs.append(get_relationships(inputs_retrieved, original_volume))
        Si.append(sobol.analyze(problem, np.asarray(total_volume_change)))
        Signs.append(get_relationships(inputs_retrieved, total_volume_change))
        Si.append(sobol.analyze(problem, np.asarray(fractional_volume_change)))
        Signs.append(get_relationships(inputs_retrieved, fractional_volume_change))

    if axis == 0:
        nrows=1
        ncols=3
        figsize=(8.5, 3)
        ileg=2
    else:
        nrows=3
        ncols=1
        figsize=(3, 8.5)
        ileg=1
    
    f, axs = plt.subplots(nrows, ncols, figsize=figsize, layout="tight", sharey=True)
    plt.rcParams.update({'font.size': 9})

    for i in range(3):
        cmap = plt.get_cmap('tab10')
        colors = cmap(np.linspace(0, 1, len(pars)))
        markers = []
        for sign in Signs[i]:
            if sign == np.sign(1.0):
                markers.append("^")
            else: 
                markers.append("v")

        for j in range(len(colors)):
            axs[i].scatter(np.max([Si[i]["S1"][j], 0]), np.max([Si[i]["ST"][j], 0]), c=colors[j], label=pars[j], marker=markers[j], edgecolors="black", linewidths=0.5)
        
        axs[i].plot(np.linspace(0, 1 ,5), np.linspace(0, 1, 5), lw=0.5, ls="--", color="grey")
        axs[i].set_box_aspect(1)  
        # axs[i//2][i%2].legend()
        axs[i].set_xlabel(r"$S$")
        if i == 0:
            axs[i].set_ylabel(r"$ST$")
        axs[i].set_ylim([0, 1])
        axs[i].set_xlim([0, 1])
    legend_elements = [Line2D([0], [0], marker='o', color="w", markerfacecolor=colors[i], label=pars[i], markeredgecolor="black", markeredgewidth=0.5) for i in range(len(pars))]
    legend_elements.append(Line2D([0], [0], marker='^', color="w", markerfacecolor="grey", label="+ effect", markeredgecolor="black", markeredgewidth=0.5))
    legend_elements.append(Line2D([0], [0], marker='v', color="w", markerfacecolor="grey", label="- effect", markeredgecolor="black", markeredgewidth=0.5))

    axs[ileg].legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.05, 0.5), draggable=True)
    # i_lower = np.tril_indices(len(pars), -1)
    # Si[2]["S2"][i_lower] = Si[2]["S2"].T[i_lower]
    # sns.heatmap(Si[2]["S2"],  vmin=0, cmap="pink_r", robust=False, annot=True, square=True, yticklabels=pars,xticklabels=pars, ax=axs[1][1], cbar=False, fmt=".2f")
    axs[0].set_title("Predevelopment OFG volume")
    axs[1].set_title("OFG volume salinized")
    axs[2].set_title("Fraction OFG salinized")
    plt.show()
    # plt.show()

def plot_effieciency_simple_effects(efficiency, total_flux, h, K, scales, lims):
    f, ax = plt.subplots(figsize=(4,4), layout="tight")
    plt.rcParams.update({'font.size': 9})
    active = np.array(total_flux) < 0
    sc = ax.scatter(h[active], K[active], c=efficiency[active], marker="o", alpha=0.4, s=5, cmap="viridis")
    sc = ax.scatter(h[~active], K[~active], c=efficiency[~active], marker="D", alpha=0.4, s=5, cmap="viridis")
    X1, X2, Y = fit_polynomial_regression(h, K, efficiency, scales, lims)
    ax.contour(X1, X2, Y, cmap="viridis") 
    ax.set_yscale("log")
    ax.set_ylabel(r"$K_b^V$ [m/day]")
    ax.set_xlabel(r"$h_f-h_{post}$ [m]")
    ax.set_box_aspect(1)  
    cb = plt.colorbar(sc, ax=ax, cmap="viridis", location="right", shrink=0.5)
    # ax.set_ylim(top=1000)
    plt.title("Abstraction Efficiency")
    plt.show()

def plot_efficiency_conductuvities(efficiency, x, y):
    f, ax = plt.subplots(figsize=(4,4), layout="tight")
    plt.rcParams.update({'font.size': 9})
    sc = ax.scatter(x, y, c=efficiency, marker="o", alpha=0.7, s=10, cmap="viridis", vmax=1) 
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"$K_a/n$")
    ax.set_ylabel(r"$K_b^V/n$")
    ax.set_box_aspect(1)  
    cb = plt.colorbar(sc, ax=ax, cmap="viridis", location="right", shrink=0.5)
    # ax.set_ylim(top=2700)
    plt.title("Abstraction Efficiency")
    plt.show()    


def Si_single(problem, result, name, pars, inputs_retrieved):

    

    Si = sobol.analyze(problem, np.asarray(result))
    Signs = get_relationships(inputs_retrieved, result)

    

    f, ax = plt.subplots(figsize=(4, 4), layout="tight")
    plt.rcParams.update({'font.size': 9})

    cmap = plt.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, len(pars)))
    markers = []
    for sign in Signs:
        if sign == np.sign(1.0):
             markers.append("^")
        else: 
            markers.append("v")

    for j in range(len(colors)):
        ax.scatter(Si["S1"][j], Si["ST"][j], c=colors[j], label=pars[j], marker=markers[j], edgecolors="black", linewidths=0.5)
        
    ax.plot(np.linspace(0, 1 ,5), np.linspace(0, 1, 5), lw=0.5, ls="--", color="grey")
    ax.set_box_aspect(1)  
        # axs[i//2][i%2].legend()
    ax.set_xlabel(r"$S$")
    ax.set_ylabel(r"$ST$")
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 1])
    legend_elements = [Line2D([0], [0], marker='o', color="w", markerfacecolor=colors[i], label=pars[i], markeredgecolor="black", markeredgewidth=0.5) for i in range(len(pars))]
    legend_elements.append(Line2D([0], [0], marker='^', color="w", markerfacecolor="grey", label="+ effect", markeredgecolor="black", markeredgewidth=0.5))
    legend_elements.append(Line2D([0], [0], marker='v', color="w", markerfacecolor="grey", label="- effect", markeredgecolor="black", markeredgewidth=0.5))

    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.05, 0.5), draggable=True)
    # i_lower = np.tril_indices(len(pars), -1)
    # Si[2]["S2"][i_lower] = Si[2]["S2"].T[i_lower]
    # sns.heatmap(Si[2]["S2"],  vmin=0, cmap="pink_r", robust=False, annot=True, square=True, yticklabels=pars,xticklabels=pars, ax=axs[1][1], cbar=False, fmt=".2f")
    ax.set_title(name)
    plt.show()


def plot_vertical_horizontal_conductivities(fraction, ka, kbv):
    # plot fraction salinized against ka and kvb
    # fit polynomial relationship
    f, ax = plt.subplots(figsize=(4,4), layout="tight")
    plt.rcParams.update({'font.size': 9})
    sc = ax.scatter(ka, kbv, c=fraction, marker="o", alpha=0.4, s=5, cmap="viridis")
    X1, X2, Y = fit_polynomial_regression(ka, kbv, fraction, ["log", "log"], [[1, 1000], [0.000001, 0.001]])
    ax.contour(X1, X2, Y, cmap="viridis") 
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"$K_a$ [m/day]")
    ax.set_ylabel(r"$K_b^V$ [m]")
    ax.set_box_aspect(1)  
    cb = plt.colorbar(sc, ax=ax, cmap="viridis", location="right", shrink=0.5)
    plt.title("Fraction OFG salinized")
    plt.show()

def plot_peclet_ka_h(h_range, k_range, pars):
    f, ax = plt.subplots()
    pass

def plot_time_scales(name, inputs, outputs):

    
    L = 30000
    D = 20
    z_aquifer = 35
    # delta_h_v = -(inputs[:, -2]-(0.025/1)*z_aquifer)
    # T_h = L**2*inputs[:, -1]/(delta_h_v*inputs[:, 0]*365)
    # T_v = -D**2*inputs[:, -1]/(inputs[:, -2]*inputs[:,1]*365) 
    # T_sal = 100/fraction_salinized
    T_h = -outputs["dl"]**2*inputs[:, -1]/(outputs["dhh"]*inputs[:, 0]*365)
    T_v = -outputs["dd"]**2*inputs[:, -1]/(outputs["dhv"]*inputs[:,1]*365)
    residual = [np.min([h, v]) for (h, v) in zip(T_h, T_v)]-outputs["T_sal"]
    active = [outputs["dhh"]<0]
    over = [residual > 0]
    f, ax = plt.subplots(figsize=(5,4), layout="tight")
    plt.rcParams.update({'font.size': 9})
    vmin = np.min(np.log10(outputs["T_sal"][active[0]]))
    vmax = np.max(np.log10(outputs["T_sal"][active[0]]))
    # norm = colors.SymLogNorm(0.5, vmin=min(np.log10(outputs["T_sal"][active[0]])), vmin=min(np.log10(outputs["T_sal"][active[0]]))vmax=1e4,)
    # colors_plot = cm.coolwarm_r(norm(outputs["T_sal"][active[0] & (inputs[:,2]<0.5)]))
    sc = ax.scatter(T_h[active[0] & over[0]] , T_v[active[0]& over[0]], c=np.log10(outputs["T_sal"][active[0]& over[0]]), cmap="coolwarm_r", norm=colors.TwoSlopeNorm(2, vmin=vmin, vmax=vmax), marker="^", s=40, edgecolors="w")
    ax.scatter(T_h[active[0]& ~over[0]] , T_v[active[0]& ~over[0]], c=np.log10(outputs["T_sal"][active[0]& ~over[0]]), cmap="coolwarm_r", norm=colors.TwoSlopeNorm(2, vmin=vmin, vmax=vmax), marker="v", s=40, edgecolors="w")

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
    cb.set_ticks([np.log10(50), np.log10(100), np.log10(200), np.log10(1000)])
    cb.set_ticklabels(["50","100","200", "1000"])
    ax.set_box_aspect(1)
    ax.set_xlim(right=1e6)
    ax.set_ylim(top=1e5)
    legend_elements = [Line2D([0], [0], marker='^', color="w", markerfacecolor="black", label="overestimate", markeredgecolor="black", markeredgewidth=0.5),
                       Line2D([0], [0], marker='v', color="w", markerfacecolor="black", label="underestimate", markeredgecolor="black", markeredgewidth=0.5)]
    ax.legend(handles=legend_elements, ncol=2)
    ax.plot(np.linspace(1e2,1e4), np.linspace(1e2,1e4), lw=1, ls="--", c="k")
    plt.show()

def plot_paleo_fraction_salinization(name, inputs, pv, total_total_volume=None):
    f, axs = plt.subplots(figsize=(4,4), layout="tight")
    plt.rcParams.update({'font.size': 9})
    vmin = 0
    vmax = np.max(pv[:, 0])
    if total_total_volume==None:
        result = pv[:, 0]
    else: 
        result = [max([0, frac]) for frac in  pv[:,0]/total_total_volume]
    sc = axs.scatter(inputs[:, 0], inputs[:, 1], c=result, alpha=0.7, s=10, norm=colors.SymLogNorm(0.001, vmax=1))
    axs.set_xscale("log")
    axs.set_yscale("log")
    X1, X2, Y = fit_polynomial_regression(inputs[:, 0], inputs[:, 1], result, ["log", "log"], [[1, 100],[0.000001, 0.001]])
    axs.contour(X1, X2, Y, cmap="viridis", norm=colors.SymLogNorm(0.001, vmax=1)) 
    axs.set_box_aspect(1)  
    cb = plt.colorbar(sc, ax=axs, cmap="viridis", location="right", shrink=0.5)
    # ax.set_ylim(top=1000)
    plt.title("Fraction Paleo OFG")
    axs.set_ylabel(r"$K_b^V$ [m/day]")
    axs.set_xlabel(r"$K_a$ [m/day]")
    # vmax = np.max(pv[:, 0])
    # axs[1].scatter(inputs[:, 0], inputs[:, 1], c=pv[:, 0]-pv[:, 1], norm=colors.SymLogNorm(100))
    # axs[1].set_xscale("log")
    # axs[1].set_yscale("log")
    # ax.set_yscale("symlog", linthresh=0.01)
    plt.show()
    return 