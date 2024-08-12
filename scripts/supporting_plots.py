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
from results_mf6 import save_results_for_saltelli, plot_color_mesh_saltelli, get_pars, get_outputs, paleo_volumes


def compare_paleo_models(names):

    f, axs = plt.subplots(1, len(names), figsize=(4, 8))
    for i in range(len(names)):
        inputs = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{names[i]}.csv", delimiter=",")
        pv = np.array(paleo_volumes(names[i], len(inputs)))
        total_total_volume_change, original_total_volume = results_mf6.extract_pre_post_volume(names[i], len(inputs))
        result = [max([0, frac]) for frac in  pv[:,0]/original_total_volume]
        sc = axs[i].scatter(inputs[:, 0], inputs[:, 1], c=result, alpha=0.7, s=10, norm=colors.SymLogNorm(0.001, vmax=1))
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")
        axs[i].set_box_aspect(1)  
        cb = plt.colorbar(sc, ax=axs[i], cmap="viridis", location="right", shrink=0.5)
        # ax.set_ylim(top=1000)
        axs[i].title.set_text(names[i])
        axs[i].set_ylabel(r"$K_b^V$ [m/day]")
        axs[i].set_xlabel(r"$K_a$ [m/day]")

    plt.show()