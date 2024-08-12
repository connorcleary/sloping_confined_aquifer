import numpy as np
from pars import load_parameters
import matplotlib.pyplot as plt
import pdb
from sklearn.linear_model import LinearRegression

def extract_results(ensemble_names):

    for name in ensemble_names:
        results = []
        for i in range(50):# len(os.listdir(f"/home/superuser/sloping_confined_aquifer/results/{name}/"))):
            output = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}.csv", delimiter=",")
            results.append(output)
            np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", np.array(results), delimiter=",")
        
        
def extract_fluxes(ensemble_names):

    
    for name in ensemble_names:
        results = []
        for i in range(50):
            pars = load_parameters(f"{name}{i}", backup=True)
            qz = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_qz.csv", delimiter=",")
            flux=[]
            #pdb.set_trace()
            for cell in pars.top_cells:
                if cell[2] >= pars.x_b/pars.dx:
                    flux.append(qz[cell[0]][cell[2]])
            results.append(flux)
        results_np = np.zeros((len(results), len(max(results, key=len))))
        for i, result in enumerate(results):
            results_np[i][:len(result)]=result
        np.savetxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all_fluxes.csv", results_np, delimiter=",")
    # pdb.set_trace()
    

def box_and_whiskers(ensemble_names):
    colors = ["blue", "orange"]
    f, ax = plt.subplots()
    averages=[[]]*len(ensemble_names)
    for n, (name, color) in enumerate(zip(ensemble_names, colors)):
        results= np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all_fluxes.csv", delimiter=",")
        # pdb.set_trace()
        for i, fluxes in enumerate(results):
            pars = load_parameters(f"{name}{i}", backup=True)
            #ax.plot(pars.dx*np.linspace(0, pars.L/pars.dx, len(fluxes)), pars.alpha_L/100*fluxes, c=color, alpha=0.3, lw=0.2)
            averages[n].append(np.average([flux for flux in fluxes if flux > 0])/pars.porosity*pars.alpha_L/100)
        
    # pdb.set_trace()
    ax.boxplot(averages, vert=True, labels=["chd", "ghb"])
    ax.set_ylabel(r"$\bar{\frac{\alpha_z q_z}{n}}$")
    ax.set_xlabel("Boundary type")
    pars = load_parameters(f"{name}{25}", backup=True)
    ax.axhline(pars.diff, c="red", ls="--", label="Molecular Diffusion")
    ax.set_yscale("log")
    ax.set_ylim(1e-6)
    ax.legend()

    plt.show()   
    
def scatter_postpre(ensemble_names):
    f, ax = plt.subplots()
    x=[]
    y=[]
    colors=[]
    colors_list = ["green", "red"]
    labels_list = ["chd", "ghb"]
    labels=[]
    
    for n, (name, color) in enumerate(zip(ensemble_names, colors_list)):
        some_x = list(np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", delimiter=",")[:,0])
        some_y = list(np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", delimiter=",")[:,1])
        ax.scatter(some_x, some_y, c=color, label=labels_list[n])
    ax.plot(np.linspace(0, 700000), np.linspace(0, 700000), ls="--", c="black", zorder=-1)
    plt.legend()
    plt.show()
    
def pairs_plot(results, scales, pars, units, lims=None, title=None):
    f, axs = plt.subplots(len(units), len(units), figsize=(6,6))
    #axs[0, 1].axis('off')
    
    for row in range(len(units)):
        for col in range(len(units)):
            if row>col:
                # pdb.set_trace()                
                scatter = axs[row, col].scatter(results[:50, col].T, results[:50, row].T, c=results[:50, -1].T, marker="o",  cmap="viridis", alpha=0.7, vmax=0)
                axs[row, col].scatter(results[50:, col].T, results[50:, row].T, c=results[50:, -1].T, marker="v",  cmap="viridis", alpha=0.7, vmax=0)
               
                if scales[row] == "log":
            	    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                axs[row, col].set_box_aspect(1)               
            elif row == col:
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")                 
                axs[row, col].set_ylabel(f"{title}")
                axs[row, col].set_box_aspect(1)               
                axs[row, col].scatter(results[:50, col].T, results[:50, -1].T, c="grey", marker="o", cmap="viridis", vmin=0, vmax=1, alpha=0.7)
                axs[row, col].scatter(results[50:, col].T, results[50:, -1].T, c="grey", marker="v", cmap="viridis", vmin=0, vmax=1, alpha=0.7)
            else:
                axs[row, col].axis('off')
            if col == 0 and row > 0:
            	axs[row, col].set_ylabel(f"{pars[row]} [{units[row]}]")
            if row == 1:
            	axs[row, col].set_xlabel(f"{pars[col]} [{units[col]}]")
           
         	
    cb = plt.colorbar(scatter, ax=axs[0, 1], cmap="viridis", extend="both", location="right")
    cb.ax.set_title(f"{title}")
    cb.ax.yaxis.set_ticks_position("left")
    #plt.tight_layout()       
   
    plt.show()
    
def make_pairs_plots(ensemble_names):
    metric_names = [r"log$_{10}(Di)$ [-]", r"$\frac{\Delta V}{V_{pre}}$ [-]", r"$\frac{V_{pre}}{V_T}$"]
    diff_metric = []
    volume_metric = []
    preservation_metric = []
    markers = []
    # pdb.set_trace()
    for name, marker_style in zip(ensemble_names, ["o", "v"]):
        results0= np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all_fluxes.csv", delimiter=",")
        results1 = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all.csv", delimiter=",")   
        for i, (fluxes, volumes) in enumerate(zip(results0, results1)):
            index = i % 50  
            pars = load_parameters(f"{name}{index}", backup=True)
            markers.append(marker_style)
            diff_metric.append(np.log10(np.average([flux for flux in fluxes if flux > 0])/pars.porosity*pars.alpha_L/100/pars.diff))
            if volumes[0] == 0:
                preservation_metric.append(0)
                volume_metric.append(0)
            else:
                volume_metric.append((volumes[0]-volumes[1])/volumes[0])
                preservation_metric.append(volumes[0]/((pars.D+pars.H)*pars.L))
    # pdb.set_trace()
    scales = ["linear", "linear", "log", "log", "linear", "linear", "linear", "linear", "linear", "linear"]
    pars = ["H", "D", r"$K_a$", "$K_b$", r"$\alpha$", "anis", r"$\beta$", "L", r"$h_{pre}$", r"$\Delta h$"]
    units = ["m", "m", "m/day", "m/day", "-", "-", "radians", "m", "m", "m"]
    
    inputs = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{ensemble_names[0]}.csv", delimiter=",")
    inputs = np.vstack([inputs, inputs])
    inputs0 = np.column_stack([inputs, np.array(diff_metric).T])
    inputs1 = np.column_stack([inputs, np.array(volume_metric).T])
    inputs2 = np.column_stack([inputs, np.array(preservation_metric).T])
    Gamma_aquifer = 0.3*120/(20000*365)*1000/(inputs[:, 2]*1025*np.tan(inputs[:, 6])**2)
    Gamma_aquitard = 0.3*120/(20000*365)*1000/(inputs[:, 3]*1025*np.tan(inputs[:, 6])**2)
    inputs0_1 = np.column_stack([inputs0[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs0[:, -1]])
    inputs1_1 = np.column_stack([inputs1[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs1[:, -1]])
    inputs2_1 = np.column_stack([inputs2[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs2[:, -1]])
    scales_1 = ["log", "log", "linear", "linear", "log", "log"]
    pars_1 = [r"$K_a$", "$K_b$", r"$\alpha$", "anis", r"$\Gamma_a$", r"$\Gamma_b$"]
    units_1 = ["m/day", "m/day", "-", "-", "-", "-",]
    
    inputs0_1 = np.column_stack([inputs0[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs0[:, -1]])
    inputs1_1 = np.column_stack([inputs1[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs1[:, -1]])
    inputs2_1 = np.column_stack([inputs2[:, [2, 3, 4, 5]], Gamma_aquifer, Gamma_aquitard, inputs2[:, -1]])
    
    K_a_alpha = inputs[:, 2]*inputs[:,4]
    K_b_alpha = inputs[:, 3]*inputs[:,4]
    
    inputs0_2 = np.column_stack([K_a_alpha, K_b_alpha, inputs0[:, 5], Gamma_aquifer, Gamma_aquitard, inputs0[:, -1]])
    inputs1_2 = np.column_stack([K_a_alpha, K_b_alpha, inputs1[:, 5], Gamma_aquifer, Gamma_aquitard, inputs1[:, -1]])
    inputs2_2 = np.column_stack([K_a_alpha, K_b_alpha, inputs2[:, 5], Gamma_aquifer, Gamma_aquitard, inputs2[:, -1]])
    
    # scales_2 = ["log", "log", "linear", "log", "log"]
    # pars_2 = [r"$K_a \alpha$", r"$K_b \alpha$", "anis", r"$\Gamma_a$", r"$\Gamma_b$"]
    # units_2 = [r"$m^2/day$", r"$m^2/day$", "-", "-", "-",]
    # f, ax = plt.subplots()
    # ax.boxplot([inputs2_2[:50, -1],inputs2_2[50:, -1]], vert=True, labels=["chd", "ghb"])
    # ax.set_ylabel(metric_names[2])
    # ax.set_xlabel("Boundary type")
    # pars = load_parameters(f"{name}{25}", backup=True)
    # ax.legend()
    
    f, ax = plt.subplots()
    sc = ax.scatter(preservation_metric[:50],volume_metric[:50], c = diff_metric[:50], cmap="viridis", marker="o", vmax=0)
    ax.scatter(preservation_metric[50:],volume_metric[50:], c = diff_metric[50:], cmap="viridis", marker="v", vmax=0)
    ax.set_xlabel(r"$\frac{V_f^{pre}}{V_T}$ [-]", fontsize="large")
    ax.set_ylabel(r"$\frac{\Delta V_{f}^{100}}{V_f^{pre}}$ [-]", fontsize="large")
    cb = plt.colorbar(sc, cmap="viridis", extend="both", location="right")
    cb.ax.set_title(r"log$_{10}\left(\frac{\alpha_v \bar{q_z}}{nD_m}\right)$ [-]")
    plt.tight_layout()
    plt.show()
    #K_a_alpha_anis = *inputs[:,4]/inputs[:, 5]
    K_b_alpha_anis_D = inputs[:, 3]*inputs[:,4]/inputs[:, 5]/inputs[:, 1]
    inputs0_3 = np.column_stack([inputs[:, 2], K_b_alpha_anis_D,inputs0[:, -1]])
    scales_3 = ["log", "log"]
    pars_3 = [r"$K_a $", r"$\frac{K_b \alpha}{anis_V D}$"]
    units_3 = [r"$m/day$", r"$m^2/day$"]
    f, ax = plt.subplots()
    sc = ax.scatter(inputs[:50, 2], K_b_alpha_anis_D[:50], c = diff_metric[:50], cmap="viridis", marker="o", vmax=0)
    ax.scatter(inputs[50:, 2], K_b_alpha_anis_D[50:], c = diff_metric[50:], cmap="viridis", marker="v", vmax=0)
    ax.set_xlabel(r"$K_a$ [m/day]", fontsize="large")
    ax.set_ylabel(r"$\frac{K_b \alpha_L}{anis_V D}$ [m$^2$/day]", fontsize="large")
    cb = plt.colorbar(sc, cmap="viridis", extend="both", location="right")
    cb.ax.set_title(r"log$_{10}\left(\frac{\alpha_v \bar{q_z}}{nD_m}\right)$ [-]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.show()
    
    # pairs_plot(inputs0_3, scales_3, pars_3, units_3, title=metric_names[0])
    # pairs_plot(inputs1_2, scales_2, pars_2, units_2, title=metric_names[1]) 
    # pairs_plot(inputs2_2, scales_2, pars_2, units_2, title=metric_names[2])      

def linear_regression(ensemble_names):

    diff_metric = []
    for name in ensemble_names:
        results0= np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/results/{name}/all_fluxes.csv", delimiter=",")
        for i, fluxes in enumerate(results0):
            index = i % 50  
            pars = load_parameters(f"{name}{index}", backup=True)
            diff_metric.append(np.log10(np.average([flux for flux in fluxes if flux > 0])/pars.porosity*pars.alpha_L/100/pars.diff))
    inputs = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{ensemble_names[0]}.csv", delimiter=",")
    inputs = np.vstack([inputs, inputs])
    Gamma_aquitard = 0.3*120/(20000*365)*1000/(inputs[:, 3]*1025*np.tan(inputs[:, 6])**2)
    inputs0 = np.column_stack([inputs, np.array(diff_metric).T])
    K_a_alpha = inputs[:, 2]*inputs[:,4]
    K_b_alpha = inputs[:, 3]*inputs[:,4]
    inputs0_0 = np.column_stack([np.log10(K_a_alpha), np.log10(K_b_alpha), inputs0[:, -1]])
    
    reg = LinearRegression().fit(inputs0_0[:, :-1], inputs0_0[:,-1])
    outputs = reg.predict(inputs0_0[:, :-1])
    residuals = outputs-inputs0_0[:,-1]
    outputs_binary = []
    for  i, (output, true) in enumerate(zip(outputs, inputs0_0[:,-1])):
        if output < -0.5 and true > 0:
            outputs_binary.append(2)
        elif output < -0.5:
            outputs_binary.append(0)
        else:
            outputs_binary.append(1)
    
    f, ax = plt.subplots()
    sc = ax.scatter(inputs0_0[:50,0], inputs0_0[:50,1], c=outputs_binary[:50], marker="o")
    ax.scatter(inputs0_0[50:,0], inputs0_0[50:,1], c=outputs_binary[50:], marker="v")# _binary_with_incorrect)
    cb = plt.colorbar(sc, ax=ax, cmap="viridis", extend="both", location="right")
    plt.show()
    print(f"D = {reg.intercept_} + a({reg.coef_[0]}K_a+{reg.coef_[1]}K_b)")
    #print(reg.score(inputs0_0[:, :-1], inputs0_0[:,-1]))
    
     
def main():
    # extract_results(["constant_c", "ghb"])
    # extract_fluxes(["constant_c", "ghb"])
    # pars=extract_pars(["constant_c", "ghb])
    # box_and_whiskers(["constant_c", "ghb"])
    #scatter_postpre(["constant_c", "ghb"])
    # corner_volume_change(results, ensemble_name)
    # corner_disp_vs_diff(results, ensemble_name)
    make_pairs_plots(["constant_c", "ghb"])
    # linear_regression(["constant_c", "ghb"])
    
if __name__=="__main__":
    main()
