from scipy.stats import qmc
import pdb
import numpy as np
import os
import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool
import post_processing as proc
import results_mf6
from results_mf6 import save_results_for_saltelli, plot_color_mesh_saltelli, get_pars, get_outputs, paleo_volumes
from SALib.sample import sobol as sbsample
from SALib.analyze import sobol
import pairs_plot
import paper_plots
import pickle
import supporting_plots as sp

def produce_samples(name, n, d, names, ranges, scales, backup=False):
    
    if backup==True:
        print("loading old problem")
        param_values_unlog = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/backups/{name}.csv", delimiter=",")
        file = open(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/backups/{name}.pkl", 'rb')
        problem = pickle.load(file)
    else:     
        print("creating new problem")
        ranges_log = []
        for i in range(len(ranges)): 
            if scales[i] == "log":
                ranges_log.append([np.log10(ranges[i][0]), np.log10(ranges[i][1])]) 
            else:
                ranges_log.append(ranges[i])

        problem = {
        'num_vars': d,
        'names': names,
        'bounds': ranges_log,
        }   

        param_values = sbsample.sample(problem, n)
        param_values_unlog = []
        for sample in param_values:
            sample_params = []
            for i in range(len(ranges)):
                if scales[i] == "log":
                    sample_params.append(10**sample[i]) 
                else:
                    sample_params.append(sample[i])
            param_values_unlog.append(sample_params)

        np.savetxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", param_values_unlog, 
    		    delimiter=",")
        file = open(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.pkl", 'wb')
        pickle.dump(problem, file)


    return param_values_unlog, problem

def model_run(inputs):
    ensemble_name, name, real = inputs
    pars = ModelParameters(name, H=real[0], D=real[1], K_aquifer=real[2], 
                           K_aquitard=real[3], alpha_L=real[4], anis_aquifer=real[5], 
                           z0=np.tan(real[6])*real[7]+real[0]+real[1], L=real[7], 
                           h_modern=real[8], porosity=real[9], Ss=real[10], outputs="all")    
    sim = mf6.build_steady_model(pars)
    succ = True
    try:
        mf6.run_model(sim)
    except:
        succ= False
    
    if succ: 
        save_results_for_saltelli(name, ensemble_name)

    mf6.remove_model_files(name)   

def model_run_sens(inputs):
    # function to run model with fixed geometry and variable hydraulic gradient
    ensemble_name, name, real = inputs
    pars = ModelParameters(name, H=20, D=20, K_aquifer=real[0], 
                        K_aquitard=real[1]*real[3], alpha_L=real[2], anis_aquifer=real[3], anis_aquitard=real[3],
                        z0=55, L=30000, x_b=4000,  
                        h_modern=real[4], porosity=real[-1], gradient=1, Ss=1e-6, outputs="all", dx=5, dz=0.5, dt=100)    
    sim = mf6.build_steady_model(pars)
    succ = True
    try:
        mf6.run_model(sim)
    except:
        succ= False
    
    # save_results_for_saltelli(name, ensemble_name)

    if succ: 
        try:
            save_results_for_saltelli(name, ensemble_name)
        except: 
            pass
        
    # mf6.remove_model_files(name)  
    #plot_color_mesh_saltelli(name, int(name.replace(ensemble_name, "")))


def rerun(name, inputs, new_name=None):
    inputs_new = []
    for i, sample in enumerate(inputs):
        if not os.path.isfile(f"/home/superuser/sloping_confined_aquifer/results/{name}/{name}{i}_metrics.csv"):
            inputs_new.append([name, f"{name}{i}", sample])
    p=Pool(processes=4)
    p.map(model_run_sens, inputs_new)
    #model_run_sens(inputs_new[0])

def main_rerun():
    name = "draft2"
    inputs_retrieved = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", delimiter=",")

def main_analysis():
    # name = "draft1g"
    # ranges = [[1, 100],
    #       [0.000001, 0.001],
    #       [0.1, 10],
    #       [1, 100],
    #       [1, 100],
    #       [-2, 2],
    #       [0.2,0.8]]
    # scales = ["log", "log", "log", "log", "log", "linear", "linear"] # logs
    # pars = [r"$K_a$", r"$K_b^V$", r"$\alpha$", r"$anis_a$", r"$anis_b$", r"$h_{post}$", r"$n$"]
    # units = ["m/day", "m/day", "m", "-", "-", "m", "-"]
    # param_values, problem = produce_samples(name, 4, 7, pars, ranges, scales)
    # name = "paper3new6"
    # ranges = [[1, 100],
    #       [1e-5, 1e-3],
    #       [0.1, 10],
    #       [1, 100],
    #       [-10, 2],
    #       [0.2,0.8]]
    # scales = ["log", "log", "log", "log", "linear", "linear"] # logs
    # pars = [r"$K_a$", r"$K_b^V$", r"$\alpha$", r"$anis$", r"$h_{post}$", r"$n$"]
    # units = ["m/day", "m/day", "m", "-", "m", "-"]
    name = "paper3new7"
    ranges = [[1, 100],
          [1e-5, 1e-3],
          [0.1, 10],
          [1, 100],
          [-10, 2],
          [0.2,0.8]]
    scales = ["log", "log", "log", "log", "linear", "linear"] 
    pars = [r"$K_a$", r"$K_b^V$", r"$\alpha$", r"$anis$", r"$h_{modern}$", r"$n$"]
    units = ["m/day", "m/day", "m", "-", "m", "-"]
    param_values, problem = produce_samples(name, 64, 6, pars, ranges, scales, backup=True)
    

    inputs_retrieved = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", delimiter=",")
    # inputs_retrieved = np.asarray(get_pars(name, 896))
    fractional_volume_change, total_volume_change, mixing, original_volume, total_flux, predev_total_flux = get_outputs(name, 896, inputs_retrieved)
    # original_total_volume, total_total_volume_change = results_mf6.extract_pre_post_volume(name, 1152)
    # h = -(inputs_reetrieved[:, -2]-(0.025/1)*35)
    # total_fraction = [total_total_volume_change[i]/original_total_volume[i] for i in range(1152)]
    # paper_plots.plot_vulnerability(inputs_retrieved, h, fractional_volume_change)
    # paper_plots.Si_combined_seperate(problem, original_volume, total_volume_change, fractional_volume_change, pars, inputs_retrieved, log=False)
    # indices = [0, 1, 3, 4]
    #pairs_plot.pairs_plot_3(inputs_retrieved[:, [0, 1, 3, 4]], fractional_volume_change, [scales[x] for x in indices], [pars[x] for x in indices], [units[x] for x in indices], [ranges[x] for x in indices])
    # pairs_plot.pairs_plot_3(inputs_retrieved[:, :4], fractional_volume_change, scales[:4], pars[:4], units[:4], ranges[:4])

    # paper_plots.fit_and_plot_regression(name, pars, inputs_retrieved, fractional_volume_change)
    # exit()
    # pv = np.array(paleo_volumes(name, 896))
    # total_total_volume_change, original_total_volume = results_mf6.extract_pre_post_volume(name, 1024)
    # paper_plots.plot_paleo_fraction_salinization(name, inputs_retrieved, pv, np.array(original_total_volume))
    # outputs = results_mf6.get_outputs_for_timescales(name, 896, inputs_retrieved)
    # paper_plots.plot_time_scales(name, inputs_retrieved, outputs)
   
    # original_total_volume, total_total_volume_change = results_mf6.extract_pre_post_volume(name, 896)
    # efficiency = (np.asarray(predev_total_flux)-np.asarray(total_flux))/(np.asarray(np.asarray(total_total_volume_change)))
    # conduc_ave = np.log10(20*inputs_retrieved[:, 0])+np.log10(30000*inputs_retrieved[:,1])
    # paper_plots.plot_effieciency_simple_effects(efficiency, total_flux, inputs_retrieved[:, -2], conduc_ave, ["linear", "linear"], [[np.min(inputs_retrieved[:,-2]), np.max(inputs_retrieved[:,-2])],[np.min(conduc_ave), np.max(conduc_ave)]])
    # exit()
    # ranges = [[1, 1000],
    #       [0.0001, 0.001],
    #       [0.1, 10],
    #       [1, 100],
    #       [-10, 1],
    #       [0.2,0.8],
    #       [1, 4]]
    # scales = ["log", "log", "log", "log", "linear", "linear", "log"] # logs
    # pars = [r"$K_a$", r"$K_b$", r"$\alpha$", r"$anis$", r"$h_{post}$", r"$n$", r"$I_{pre}$"]
    # name = "draft4"
    # units = ["m/day", "m/day", "m", "-", "m", "-", "-"]
    # param_values, problem = produce_samples(name, 64, 7, pars, ranges, scales)

    # units = ["m/day", "m/day", "m", "-", "-", "m", "-"]
    # param_values, problem = produce_samples(name, 64, 7, pars, ranges, scales) # 8
    # inputs = [[name, f"{name}{i}", sample] for i, sample in enumerate(param_values)]
    # inputs_retrieved = np.asarray(get_pars(name, 1024))    # rerun(name, inputs_retrieved)
    # fractional_volume_change, total_volume_change, mixing, original_volume, total_flux, predev_total_flux = get_outputs(name, 1024, inputs_retrieved)
    # fractional_volume_change = np.array(fractional_volume_change)
    # total_volume_change = np.array(total_volume_change)
    # fractional_volume_change[fractional_volume_change <= 0] = 0.0001 
    # total_volume_change[total_volume_change <= 0] = 0.0001
    # pv = np.array(paleo_volumes(name, 1024))
    # fraction_paleo = 1-(original_volume-pv[:,0])/original_volume
    # volume_paleo = pv[:,0]
    # paleo_fraction_salinized = np.array((pv[:,0] - pv[:,1])/pv[:,0])
    # paleo_fraction_salinized[np.isnan(paleo_fraction_salinized)] = 0
    # total_total_volume_change, original_total_volume = results_mf6.extract_pre_post_volume(name, 1024)
    # paper_plots.plot_paleo_fraction_salinization(name, inputs_retrieved, pv, np.array(original_total_volume))
    # exit()
    # outputs = results_mf6.get_outputs_for_timescales(name, 1024, inputs_retrieved)
    # paper_plots.plot_time_scales(name, inputs_retrieved, outputs)
    pv = np.array(paleo_volumes(name, 896))
    fraction_paleo = 1-(original_volume-pv[:,0])/original_volume
    volume_paleo = pv[:,0]
    paleo_fraction_salinized = np.array((pv[:,0] - pv[:,1])/pv[:,0])
    paleo_fraction_salinized[np.isnan(paleo_fraction_salinized)] = 0
    #aper_plots.Si_combined_seperate(problem, original_volume, total_volume_change, fractional_volume_change, pars, inputs_retrieved)
    paper_plots.Si_combined_seperate(problem, fraction_paleo, volume_paleo, paleo_fraction_salinized, pars, inputs_retrieved, log=False)
    # paper_plots.Si_single(problem, (original_volume-pv[:,0])/original_volume, r"Fraction $OFG_{paleo}$ salinized", pars, inputs_retrieved)
    # plot_color_mesh_saltelli(name, 1024)
    # original_total_volume, total_total_volume_change = np.asarray(results_mf6.extract_pre_post_volume(name, 896))
    # efficiency = 1-((np.asarray(total_flux))+(np.asarray(total_total_volume_change)))/(np.asarray(total_total_volume_change))
    # efficiency = (np.asarray(predev_total_flux)-np.asarray(total_flux))/(np.asarray(total_total_volume_change))
    # efficiency[efficiency != efficiency] = 0
    # h = -(inputs_retrieved[:, -2]-(0.025/1)*30)
    # paper_plots.plot_effieciency_simple_effects(efficiency, total_flux, inputs_retrieved[:, -2], inputs_retrieved[:, 0], ["linear", "log"], [[-10, 2],[1, 100]])


def main_sensitivity():
    # Main function for running sensitivity analysis.
    # The difference between this and other main function is that aquifer geometry is constant
    # An added sensitivity is to the predevelopment hydraulic gradient
    name = "paper5"
    ranges = [[1, 100],
          [1e-5, 1e-3],
          [0.1, 10],
          [1, 100],
          [-10, 2],
          [0.2,0.8]]
    scales = ["log", "log", "log", "log", "linear", "linear"] 
    pars = [r"$K_a$", r"$K_b^V$", r"$\alpha$", r"$anis$", r"$h_{post}$", r"$n$"]
    units = ["m/day", "m/day", "m", "-", "m", "-"]
    param_values, problem = produce_samples(name, 4, 6, pars, ranges, scales) # 8
    inputs = [[name, f"{name}{i}", sample] for i, sample in enumerate(param_values)]
    # inputs_retrieved = np.asarray(get_pars(name, 1024))    # rerun(name, inputs_retrieved)
    # p=Pool(processes=2)
    # p.map(model_run_sens, inputs)
    model_run_sens(["test5", "test50", [100, 1e-5, 10, 10, 0, 0.5]])
    # inputs_retrieved = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", delimiter=",")
    # rerun(name, inputs_retrieved)
    # model_run_sens(inputs[0])
    # plot_color_mesh_saltelli(name, 1024)
    # paper_plots.plot_rate_of_intrusion(name, 1024)
    # pass
    # fractional_volume_change, total_volume_change, mixing, original_volume, total_flux = get_outputs(name, 1024, inputs_retrieved)
    # paper_plots.fit_and_plot_regression(name, pars, inputs_retrieved, fractional_volume_change)
    # total_total_volume_change = np.asarray(results_mf6.extract_pre_post_volume(na, 1024))
    # pv = np.array(paleo_volumes(name, 1024))
    # pv_change = (pv[:, 0]-pv[:, 1])/np.array([max([pvi[0], 0.001]) for pvi in pv])

    # paper_plots.plot_total_vs_paleo_abstraction(inputs_retrieved, pv, total_volume_change, fractional_volume_change, (pv[:,0]-pv[:, 1])/pv[:, 0])
    # plot_paleo_volume_and_changes(pv, inputs_retrieved, pars, scales)
    # toe = toe_position(name, 256)
    # plot_original_volume_and_total_change(original_volume, total_volume_change, inputs_retrieved, pars, scales)
    # plot_color_mesh_saltelli(name, 1024)
    # paper_plots.Si_combined_seperate(problem, original_volume, total_volume_change, fractional_volume_change, pars, inputs_retrieved)
    # inputs_retrieved_change = inputs_retrieved[:, [0, 1, 3, 4]]
    # inputs_retrieved_eff =  inputs_retrieved[:, [0, 1, 4, 6]]
    # horizontal conductivity times dispersivity
    # inputs_retrieved_eff[:, 0] = inputs_retrieved_eff[:, 0]/inputs_retrieved[:, -2]/inputs_retrieved[:, 2]*(inputs_retrieved_eff[:, -1]-(0.025/1)*30)*inputs_retrieved[:, 3]/inputs_retrieved[:, 1]
    # forces opposing free convection in aquifer
    # inputs_retrieved_eff[:, 1] = inputs_retrieved[:, 1]/inputs_retrieved[:, 3]
    # total_total_volume_change = total_total_volume_change*inputs_retrieved[:, -2]
    # pars_eff = [r"$K_a$", r"$K_b^V$", r"$h_{post}$", r"$I$"]
    # scales_eff = ["log", "log", "linear", "log"]
    # units_eff = ["m/day", "m/day", "m", "-"]
    # ranges_eff = [[1, 1000], [0.000001, 0.001],[-10, 1], [1,4]]
    # efficiency = 1-((np.asarray(total_flux))+(np.asarray(total_total_volume_change)))/(np.asarray(total_total_volume_change))
    # results = np.column_stack([inputs_retrieved_eff, efficiency])
    # h = -(inputs_retrieved[:, -3]-(0.025/1)*30)
    # hori = inputs_retrieved[:,0]/inputs_retrieved[:,-2]
    # vert = inputs_retrieved_eff[:, 1]/inputs_retrieved[:, -2]
    # hori_disp = hori*inputs_retrieved[:,2]
    # vert_disp = vert*inputs_retrieved[:,2]
    # x = hori
    # y = vert
    # lims = [[np.min(x), np.max(x)],[np.min(y), np.max(y)]]
    # paper_plots.fit_polynomial_regression(x, y, efficiency, ["log","log"], lims)
    # find_plot_clusters("draft4", 1024, scales, pars, n_clusters=4)
    # all_single_similarity("draft4", 1024, scales, pars, 286, 100)
    # paper_plots.plot_effieciency_simple_effects(efficiency, h, inputs_retrieved[:,1]/inputs_retrieved[:,3], ["linear", "log"], [[-0.15, 11],[0.001, 0.000001]])
    # paper_plots.plot_efficiency_conductuvities(efficiency, hori, vert)
    # pairs_plot.pairs_plot_3(results, scales_eff, pars_eff, units_eff, ranges_eff)
    # inputs_retrieved_change[:,1] = inputs_retrieved_change[:,1]/inputs_retrieved_change[:,2]
    # pars_change = [r"$K_a$", r"$K_b^V$", r"$anis$", r"$h_{post}$"]
    # scales_change = ["log", "log", "log", "linear"]
    # units_change = ["m/day", "m/day", "-", "m"]
    # ranges_change = [[1, 1000], [0.000001, 0.001],[1, 100],[-10, 1]]
    # results = np.column_stack([np.asarray(inputs_retrieved_change), np.asarray(fractional_volume_change).T])
    # pairs_plot.pairs_plot_3(results, scales_change, pars_change, units_change, lims=ranges_change, title=None)
    # plot_original_volume_vs_volume_change(original_volume, fractional_volume_change, np.asarray(inputs_retrieved), pars, scales)
    # paper_plots.plot_abstraction_efficiency(total_flux, total_total_volume_change, inputs_retrieved_eff, pars_eff, scales_eff, units_eff)
    # paper_plots.plot_vertical_horizontal_conductivities(fractional_volume_change, inputs_retrieved[:, 0], inputs_retrieved[:, 1]/inputs_retrieved[:, 3])
    pass
        
def main():
    # name = "64saltelli"

    name = "32again"
    ranges = [[10, 50], 
          [10, 20], 
          [1, 1000], 
          [0.0001, 0.001], 
          [0.1, 10],
          [1, 100], 
          [1e-4, 2e-3], 
          [10000, 20000],
          [-10, 0.5],
          [0.2,0.8]]
    
    # scales = ["linear", "linear", "log", "log", "linear", "linear", "linear", "linear", "linear", "linear"]
    # pars = ["H", "D", r"$K_a$", "$K_b$", r"$\alpha$", "anis", r"$\beta$", "L", r"$h_{post}$", "n"]
    # units = ["m", "m", "m/day", "m/day", "m", "-", "radians", "m", "m", "-"]
    # param_values, problem = produce_samples(name, 32, 10, pars, ranges, scales)
    # inputs = [[name, f"{name}{i}", sample] for i, sample in enumerate(param_values)]
    # inputs_retrieved = get_pars(name, len(inputs))
    # # inputs_retrieved = np.genfromtxt(f"/home/superuser/sloping_confined_aquifer/ensemble_pars/{name}.csv", delimiter=",")

    # # rerun(name, inputs_retrieved)
    # # fresh_volume_change, mixing = get_outputs(name, len(inputs))
    # # sfit_and_plot_regression(name, pars, inputs_retrieved, fresh_volume_change)
    # # get_rate_of_intrusion(name, len(fresh_volume_change))
    # # Si = sobol.analyze(problem, np.log10(np.asarray(fresh_volume_change)))
    # # scatter_and_heatmap(Si, pars)
    # # Si2 = sobol.analyze(problem, np.asarray(mixing))
    # # scatter_and_heatmap(Si2, pars)
    # # p=Pool(processes=8)
    # # p.map(model_run, inputs)
    # # model_run(inputs[0])
    # # plot_color_mesh_saltelli(name, len(param_values))
    # # inputs_arr_raw = np.asarray(inputs_retrieved)
    # #inputs_arr = inputs_arr_raw[:, [1, 2, 4, 8]]
    # inputs_arr = np.column_stack([inputs_arr, inputs_arr_raw[:, 3]/inputs_arr_raw[:, 5]])
    # results = np.column_stack([inputs_arr, np.asarray(fresh_volume_change).T])
    # # results2 = np.column_stack([np.asarray(inputs_retrieved)[:, :], np.asarray(mixing).T])
    # scales = ["linear", "log", "linear", "linear", "log"]
    # pars = ["D", r"$K_a$", r"$\alpha$", r"$h_{post}$", r"$\frac{K_b}{anis}$"]
    # units = ["m", "m/day", "m", "m", "1/day"]
    # ranges = [[10, 20], 
    #       [1, 1000],  
    #       [0.1, 10],
    #       [-10, 2],
    #       [10**(-6), 10**(-3)]]
    # pairs_plot.pairs_plot(results, scales, pars, units, lims=ranges, title=None)
    # # pairs_plot.pairs_plot(results2, scales, pars, units, lims=None, title=None)

if __name__=="__main__":
    plot_color_mesh_saltelli("test5", 1)
    #main_sensitivity()
    # main_analysis()
    # sim = mf6.build_model_start_from_modern("test4", "1", "1b", new_head=-2)
    # mf6.run_model(sim)