import sloping_confined_aquifer as sca
from pars import ModelParameters, load_parameters
import results
import mf6
import results_mf6

def create_run_plot_model(name, **kwargs):
    """
        Create, run and plot a new scenario
        
        Inputs:5
            name: name for the model
            **kwargs: any key word arguments to be passed to the ModelParameters class
        Outputs:
            None
    """
    pars = ModelParameters(name, **kwargs)
    # # # swt = sca.build_steady_model(pars)
    # # # swt = sca.build_modern_model(name)
    swt = mf6.build_steady_model(pars)
    # print("running")
    mf6.run_model(swt)
    results_mf6.all_results(name)
        #results.plot_last_step(name)

    # concentration, head, qx, qy, qz = sca.extract_results(name)
    # results.plot_results(name)
    # results.plot_gif(name)
    # results.metrics(name)
    # results.plot_metrics(name)

    #results.plot_last_step(name+"_modern")
     
def main():
    create_run_plot_model("b2", L=10000, T=1000, dt=1000, z0=44, x_b=3500, h_b=1.9, dx=10, dz=1, h_modern=-2, H=36, D=5, T_modern=365*200, K_aquifer=3, dt_modern=365*200, alpha_L=5)
    # create_run_plot_model("canterbury_mid_test_fixed_head_tsmult", L=10000, T=1000, dt=1000, z0=44, x_b=3500, h_b=1.9, dx=10, dz=1, h_modern=-1, H=36, D=5, T_modern=365*200, K_aquifer=3, dt_modern=365*200)
    #create_run_plot_model("canterbury_full_dx_100_alpha_1_hb_0", dz=1, L=50000, z0=76,x_b = 1900, h_b=0, dx=100, dt=1e3, T=1e3, alpha_L=1)
    #create_run_plot_model("fine_diff", dx=200, dz=1, L=50000, diff=8.64e-6)
    #create_run_plot_model("canterbury_steady_low", steady=True, sea_level=-20)
    #create_run_plot_model("suriname", L=130000, H=85, z0=185, K_aquifer=25, x_b=5, h_b=3, K_aquitard=0.04, D=30, T_init=2e5, anis_aquifer=100)
    # create_run_plot_model("kooi_case15b", L=400, H=9, D=1, x_b=1, h_b = 0, z0=10.4, dx=0.5, dz=0.1, dt=1e3,T=1e4, K_aquifer=1e-1, K_aquitard=1e-5,
    #     alpha_L=1e-3, alpha_anisT=1, alpha_anisV=1, v_slr=0.1/1000/365, anis_aquifer=1)
    # create_run_plot_model("kooi_case17", L=400, H=9.5, D=0.6, x_b=1, h_b = 0, z0=10.5, dx=0.5, dz=0.1, dt=1e2,T=1e2, K_aquifer=0.864e-1, K_aquitard=0.864e-5,
        # alpha_L=1e-3, alpha_anisT=1, alpha_anisV=1, porosity=0.3, v_slr=1/1000/365, anis_aquifer=1, sea_level=-0.35, T_init=1e5)
    #create_run_plot_model("Northern_Florida", T_init=2e5, K_aquifer=42.5, z0=310, H=160, h_b=18, x_b=5, D=30, K_aquitard=0.03, L=110000)

if __name__=="__main__":
    main()
