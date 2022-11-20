import sloping_confined_aquifer as sca
from pars import ModelParameters, load_parameters
import results

def create_run_plot_model(name, **kwargs):
    """
        Create, run and plot a new scenario
        
        Inputs:
            name: name for the model
            **kwargs: any key word arguments to be passed to the ModelParameters class
        Outputs:
            None
    """

    pars = ModelParameters(name, **kwargs)
    swt = sca.build_steady_model(pars)
    sca.run_model(swt)
    concentration, head, qx, qy, qz = sca.extract_results(name)
    # results.plot_results(name)
    results.plot_gif(name)
    # results.plot_evolutions(name)
    # results.save_metrics(name, fraction=0.05)
    # results.plot_boundary_concentration(name)


def main():
    #create_run_plot_model("canterbury_steady")
    #create_run_plot_model("canterbury_steady", steady=True)
    #create_run_plot_model("canterbury_steady_low", steady=True, sea_level=-20)
    #create_run_plot_model("suriname", L=130000, H=85, z0=185, K_aquifer=25, x_b=5, h_b=3, K_aquitard=0.04, D=30, T_init=2e5, anis_aquifer=100)
    # create_run_plot_model("kooi_case15b", L=400, H=9, D=1, x_b=1, h_b = 0, z0=10.4, dx=0.5, dz=0.1, dt=1e3,T=1e4, K_aquifer=1e-1, K_aquitard=1e-5,
    #     alpha_L=1e-3, alpha_anisT=1, alpha_anisV=1, v_slr=0.1/1000/365, anis_aquifer=1)
    create_run_plot_model("kooi_case17_finer_again", L=400, H=9.5, D=0.6, x_b=1, h_b = 0, z0=10.5, dx=0.25, dz=0.1, dt=1e2,T=1e2, K_aquifer=0.864e-1, K_aquitard=0.864e-5,
        alpha_L=1e-3, alpha_anisT=1, alpha_anisV=1, porosity=0.3, v_slr=1/1000/365, anis_aquifer=1, sea_level=-0.35, T_init=1e5)
    #create_run_plot_model("Northern_Florida", T_init=2e5, K_aquifer=42.5, z0=310, H=160, h_b=18, x_b=5, D=30, K_aquitard=0.03, L=110000)

if __name__=="__main__":
    main()