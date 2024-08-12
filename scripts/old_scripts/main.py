import scripts.old_scripts.sloping_confined_aquifer as sca
from pars import ModelParameters, load_parameters
import results
import mf6
import results_mf6

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
    # # # swt = sca.build_steady_model(pars)
    # # # swt = sca.build_modern_model(name)
    swt = mf6.build_steady_model(pars)
    # print("running")
    mf6.run_model(swt)
    results_mf6.plot_gif(name)
        #results.plot_last_step(name)

    # concentration, head, qx, qy, qz = sca.extract_results(name)
    # results.plot_results(name)
    # results.plot_gif(name)
    # results.metrics(name)
    # results.plot_metrics(name)

    #results.plot_last_step(name+"_modern")
     
def main():
    create_run_plot_model("test", H=10, D=10, K_aquifer=3, K_aquitard=0.001, alpha_L=5, anis_aquifer=10, z0=25, L=10000, h_over=0, delta_h=-5, outputs="all")


if __name__=="__main__":
    main()
