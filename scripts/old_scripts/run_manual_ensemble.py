import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool

def setup_pars():
    # ModelParameters("d0")
    #ModelParameters("d1", K_aquifer=1)
    #ModelParameters("d2", K_aquifer=100)
    # ModelParameters("c3", K_aquitard=0.001)
    # ModelParameters("c4", K_aquitard=0.1)
    # ModelParameters("c5", anis_aquifer=10)
    # ModelParameters("c6", anis_aquifer=1000)    
    # ModelParameters("c7", anis_aquitard=10)
    # ModelParameters("c8", anis_aquitard=1000)    
    #ModelParameters("d9", h_modern=-2)
    # ModelParameters("d10", h_modern=-5)
    # ModelParameters("c11", H=36, D=5)
    # ModelParameters("c12", H=21, D=20)
    
    # ModelParameters("l0", L=20000, z0=47, h_modern=-2)
    # ModelParameters("l1", L=20000, z0=47, h_modern=-2, K_aquifer=5)
    # ModelParameters("l2", L=20000, z0=47, h_modern=-2, K_aquifer=20)
    
    ModelParameters("test", L=20000, z0=47, h_modern=-2, K_aquifer=20)
    return [f"test{i}" for i in range(1)]

def run_model(name):
    pars = load_parameters(name)
    sim = mf6.build_steady_model(pars)
    try:
        mf6.run_model(sim)
        print(f"{name} success")
    except:
        print(f"{name} failure")

def run_ensemble(names):
    p = Pool(processes=8)
    p.map(run_model, names)

def main():
    names = setup_pars()
    run_ensemble(names)

if __name__=="__main__":
    main()

