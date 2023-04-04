import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool

def setup_pars():
    ModelParameters("c0")
    ModelParameters("c1", K_aquifer=1)
    ModelParameters("c2", K_aquifer=100)
    ModelParameters("c3", K_aquitard=0.001)
    ModelParameters("c4", K_aquitard=0.1)
    ModelParameters("c5", anis_aquifer=10)
    ModelParameters("c6", anis_aquifer=1000)    
    ModelParameters("c7", anis_aquitard=10)
    ModelParameters("c8", anis_aquitard=1000)    
    ModelParameters("c9", h_modern=-2)
    ModelParameters("c10", h_modern=-5)
    ModelParameters("c11", H=36, D=5)
    ModelParameters("c12", H=21, D=20)
    
    return [f"c{i}" for i in range(13)]

def run_model(name):
    pars = load_parameters(name)
    sim = mf6.build_steady_model(pars)
    try:
        mf6.run_model(sim)
        print(f"{name} success")
    except:
        print(f"{name} failure")

def run_ensemble(names):
    p = Pool(processes=16)
    p.map(run_model, names)

def main():
    names = setup_pars()
    run_ensemble(names)

if __name__=="__main__":
    main()

