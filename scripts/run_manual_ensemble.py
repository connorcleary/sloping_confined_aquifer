import mf6
from pars import ModelParameters, load_parameters
from multiprocessing import Pool

def setup_pars():
    ModelParameters("a0")
    ModelParameters("a1", K_aquifer=1)
    ModelParameters("a2", K_aquifer=100)
    ModelParameters("a3", K_aquitard=0.001)
    ModelParameters("a4", K_aquitard=0.1)
    ModelParameters("a5", anis_aquifer=10)
    ModelParameters("a6", anis_aquifer=1000)    
    ModelParameters("a7", anis_aquitard=10)
    ModelParameters("a8", anis_aquitard=1000)    
    ModelParameters("a9", h_modern=-2)
    ModelParameters("a10", h_modern=-5)
    ModelParameters("a11", H=36, D=5)
    ModelParameters("a12", H=21, D=20)
    
    return [f"a{i}" for i in range(13)]

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

