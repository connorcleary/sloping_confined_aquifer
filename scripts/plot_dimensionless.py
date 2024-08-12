from results_smaller import load_smaller
from pars import ModelParameters, load_parameters
import mf6
import numpy as np
import matplotlib.pyplot as plt

def raleigh(qx, pars, pars2):

    qxs = [qx[cell[0]][cell[1]][cell[2]] for cell in pars2.top_cells[3000:]]
    return [np.abs(pars.porosity*0.025*pars.alpha_L*pars.K_aquifer/pars.anis_aquifer/(1*pars.diff+np.abs(pars.alpha_L*qxi))) for qxi in qxs]

def back_dispersion(qz, pars, pars2):
    qzs = [qz[int(cell[0]-40)][cell[1]][cell[2]] for cell in pars2.top_cells[3000:]]
    return [(pars.alpha_L*pars.alpha_anisV*qzi)/pars.porosity/pars.diff for qzi in qzs]

def peclet(q, pars, direction):
    if direction == 0: # horizontal
        return np.abs(q)/pars.porosity*pars.dx/(pars.diff+pars.alpha_L*np.abs(q)/pars.porosity)
    else: # veritcal
        return np.abs(q)/pars.porosity*pars.dz/(pars.diff+pars.alpha_L*pars.alpha_anisV*np.abs(q)/pars.porosity)
    
def plot_peclet(qx, qz, pars):
    
    Pe = np.zeros_like(qx)
    # for each row
    for col in range(pars.ncol):
        # for each row
        for lay in range(pars.nlay):
            pe = 0
            # calculate if PeH > 2
            if peclet(qx[lay, 0, col], pars, 0) > 2:
                if peclet(qz[lay, 0, col], pars, 1) > 2:
                    pe = 3
                else: 
                    pe = 1
            else:
                if peclet(qz[lay, 0, col], pars, 1) > 2:
                    pe = 2
                else: 
                    pe = 0
            Pe[lay, 0 , col] = pe
            # calculate if PeV > 2
    pass

def main():

    ensemble_name="base30"
    Kas = [1, 10]
    KbVs = [1e-4, 1e-3]
    alpha = 10
    anis = 100
    n = 0.5
    h_modern = 1
    epsilons = []
    raleighs = []

    pars2 = ModelParameters("...", H=20, D=20, K_aquifer=1, 
                        K_aquitard=1, alpha_L=10, anis_aquifer=100, anis_aquitard=100,
                        z0=50, L=30000, x_b=3000,  
                        h_modern=0, porosity=0.5, gradient=1, Ss=1e-6, outputs="all", dx=1, dz=0.5, dt=100) 
       
    pars2 = mf6.build_steady_model(pars2, True)

    for Ka in Kas: 
        for KbV in KbVs:
            name = f"{Ka}_{KbV}_{h_modern}"
            results = load_smaller(name)
            plot_peclet(results["qx"][1], results["qz"][1], load_parameters(name, backup = True))
            # raleighi = raleigh(results["qx"][1], load_parameters(name, backup = True), pars2)
            # raleighs.append(raleighi)
            # epsilon = back_dispersion(results["qz"][1], load_parameters(name, backup = True), pars2)
            # epsilons.append(epsilon)

    f, axs = plt.subplots(2, 1, figsize=(4, 4))
    epsilons = np.array(epsilons)
    raleighs = np.array(raleighs)
    # eps_ave = np.mean(epsilons.reshape(4, 300, -1), axis=2)
    rals_ave = np.mean(raleighs.reshape(4, 300, -1), axis=2)
    axs[0].plot(epsilons.T)
    axs[1].plot(np.linspace(0, 30000, 300), rals_ave.T)
    # axs[0].set_yscale("log")
    axs[1].set_yscale("log")
    
    plt.show()
    pass

if __name__ == "__main__":
    main()