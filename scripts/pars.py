import os
import pickle
import numpy as np
import pdb

class ModelParameters:
    """
        Class to store model parameters

        Inputs:
            name: The name of the model/scenario/realization
            L: Length of the aquifer model [m]
            H: Aquifer Thickness
            D: Aquitard Thickness
            offshore_proportion: portion of the aquifer which is submarine 
            dx: column width
            dz: layer depth
            K_aquifer: Horizontal hydraulic conductivity [m/day]
            K_aquitard: Horizontal hydraulic conductivity [m/day]
            anis_aquifer: ratio between K and vertical hydraulic conductivity
            anis_aquitard: ratio between K and vertical hydraulic conductivity
            sy_aquifer: specific yield 
            sy_aquitard: specific yield 
            ss_aquifer: specific storage [m^-1]
            ss_aquitard: specific storage [m^-1]
            n_aquifer: porosity
            n_aquitard: porosity
            alpha_L: longitudinal dispersivity
            alpha_anisT: ratio between longitudinal and transverse dispersivity
            alpha_anisV: ratio between longitudinal and vertical dispersivity
                !!! I've set this as 0.1, but which werner says is the transverse ratio:
                    but there is no mention of vertical dispersivity !!!
            diff: molecular diffusion coefficient
            dt: timestep [days]
            h_b: inland boundary head with respect to msl
            rho_f: fresh water density [kg/m^3]
            rho_s: salt water density [kg/m^3]
            exe_path: path to the seawat executable
            frequency: frequency of timesteps save for results
    """
    
    def __init__(self, name="none", L=50000, H=16, D=25, z0=76, offshore_proportion=1, 
            sea_level=0, dx=25, dz=1, x_b=4000, K_aquifer=300, K_aquitard=0.01, anis_aquifer=100, anis_aquitard=100,
                sy=0.24, Ss=1e-6,
                n_aquifer=0.3, n_aquitard=0.3, alpha_L=5, alpha_anisT=0.1, alpha_anisV=0.01, 
                diff=8.64e-5, dt=5000, T=5000, T_init=100*10000, v_slr=120/(20000*365), h_b=8.5, rho_f=1000, rho_s=1025, 
                exe_file=r"./exe_path.txt", frequency=1, porosity=0.3, steady=False, T_modern=365*100, h_modern=-1, flux_modern=None, dt_modern=1000, outputs="all", delta_h=None, h_over=None, gradient=1):

        self.name=name
        self.porosity=porosity
        self.L=L
        self.H=H
        self.D=D
        self.z0=z0
        self.x_b=x_b
        self.offshore_proportion=offshore_proportion
        self.sea_level=sea_level
        self.dx=dx
        self.dz=dz
        self.K_aquifer=K_aquifer
        self.K_aquitard=K_aquitard
        self.anis_aquifer=anis_aquifer
        self.anis_aquitard=anis_aquitard
        self.sy=sy
        self.ss=Ss
        self.n_aquifer=porosity
        self.n_aquitard=porosity
        self.alpha_L=alpha_L
        self.alpha_anisT=alpha_anisT
        self.alpha_anisV=alpha_anisV
        self.diff=diff
        self.dt=dt
        self.T=T	
        self.T_init=100*10000
        self.v_slr=v_slr
        self.rho_f=rho_f
        self.rho_s=rho_s
        self.beta = np.arctan((self.z0-self.D-self.H)/self.L)
        self.Lx= self.L+self.x_b
        self.Lz = self.z0+self.x_b*np.tan(self.beta)
        self.exe_path = "/home/superuser/mf6"
        self.frequency=frequency
        self.nrow=1
        self.ncol = int(np.floor(self.Lx/self.dx))
        self.nlay = int(np.floor(self.Lz/self.dz))
        self.steady = steady
        self.T_modern = T_modern
        if h_over == None:
            self.h_b=h_b
        else: 
            self.h_b=self.Lz-self.z0+h_over
        if delta_h == None:
            self.h_modern = h_modern
        else:
            self.h_modern = self.h_b+delta_h
        self.flux_modern = flux_modern
        self.aquifer_cells = None
        self.aquitard_cells = None
        self.top_cells = None
        self.bottom_cells = None
        self.dt_modern = dt_modern
        self.outputs = outputs
        self.gradient = gradient
        self.Ss = Ss
        self.save_parameters()
        # pdb.set_trace()

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d

    def save_parameters(self):
        """
            Save object
        """
        model_ws = f"./model_files/{self.name}"
        if not os.path.exists(model_ws):
            os.makedirs(model_ws)
        model_ws_b = f"./model_files_backup/{self.name}"
        if not os.path.exists(model_ws_b):
            os.makedirs(model_ws_b)

        f = open(f"{model_ws}/pars", 'wb')
        pickle.dump(self, f)
        f.close()
        
        f = open(f"{model_ws_b}/pars", 'wb')
        pickle.dump(self, f)
        f.close()

def load_parameters(name, backup=False):

    if not backup:
        model_ws = f"./model_files/{name}"
    else: 
        model_ws = f"./model_files_backup/{name}"
    	
    f = open(f"{model_ws}/pars", 'rb')
    temp = pickle.load(f)
    f.close()          
    return temp

