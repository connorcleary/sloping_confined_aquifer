import numpy as np

def raleigh():
    return 

def peclet(dz, mu, n, K, D=8.64e-5):
     return dz*(25*9.81*k(K)/(mu*n))/D

def vertical_courant(K, n, dz):
    return dz/(vz(K, n))

def k(K): # permeability 
    return K*0.001/1000/9.81

def vz(K, n): # velocity of descending fingers
    return 25*9.81*k(K)/(4*0.001*n)

def v_adv(K, n, z, L):
    rho_f = 1000
    rho_s = 1025
    beta = np.arctan(z/L)
    return K*rho_s/rho_f/n*np.tan(beta)

def vertical_peclet(K, n, D=8.64e-5): 
    return 2*D/abs(vz(K, n)*n)

def vertical_rayleigh(K, Ra_crit=7, D=8.64e-5):
    return Ra_crit*0.001*D/(25*9.81*k(K))

def vertical_rayleigh2(K, n, alpha_l=1, Ra_crit=7, D=8.64e-5):
    return Ra_crit*0.001*(D+vz(K,n)*alpha_l)/(alpha_l*25*9.81*k(K))

def horizontal_rayleigh(K, Ra_crit=7, D=8.64e-5):
    return 2*Ra_crit*0.001*D/(25*9.81*k(K))

def horizontal_rayleigh2(K, n, alpha_l=1, Ra_crit=7, D=8.64e-5):
    return 2*Ra_crit*0.001*(D+vz(K,n)*alpha_l)/(alpha_l*25*9.81*k(K))

def horizontal_courant(K, n, z, L, dx):
    return dx/(v_adv(K, n, z, L))

def horizontal_peclet(K, n,  z, L, D=8.64e-5):
    # return np.inf
    return 2*D/abs(v_adv(K, n, z, L))

def main(Kz, Kx, n, alpha_l=0, z=35, L=50000):
    if alpha_l == 0:
        dz = min(vertical_peclet(Kz, n), vertical_rayleigh(Kz))
        dx = min(horizontal_rayleigh(Kz, n), horizontal_peclet(Kx, n,  z, L))
    else: 
        dz = min(vertical_peclet(Kz, n), vertical_rayleigh2(Kz, n, alpha_l))
        dx = min(horizontal_rayleigh2(Kz, n, alpha_l), horizontal_peclet(Kx, n,  z, L))

    dt = min(vertical_courant(Kz, n, dz), horizontal_courant(Kx, n, z, L, dx))

    print(f"dz = {dz}m, dx = {dx}m, dt = {dt}days")

if __name__=="__main__":
    print(f"dt= {v_adv(30, 0.3, 3, 10000)}")
    print(f"dt= {horizontal_courant(30, 0.3, 3, 10000, horizontal_peclet(30, 0.3, 44, 10000))}")
    print(f"dx= {horizontal_peclet(30, 0.3, 3, 10000)}")


