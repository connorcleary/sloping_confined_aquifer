from pars import ModelParameters, load_parameters
import mf6
import flopy
import matplotlib.pyplot as plt
import pdb
import numpy as np

def main():
   pars = ModelParameters("to_plot", L=10000, x_b=3000, z0=43.33, H=20, D=20)
   swt = mf6.build_steady_model(pars)
   f,ax = plt.subplots()
   #pdb.set_trace()
   cs = flopy.plot.PlotCrossSection(model=swt.to_plot,ax= ax, line={"row":0})
   cs.plot_grid(lw=0.1, color="0.5")
   cs.plot_ibound()
   x = [0, 3000, 8000,13000]
   xticks =[int(-3000), int(0), int(5000),int(10000)]
   plt.xticks(x, xticks)
   ax.set_xlabel("Distance offshore [m]")
   ax.set_ylabel("Depth [m]")
   plt.show()
   
if __name__=="__main__":
    main()
