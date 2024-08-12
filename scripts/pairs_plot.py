import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib as mpl
import matplotlib.colors as colors

def fit_polynomial_regression(x1, x2, y, scales, lims):

    # poly = PolynomialFeatures(degree=2)
    poly = PolynomialFeatures(interaction_only=True,include_bias = False)
    if scales[0]=="log":
        x1 = np.log10(x1)
        x1_mesh = np.log10(np.logspace(np.log10(lims[0][0])/np.log10(10),np.log10(lims[0][1]/np.log10(10))))
    else: 
        x1_mesh = np.linspace(lims[0][0],lims[0][1])
    if scales[1]=="log":
        x2 = np.log10(x2)
        x2_mesh = np.log10(np.logspace(np.log10(lims[1][0])/np.log10(10),np.log10(lims[1][1]/np.log10(10))))
    else:
        x2_mesh = np.linspace(lims[1][0],lims[1][1])
    	
    X1, X2 = np.meshgrid(x1_mesh, x2_mesh)
    
    x = np.dstack([x1, x2])[0]
    x_ = poly.fit_transform(x)
    X = np.dstack([X1.flatten(),X2.flatten()])[0]
    X_ = poly.fit_transform(X)
    
    clf = LinearRegression()
    clf.fit(x_, y)
    
    Y = clf.predict(X_)
    if scales[0]=="log":
        X1 = np.power(10, X1)
    if scales[1]=="log":
        X2 = np.power(10, X2)
    
    print(clf.score(x_, y))
    return X1, X2, Y.reshape(50, -1)

def pairs_plot(results, scales, pars, units, lims=None, title=None):
    f, axs = plt.subplots(len(units), len(units), figsize=(6,6), layout="tight")
    #axs[0, 1].axis('off')
    marker_s = 0.1
    for row in range(len(units)):
        for col in range(len(units)):
            if row>col:
                # pdb.set_trace()                           
                if scales[row] == "log":
                    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                # X1, X2, Y = fit_polynomial_regression(results[:, col].T, results[:, row].T, np.log10(results[:, -1].T), [scales[i] for i in [col, row]], [lims[i] for i in [col, row]])
                scatter = axs[row, col].scatter(results[:, col].T, results[:, row].T, c=results[:, -1].T, marker="o", alpha=0.4, s=2, cmap="viridis") 
                # axs[row, col].contour(X1, X2, Y, cmap="viridis") 
                # if col > 0:
                #     axs[row,col].sharey(axs[row][0])
            elif row == col:
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")  
                # axs[row, col].set_yscale("log")           
                axs[row, col].scatter(results[:, col].T, results[:, -1].T, c="grey", marker="o", alpha=0.7, s=2)
                # if row > 0:
                #     axs[row, col].sharey(axs[0][0])
            else:
                axs[row, col].axis('off')

            
            if col == 0 and row == len(units)-1:
                axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
            elif col == 0:
                if row == 0:
                    axs[row, col].set_ylabel(r"$\frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
                    # axs[row, col].set_ylabel(r"$\Delta V^{100}_{ofg}$")

                else:# [{units[row]}]")
                    axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                # plt.setp(axs[row, col].get_xticklabels(), visible=False)
            elif row == len(units)-1:
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
                # axs[row, col].set_yticklabels([])
            else:
                pass
                # axs[row, col].set_xticklabels([])
                # axs[row, col].set_yticklabels([])
            # if col > 0:
            axs[row, col].set_box_aspect(1)  
    cb = plt.colorbar(scatter, ax=axs[1, 2], cmap="viridis", extend="both", location="right")
    cb.ax.set_title(r"$\frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
    cb.ax.yaxis.set_ticks_position("left")
          
   
    plt.show()


def pairs_plot_2(results, scales, pars, units, lims=None, title=None):
    # double pairs plot for multiple factors
    f, axs = plt.subplots(len(units), len(units), figsize=(6,6))
    #axs[0, 1].axis('off')
    marker_s = 0.1
    for row in range(len(units)):
        for col in range(len(units)):
            if row>col:
                # pdb.set_trace()                           
                if scales[row] == "log":
                    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                X1, X2, Y = fit_polynomial_regression(results[:, col].T, results[:, row].T, np.log10(results[:, -2].T), [scales[i] for i in [col, row]], [lims[i] for i in [col, row]])
                scatter = axs[row, col].scatter(results[:, col].T, results[:, row].T, c=np.log10(results[:, -2].T), marker="o", alpha=0.7, s=2) 
                axs[row, col].contour(X1, X2, Y) 
           
            elif row == col:
                axs[row, col].set_yscale("log")
                model1 = LinearRegression()
                model2 = LinearRegression() 
                x = np.linspace(lims[col][0], lims[col][1], 50)
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")  
                    X_ = np.log10(results[:, col].T)
                    x = np.log10(x)

                fit1 = model1.fit(X_.reshape(-1, 1), -np.log10(results[:, -2]).reshape(-1, 1))
                fit2 = model2.fit(X_.reshape(-1, 1), np.log10(results[:, -1]).reshape(-1, 1))
                pred1 = model1.predict(x.T.reshape(-1, 1))
                pred2 = model2.predict(x.T.reshape(-1, 1))
                twin=axs[row, col].twinx()
                axs[row, col].plot(x, pred1, c="blue")
                twin.plot(x, pred2, c="red")
            else:
                # axs[row, col].axis('off')
                # pdb.set_trace()                           
                if scales[row] == "log":
                    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                X1, X2, Y = fit_polynomial_regression(results[:, col].T, results[:, row].T, np.log10(results[:, -1].T), [scales[i] for i in [col, row]], [lims[i] for i in [col, row]])
                scatter = axs[row, col].scatter(results[:, col].T, results[:, row].T, c=np.log10(results[:, -1].T), marker="o", alpha=0.7, s=2, cmap="inferno") 
                axs[row, col].contour(X1, X2, Y, cmap="inferno") 
  

            
            if col == 0 and row == len(units)-1:
                axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
            elif col == 0:
                if row == 0:
                    axs[row, col].set_ylabel(r"$\frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
                else:# [{units[row]}]")
                    axs[row, col].set_ylabel(f"{pars[row]}")# [{units[row]}]")
                plt.setp(axs[row, col].get_xticklabels(), visible=False)
            elif row == len(units)-1:
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
                axs[row, col].set_yticklabels([])
            else:
                axs[row, col].set_xticklabels([])
                axs[row, col].set_yticklabels([])
    
            axs[row, col].set_box_aspect(1)  
    plt.tight_layout() 
    cb = plt.colorbar(scatter, ax=axs[1, 2], cmap="viridis", extend="both", location="right")
    cb.ax.set_title(r"$log_{10} \frac{\Delta V^{100}_{ofg}}{V^{pre}_{ofg}}$")
    cb.ax.yaxis.set_ticks_position("left")
          
   
    plt.show()

def pairs_plot_3(inputs, results, scales, pars, units, lims=None, title=None):
    # pairs plot wihout individual effects
    f, axs = plt.subplots(len(units)-1, len(units)-1, figsize=(9,9), sharex='col', sharey='row', layout="tight")
    plt.rcParams.update({'font.size': 9})

    for row in range(len(units)-1):
        for col in range(len(units)-1):
            if row<col:
                axs[row, col].axis('off')
            else:
                axs[row, col].set_box_aspect(1) 
                if scales[row+1] == "log":
                    axs[row, col].set_yscale("log")
                #axs[row, col].set_ylim(lims[row])
                #axs[row, col].set_xlim(lims[col])
                if scales[col] == "log":
                    axs[row, col].set_xscale("log")
                # else: 
                    #axs[row, col].set_xlim([0, -10])
                X1, X2, Y = fit_polynomial_regression(inputs[:, col].T, inputs[:, row+1].T, results, [scales[i] for i in [col, row+1]], [lims[i] for i in [col, row+1]])
                scatter = axs[row, col].scatter(inputs[:, col].T, inputs[:, row+1].T, c=results, norm=colors.LogNorm(vmin=0.01), marker="o", alpha=1, s=1) 
                axs[row, col].contour(X1, X2, Y) 
            if col == 0 and row == len(units)-2:
                axs[row, col].set_ylabel(f"{pars[row+1]}")# [{units[row]}]")
                axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
            elif col == 0:
                axs[row, col].set_ylabel(f"{pars[row+1]}")# [{units[row]}]")
            elif row == len(units)-2:
                 axs[row, col].set_xlabel(f"{pars[col]}")# [{units[col]}]")
                 # axs[row, col].hlines(0.75, lims[col][0], lims[col][1], linestyles="dashed", color="red", linewidth=0.5)


    cb = plt.colorbar(scatter, ax=axs[0, 1:], cmap="viridis", location="top", shrink=0.75, ticklocation="bottom")
    # cb.ax.set_xticklabels(['1', '0.5', '0.1'])
    cb.ax.set_title("Fraction OFG Salinized")
    plt.show()
