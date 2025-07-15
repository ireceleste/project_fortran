import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
import os


#########################################################
#  Function to plot the results of the fit
#########################################################

def plotFit(results, results_path='output/'):

    
    if not os.path.exists(results_path):
        os.makedirs(results_path)


    data = np.loadtxt(results, comments='#')
    n_columns = data.shape[1]


    x_raw = data[1:, 0]
    y_hist_raw = data[1:, 1]
    err_y_hist_raw = data[1:, 2]
    y_fcn_raw = data[1:, 3]

    min_fcn = data[0, 0]
    ndof = data[0, 1]

    pars = []
    err_pars = []


    # Fix input data for plotting reasons

    # sort on x
    indexes = np.argsort(x_raw)
    x = x_raw[indexes]
    y_hist = y_hist_raw[indexes]
    err_y_hist= err_y_hist_raw[indexes]
    y_fcn = y_fcn_raw[indexes]

    # create a smooth version of the model y_fcn using an spline
    x_dense = np.linspace(x.min(), x.max(), 1000)
    spl = make_interp_spline(x, y_fcn, k=3)
    y_fcn_dense = spl(x_dense)



    # Extract parameters and their errors
    for i in range(2, n_columns, 2):
        pars.append(data[0, i])
        err_pars.append(data[0, i+1])
    pars_names = [f'$p_{i-1}$' for i in range(2, n_columns)]

    

    # Compute pulls
    pulls = np.zeros_like(y_hist)
    mask = err_y_hist != 0
    pulls[mask] = (y_hist[mask] - y_fcn[mask]) / err_y_hist[mask]


    # Create figure with two subplots: top for fit, bottom for pulls

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6), gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0.03)

    # --- Plot histogram and fitted function ---

    ax1.plot(x_dense, y_fcn_dense, label='Fitted pdf', color='red')
    ax1.errorbar(x, y_hist, yerr=err_y_hist, fmt='.', label='Data', capsize=2, color='blue')

    fit_info = (
        f"$\chi^2$/ndof = {min_fcn:.1f} / {ndof:.0f}\n"
    )
    for i in range(len(pars)):
        if err_pars[i] > 1:
            fit_info  += f"{pars_names[i]} = {pars[i]:.0f} $\pm$ {err_pars[i]:.0f}\n"
        elif err_pars[i] > 0.1:
            fit_info += f"{pars_names[i]} = {pars[i]:.1f} $\pm$ {err_pars[i]:.1f}\n"
        elif err_pars[i] > 0.01:
            fit_info += f"{pars_names[i]} = {pars[i]:.2f} $\pm$ {err_pars[i]:.2f}\n"
        elif err_pars[i] > 0.001:
            fit_info += f"{pars_names[i]} = {pars[i]:.3f} $\pm$ {err_pars[i]:.3f}\n"
        else:
            fit_info += f"{pars_names[i]} = {pars[i]:.4f} $\pm$ {err_pars[i]:.4f}\n"
    ax1.plot([], [], ' ', label=fit_info)  


    ax1.legend(loc='best', fontsize=11)
    ax1.set_ylabel("y")
    ax1.grid(True)


    # --- Plot pulls ---

    ax2.axhline(0, color='black', linewidth=1, linestyle='--')
    ax2.errorbar(x, pulls, yerr=1, fmt='.',  label='Pulls', capsize=2, color = 'blue')
    ax2.set_ylabel('Pulls')
    ax2.grid(True)
    ax2.set_xlabel("x")

    plt.tight_layout()
    plt.savefig(results_path + "fit_results.png")

    print(f"Plot saved to {results_path}fit_results.png")


#########################################################
#  Main function to run the script
#########################################################

if __name__ == "__main__":

    import sys
    if len(sys.argv) > 2:
        raise ValueError("Usage: python plot.py <path_to_results.txt>")
    if(len(sys.argv) == 1):
        results = "output/fit_results.txt"
    else:
        results = sys.argv[1]
    plotFit(results, results_path='output/')

