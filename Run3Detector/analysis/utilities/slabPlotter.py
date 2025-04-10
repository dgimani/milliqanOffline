import uproot
import awkward as ak
import numpy as np
import argparse
import hist
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mplhep as hep

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=int)
    parser.add_argument("--plot_var", required=True)
    parser.add_argument("--plot_range", required=True, help="e.g. 0-100")
    # Allow arbitrary many cuts in the format var:low-high
    parser.add_argument("--cuts", nargs="*", default=[],
                        help='List of cuts, e.g. --cuts "height:10-50" "time:100-200"')
    parser.add_argument("--do_fit", default=False, help="Whether to perform a gaussian fit")
    parser.add_argument("--fit_range", help="e.g. 50-70")
    return parser.parse_args()

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-0.5 * ((x - mean)/sigma)**2)

def fit(func, hist, fit_range):
    fit_low, fit_high = (float(x) for x in fit_range.split("_"))
    bins = hist.axes[0].centers
    values = hist.values()
    var = hist.variances()
    fit_mask = (bins >= fit_low) & (bins <= fit_high)
    p0 = [max(values[fit_mask]),(fit_high+fit_low)/2,(fit_high-fit_low)]
    popt, pcov = curve_fit(func, bins[fit_mask], values[fit_mask], p0=p0, sigma=var[fit_mask])
    x_fit=np.linspace(fit_low, fit_high, 1000)
    y_fit=gaussian(x_fit,*popt)
    return popt, pcov, x_fit, y_fit


def make_plot(run, plot_var, plot_range, cuts, do_fit=False, fit_range=None, fit_func=gaussian):
    # parse the plot range
    rlow, rhigh = (float(x) for x in plot_range.split("_"))

    # parse the cut strings: each is like "var:low-high"
    cut_specs = []
    for c in cuts:
        var_str, range_str = c.split(":")
        cut_low, cut_high = (float(x) for x in range_str.split("_"))
        if var_str == "chan": chan = int(cut_low)
        cut_specs.append((var_str, cut_low, cut_high))

    # gather all needed branches: the plot_var plus any cut variables
    branches_needed = [plot_var]
    for (cutvar, _, _) in cut_specs:
        if cutvar not in branches_needed:
            branches_needed.append(cutvar)

    # build the histogram
    h = hist.Hist(
        hist.axis.Regular(80, rlow, rhigh, label=plot_var)
    )

    file_pattern = f"/net/cms18/cms18r0/milliqan/Run3Offline/v36/slab/MilliQanSlab_Run{run}.*_v36.root:t"
    for branches in uproot.iterate(file_pattern, branches_needed, step_size=1000):

        # Start with a mask of “select everything”
        cut_mask = ak.ones_like(branches[plot_var], dtype=bool)

        # apply each cut
        for (cutvar, cut_low, cut_high) in cut_specs:
            cut_mask = cut_mask & ((branches[cutvar] >= cut_low) & (branches[cutvar] <= cut_high))

        # fill histogram with masked values
        arr_plot = branches[plot_var][cut_mask]
        h.fill(ak.ravel(arr_plot)) 


    # draw the histogram
    hep.style.use("CMS")
    h.plot(yerr=False)
    if do_fit:
        popt, pcov, xfit, yfit = fit(fit_func, h, fit_range)
        plt.plot(xfit,yfit)
    plt.title(f"milliQan Preliminary",loc="left")
    plt.title(f"Slab Detector",loc="right")
    plt.grid(True)
    plt.yscale('log')
    plt.ylabel('Number of Pulses')
    plt.ylim(10,1000)
    plt.text(0.73,0.93,f"Run {run}", transform=plt.gca().transAxes)
    if 'chan' in locals(): plt.text(0.73,0.88,f"Channel {chan}",transform=plt.gca().transAxes)
    if do_fit:
        plt.text(0.73,0.83,rf'$\mu$={int(popt[1])}$\pm${int(np.sqrt(pcov[1][1]))}',transform=plt.gca().transAxes)
        plt.text(0.73,0.78,rf'$\sigma$={int(popt[2])}$\pm${int(np.sqrt(pcov[2][2]))}',transform=plt.gca().transAxes)
    if(plot_var == "height"): plt.xlabel("Pulse Height [mV]",loc='center')
    if(plot_var == "area"): plt.xlabel("Pulse Area [nVs]",loc='center')
    else: plt.xlabel(plot_var,loc='center')
    #plt.xlabel("Pulse Height [mV]",loc='center')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    make_plot(args.run, args.plot_var, args.plot_range, args.cuts, args.do_fit, args.fit_range)

