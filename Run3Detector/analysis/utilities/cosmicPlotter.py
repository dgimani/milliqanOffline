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
    p0 = [max(values[fit_mask]),60,20]
    popt, pcov = curve_fit(func, bins[fit_mask], values[fit_mask], p0=p0, sigma=var[fit_mask])
    x_fit=np.linspace(fit_low, fit_high, 1000)
    y_fit=gaussian(x_fit,*popt)
    return popt, pcov, x_fit, y_fit

def straightLinePath(pulses):  #NOTE: as of now this function returns 4 hits in a row but includes all after pulses and also includes multiple paths if any exist
	local_index = ak.local_index(pulses)

	l0_mask = pulses.layer == 0
	l1_mask = pulses.layer == 1
	l2_mask = pulses.layer == 2
	l3_mask = pulses.layer == 3

	l0_local_indices = local_index[l0_mask]
	l1_local_indices = local_index[l1_mask]
	l2_local_indices = local_index[l2_mask]
	l3_local_indices = local_index[l3_mask]

	l0_pulses = pulses[l0_mask]
	l1_pulses = pulses[l1_mask]
	l2_pulses = pulses[l2_mask]
	l3_pulses = pulses[l3_mask]

	l0l1_pairs = ak.cartesian([l0_pulses,l1_pulses],axis=1)
	l0l2_pairs = ak.cartesian([l0_pulses,l2_pulses],axis=1)
	l0l3_pairs = ak.cartesian([l0_pulses,l3_pulses],axis=1)
	l1l2_pairs = ak.cartesian([l1_pulses,l2_pulses],axis=1)
	l1l3_pairs = ak.cartesian([l1_pulses,l3_pulses],axis=1)
	l2l3_pairs = ak.cartesian([l2_pulses,l3_pulses],axis=1)

	l0l1_argpairs = ak.argcartesian([l0_pulses,l1_pulses],axis=1)
	l0l2_argpairs = ak.argcartesian([l0_pulses,l2_pulses],axis=1)
	l0l3_argpairs = ak.argcartesian([l0_pulses,l3_pulses],axis=1)
	l1l2_argpairs = ak.argcartesian([l1_pulses,l2_pulses],axis=1)
	l1l3_argpairs = ak.argcartesian([l1_pulses,l3_pulses],axis=1)
	l2l3_argpairs = ak.argcartesian([l2_pulses,l3_pulses],axis=1)

	l0l1_l0half, l0l1_l1half = ak.unzip(l0l1_pairs)
	l0l2_l0half, l0l2_l2half = ak.unzip(l0l2_pairs)
	l0l3_l0half, l0l3_l3half = ak.unzip(l0l3_pairs)
	l1l2_l1half, l1l2_l2half = ak.unzip(l1l2_pairs)
	l1l3_l1half, l1l3_l3half = ak.unzip(l1l3_pairs)
	l2l3_l2half, l2l3_l3half = ak.unzip(l2l3_pairs)

	l0l1_argl0half, l0l1_argl1half = ak.unzip(l0l1_argpairs)
	l0l2_argl0half, l0l2_argl2half = ak.unzip(l0l2_argpairs)
	l0l3_argl0half, l0l3_argl3half = ak.unzip(l0l3_argpairs)
	l1l2_argl1half, l1l2_argl2half = ak.unzip(l1l2_argpairs)
	l1l3_argl1half, l1l3_argl3half = ak.unzip(l1l3_argpairs)
	l2l3_argl2half, l2l3_argl3half = ak.unzip(l2l3_argpairs)

	l0l1_timeMask = np.abs(l0l1_l0half.time - l0l1_l1half.time) < 400
	l0l2_timeMask = np.abs(l0l2_l0half.time - l0l2_l2half.time) < 400
	l0l3_timeMask = np.abs(l0l3_l0half.time - l0l3_l3half.time) < 400
	l1l2_timeMask = np.abs(l1l2_l1half.time - l1l2_l2half.time) < 400
	l1l3_timeMask = np.abs(l1l3_l1half.time - l1l3_l3half.time) < 400
	l2l3_timeMask = np.abs(l2l3_l2half.time - l2l3_l3half.time) < 400

	l0l1_rowMask = l0l1_l0half.row == l0l1_l1half.row
	l0l2_rowMask = l0l2_l0half.row == l0l2_l2half.row
	l0l3_rowMask = l0l3_l0half.row == l0l3_l3half.row
	l1l2_rowMask = l1l2_l1half.row == l1l2_l2half.row
	l1l3_rowMask = l1l3_l1half.row == l1l3_l3half.row
	l2l3_rowMask = l2l3_l2half.row == l2l3_l3half.row

	l0l1_colMask = l0l1_l0half.column == l0l1_l1half.column
	l0l2_colMask = l0l2_l0half.column == l0l2_l2half.column
	l0l3_colMask = l0l3_l0half.column == l0l3_l3half.column
	l1l2_colMask = l1l2_l1half.column == l1l2_l2half.column
	l1l3_colMask = l1l3_l1half.column == l1l3_l3half.column
	l2l3_colMask = l2l3_l2half.column == l2l3_l3half.column

	l0l1pointingMask = l0l1_colMask & l0l1_rowMask & l0l1_timeMask
	l0l2pointingMask = l0l2_colMask & l0l2_rowMask & l0l2_timeMask
	l0l3pointingMask = l0l3_colMask & l0l3_rowMask & l0l3_timeMask
	l1l2pointingMask = l1l2_colMask & l1l2_rowMask & l1l2_timeMask
	l1l3pointingMask = l1l3_colMask & l1l3_rowMask & l1l3_timeMask
	l2l3pointingMask = l2l3_colMask & l2l3_rowMask & l2l3_timeMask

	l0_common_indices = ak.Array([
	    np.intersect1d(np.intersect1d(l0l1, l0l2), l0l3)
	    for l0l1, l0l2, l0l3 in zip(
		ak.to_list(l0l1_argl0half[l0l1pointingMask]),
		ak.to_list(l0l2_argl0half[l0l2pointingMask]),
		ak.to_list(l0l3_argl0half[l0l3pointingMask])
	    )
	])

	l1_common_indices = ak.Array([
	    np.intersect1d(np.intersect1d(l0l1, l1l2), l1l3)
	    for l0l1, l1l2, l1l3 in zip(
		ak.to_list(l0l1_argl1half[l0l1pointingMask]),
		ak.to_list(l1l2_argl1half[l1l2pointingMask]),
		ak.to_list(l1l3_argl1half[l1l3pointingMask])
	    )
	])

	l2_common_indices = ak.Array([
	    np.intersect1d(np.intersect1d(l0l2, l1l2), l2l3)
	    for l0l2, l1l2, l2l3 in zip(
		ak.to_list(l0l2_argl2half[l0l2pointingMask]),
		ak.to_list(l1l2_argl2half[l1l2pointingMask]),
		ak.to_list(l2l3_argl2half[l2l3pointingMask])
	    )
	])

	l3_common_indices = ak.Array([
	    np.intersect1d(np.intersect1d(l0l3, l1l3), l2l3)
	    for l0l3, l1l3, l2l3 in zip(
		ak.to_list(l0l3_argl3half[l0l3pointingMask]),
		ak.to_list(l1l3_argl3half[l1l3pointingMask]),
		ak.to_list(l2l3_argl3half[l2l3pointingMask])
	    )
	])

	l0_pointingPulses = l0_pulses[l0_common_indices]
	l1_pointingPulses = l1_pulses[l1_common_indices]
	l2_pointingPulses = l2_pulses[l2_common_indices]
	l3_pointingPulses = l3_pulses[l3_common_indices]

	l0_pointingIdx = l0_local_indices[l0_common_indices]
	l1_pointingIdx = l1_local_indices[l1_common_indices]
	l2_pointingIdx = l2_local_indices[l2_common_indices]
	l3_pointingIdx = l3_local_indices[l3_common_indices]

	straightLineMask = ak.concatenate([l0_pointingIdx,l1_pointingIdx,l2_pointingIdx,l3_pointingIdx],axis=1)
	return pulses[straightLineMask]


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
    branches_needed.append("height")
    branches_needed.append("area")
    branches_needed.append("layer")
    branches_needed.append("row")
    branches_needed.append("column")
    branches_needed.append("time")
    branches_needed.append("chan")
    branches_needed.append("ipulse")

    # build the histogram
    h = hist.Hist(
        hist.axis.Regular(80, rlow, rhigh, label=plot_var)
    )

    file_pattern = f"/net/cms18/cms18r0/milliqan/Run3Offline/v36/slab/MilliQanSlab_Run{run}.*_v36.root:t"
    for branches in uproot.iterate(file_pattern, branches_needed, step_size=1000):

        pulses = ak.zip(
			{
				"height": branches["height"],
				"area": branches["area"],
				"row": branches["row"],
				"column": branches["column"],
				"layer": branches["layer"],
				"chan": branches["chan"],
				"time": branches["time"],
				"ipulse": branches["ipulse"],
			}
	)

	#pairs = ak.combinations(pulses, 2)
	#pulse1, pulse2 = ak.unzip(pairs)

        straightLinePulses = straightLinePath(pulses)

        # Start with a mask of “select everything”
        cut_mask = ak.ones_like(straightLinePulses[plot_var], dtype=bool)

        # apply each cut
        for (cutvar, cut_low, cut_high) in cut_specs:
            cut_mask = cut_mask & ((straightLinePulses[cutvar] >= cut_low) & (straightLinePulses[cutvar] <= cut_high))

        # fill histogram with masked values
        arr_plot = straightLinePulses[plot_var][cut_mask]
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
    plt.text(0.73,0.93,f"Run {run}", transform=plt.gca().transAxes)
    if 'chan' in locals(): plt.text(0.73,0.88,f"Channel {chan}",transform=plt.gca().transAxes)
    if do_fit:
        plt.text(0.73,0.83,rf'$\mu$={round(popt[1],1)}',transform=plt.gca().transAxes)
        plt.text(0.73,0.78,rf'$\sigma$={round(popt[2],1)}',transform=plt.gca().transAxes)
    plt.xlabel(plot_var,loc='center')
    #plt.xlabel("Pulse Height [mV]",loc='center')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    make_plot(args.run, args.plot_var, args.plot_range, args.cuts, args.do_fit, args.fit_range)

