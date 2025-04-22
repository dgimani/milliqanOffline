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
    parser.add_argument("--runs", required=True, nargs="*", default=[], type=int)
    parser.add_argument("--plot_var", help="A variable you wish to plot")
    parser.add_argument("--layer", type=int, help="The layer of the detector (other than 0)")
    parser.add_argument("--row", type=int, help="The row of the detector")
    parser.add_argument("--col", type=int, help="The column of the detector")
    parser.add_argument("--chan_type", type=int, help="0=even,1=odd channels")
    parser.add_argument("--plot_range", required=True, help="e.g. 0-100")
    # Allow arbitrary many cuts in the format var:low-high
    parser.add_argument("--cuts", nargs="*", default=[],
                        help='List of cuts, e.g. --cuts "height:10-50" "time:100-200"')
    parser.add_argument("--do_fit", default=False, help="Whether to perform a gaussian fit")
    parser.add_argument("--fit_range", help="e.g. 50-70")
    parser.add_argument("--save_all", default=False, help="Whether to save all figures")
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

def slopes(pulse1, pulse2):
	#By convention pulse 1 should always be of a higher layer than pulse 2
	#All lengths in cm
	#Slabs are 40x60cm, ~65cm vertical between layers, ~95cm offset between fronts of two consecutive layers
	inv_z_slope = (60*(pulse1.row - pulse2.row) + 95*(pulse1.layer-pulse2.layer)) / (65*(pulse1.layer - pulse2.layer))
	inv_y_slope = (40*(pulse1.column - pulse2.column)) / (65 * (pulse1.layer - pulse2.layer)) 
	return inv_z_slope, inv_y_slope #Returns 1/slope


def straightLinePath_with_dt(pulses,chan_type):
    local_index = ak.local_index(pulses)
    # Separate out each layer’s pulses
    l0_mask = (pulses.layer == 0)
    l1_mask = (pulses.layer == 1)
    l2_mask = (pulses.layer == 2)
    l3_mask = (pulses.layer == 3)

    l0_pulses, l0_local_idx = pulses[l0_mask], local_index[l0_mask]
    l1_pulses, l1_local_idx = pulses[l1_mask], local_index[l1_mask]
    l2_pulses, l2_local_idx = pulses[l2_mask], local_index[l2_mask]
    l3_pulses, l3_local_idx = pulses[l3_mask], local_index[l3_mask]
    
    # Build the 4D Cartesian product over layers 0,1,2,3
    # Each event now has an array of quadruples (i0, i1, i2, i3) 
    # that reference l0_pulses[i0], l1_pulses[i1], ...
    quadruples = ak.argcartesian([l0_pulses, l1_pulses, l2_pulses, l3_pulses], axis=1)
    i0, i1, i2, i3 = ak.unzip(quadruples)

    # Pull out the “candidate” pulses in each layer
    p0 = l0_pulses[i0]
    p1 = l1_pulses[i1]
    p2 = l2_pulses[i2]
    p3 = l3_pulses[i3]

    # Now define a “straight line” mask. For a simple “exact match”
    # in row & column, do:
    row_match = (p0.row == p1.row) & (p1.row == p2.row) & (p2.row == p3.row)
    col_match = (p0.column == p1.column) & (p1.column == p2.column) & (p2.column == p3.column)
    chan_match = ((p0.chan%2) == chan_type) & ((p1.chan%2) == chan_type) & ((p2.chan%2) == chan_type) & ((p3.chan%2) == chan_type)

    mask = row_match & col_match & chan_match
    
    # Keep only the quadruples that pass
    p0 = p0[mask]
    p1 = p1[mask]
    p2 = p2[mask]
    p3 = p3[mask]
    
    # If you also want an index mask back to the original pulses array:
    # By construction, i0, i1, i2, i3 are indices in each sub-array, so
    # you can zip them back into a single array of shape (N,4).
    # (Then pulses[straightLineMask] would pick out just those hits.)
    i0f = i0[mask]
    i1f = i1[mask]
    i2f = i2[mask]
    i3f = i3[mask]

    g0, g1, g2, g3 = l0_local_idx[i0f], l1_local_idx[i1f], l2_local_idx[i2f], l3_local_idx[i3f],
    all_idx  = ak.concatenate([g0, g1, g2, g3], axis=1)
    uniq_idx = ak.Array([np.unique(evt) for evt in ak.to_list(all_idx)])
    
    straightLinePulses = pulses[uniq_idx]

    return uniq_idx, straightLinePulses

def add_all(straightLinePulses, chanmask, l0mask, dt_hists, time_2Dhists, dt_2Dhists):
	for i, row_ in enumerate(range(0,4)):
		for j, col_ in enumerate(range(0,3)):
			rowmask = straightLinePulses.row == row_
			colmask = straightLinePulses.column == col_
			dts = []
			for k, layer_ in enumerate(range(1,4)):
				layermask = straightLinePulses.layer == layer_
				p0 = straightLinePulses[rowmask & colmask & chanmask & l0mask]
				pL = straightLinePulses[rowmask & colmask & chanmask & layermask]
				dt = p0.timeFit - pL.timeFit
				#dt = ak.firsts(p0.timeFit) - ak.firsts(pL.timeFit)
				#dt = ak.drop_none(dt)
				print(ak.ravel(dt))
				dts.append(dt)
				idx = k + 3*j + 9*i 
				dt_hists[idx].fill(ak.ravel(dt))
				both = len(ak.drop_none(ak.firsts(p0.timeFit))) == len(ak.drop_none(ak.firsts(pL.timeFit)))
				#print(both)
				time_2Dhists[idx].fill(ak.ravel(ak.drop_none(ak.firsts(p0.timeFit))), ak.ravel(ak.drop_none(ak.firsts(pL.timeFit))))
				#if(both): time_2Dhists[idx].fill(ak.ravel(ak.drop_none(ak.firsts(p0.timeFit))), ak.ravel(ak.drop_none(ak.firsts(pL.timeFit))))
			dt_2Dhists[idx-2].fill(ak.ravel(dts[0]), ak.ravel(dts[1]))
			dt_2Dhists[idx-1].fill(ak.ravel(dts[0]), ak.ravel(dts[2]))
			dt_2Dhists[idx].fill(ak.ravel(dts[1]),ak.ravel(dts[2]))
			#if(ak.any(dts[0]) & ak.any(dts[1]) & (len(dts[0]) == len(dts[1])) ): dt_2Dhists[idx-2].fill(ak.ravel(dts[0]), ak.ravel(dts[1]))
			#if(ak.any(dts[0]) & ak.any(dts[2]) & (len(dts[0]) == len(dts[2])) ): dt_2Dhists[idx-1].fill(ak.ravel(dts[0]), ak.ravel(dts[2]))
			#if(ak.any(dts[1]) & ak.any(dts[2]) & (len(dts[1]) == len(dts[2])) ): dt_2Dhists[idx].fill(ak.ravel(dts[1]),ak.ravel(dts[2]))


def draw_all(hists, runs, chan_type, savedir):
	hep.style.use("CMS")
	for idx, hist in enumerate(hists):
		fig, ax = plt.subplots()
		row = int(idx/9)
		col = int(idx/3) - 3*row
		layer = idx%3 + 1
		chan = 24*layer + 2*row + 8*(2-col) + chan_type
		l0chan = 2*row + 8*(2-col) + chan_type
		# draw the histogram
		hist.plot(yerr=False)
		ax.set_title(f"milliQan Preliminary",loc="left")
		ax.set_title(f"Slab Detector",loc="right")
		ax.grid(True)
		#plt.yscale('log')
		ax.set_ylabel('Number of Pulses')
		if len(runs) > 1: ax.text(0.69,0.93,f"Runs {runs[0]}-{runs[-1]}", transform=ax.transAxes)
		else: ax.text(0.69,0.93,f"Run {runs[0]}", transform=ax.transAxes)
		ax.text(0.69,0.88,rf'Channels {l0chan}/{chan}',transform=ax.transAxes)
		#if(chan_type == 0): ax.text(0.69,0.88,rf'Even Channels',transform=ax.transAxes)
		#if(chan_type == 1): ax.text(0.69,0.88,rf'Odd Channels',transform=ax.transAxes)
		ax.text(0.69,0.83,rf'Row {row}',transform=ax.transAxes)
		ax.text(0.69,0.78,rf'Column {col}',transform=plt.gca().transAxes)
		ax.set_xlabel(f"Time Layer 0 - Time Layer {layer} [ns]",loc='center')
		fig.tight_layout()
		figname = f"run{runs[0]}_dt_layer{layer}_row{row}_col{col}_chantype{chan_type}"
		fig.savefig(f"{savedir}/{figname}.png")
		plt.close(fig)

def write_to_root(file, runs, hists, name_type, chan_type):
	for idx, hist in enumerate(hists):
		row = int(idx/9)
		col = int(idx/3) - 3*row
		layer = idx%3 + 1
		chan = 24*layer + 2*row + 8*(2-col) + chan_type
		l0chan = 2*row + 8*(2-col) + chan_type
		if name_type == 0: name = f"run{runs[0]}_dt_layer{layer}_row{row}_col{col}_chantype{chan_type}"
		if name_type == 1: name = f"run{runs[0]}_time2D_l0l{layer}_row{row}_col{col}_chantype{chan_type}"
		if name_type == 2: name = f"run{runs[0]}_dt2D_type{layer}_row{row}_col{col}_chantype{chan_type}"
		file[name] = hist

def make_plot(runs, plot_var, plot_range, layer, row, col, chan_type, cuts, do_fit=False, save_all=False, fit_range=None, fit_func=gaussian):
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
    if(plot_var): branches_needed = [plot_var]
    else: branches_needed = []
    for (cutvar, _, _) in cut_specs:
        if cutvar not in branches_needed:
            branches_needed.append(cutvar)
    branches_needed.append("height")
    branches_needed.append("area")
    branches_needed.append("layer")
    branches_needed.append("row")
    branches_needed.append("column")
    branches_needed.append("timeFit")
    branches_needed.append("chan")
    branches_needed.append("ipulse")
    branches_needed.append("event")
    branches_needed.append("fileNumber")

    # build the histogram(s)
    dt_axis = hist.axis.Regular(50, rlow, rhigh, label=" ")
    time_axis = hist.axis.Regular(80, 1200, 1400, label=" ")
    if(save_all):
	    dt_hists = [hist.Hist(dt_axis, storage=hist.storage.Weight()) for _ in range(36)]
	    time_2Dhists = [hist.Hist(time_axis, time_axis, storage=hist.storage.Weight()) for _ in range(36)]
	    dt_2Dhists = [hist.Hist(dt_axis, dt_axis, storage=hist.storage.Weight()) for _ in range(36)]
    else: h = hist.Hist(dt_axis)

    run_list = []
    for run in runs:
	    if(run>=1105): run_list.append(f"/net/cms18/cms18r0/milliqan/Run3Offline/v36/slab/MilliQanSlab_Run{run}.*_v36.root:t")
	    else: run_list.append(f"/net/cms18/cms18r0/milliqan/Run3Offline/v35/slab/MilliQanSlab_Run{run}.*_v35.root:t")


    for branches in uproot.iterate(run_list, branches_needed, step_size=1000):

        pulses = ak.zip(
			{
				"height": branches["height"],
				"area": branches["area"],
				"row": branches["row"],
				"column": branches["column"],
				"layer": branches["layer"],
				"chan": branches["chan"],
				"timeFit": branches["timeFit"],
				"ipulse": branches["ipulse"],
			}
	)
        events = branches["event"]
        fileNumber = branches["fileNumber"]

        #Correct for mismatched columns/rows if necessary
        ch29_mask = (pulses.chan == 29)
        ch29_corrected = ak.where(ch29_mask, pulses.column + 1, pulses.column)
        pulses = ak.with_field(pulses, ch29_corrected, "column")
        ch64_mask = (pulses.chan == 64)
        ch64_corrected = ak.where(ch64_mask, pulses.row - 2, pulses.row)
        pulses = ak.with_field(pulses, ch64_corrected, "row")
        ch65_mask = (pulses.chan == 65)
        ch65_corrected = ak.where(ch65_mask, pulses.row - 2, pulses.row)
        pulses = ak.with_field(pulses, ch65_corrected, "row")

        cleaning_cuts = (pulses.height > 50) & (pulses.timeFit > 1000) & (pulses.timeFit < 1500) & (pulses.ipulse==0)#Start off with some cleaning cuts 
        straightLineMask, straightLinePulses= straightLinePath_with_dt(pulses[cleaning_cuts], chan_type) #Returns pulses which are in a straight line 

        l0mask = straightLinePulses.layer == 0
        l1mask = straightLinePulses.layer == 1
        l2mask = straightLinePulses.layer == 2
        l3mask = straightLinePulses.layer == 3
        chanmask = (straightLinePulses.chan % 2) == chan_type
        if(save_all):
                add_all(straightLinePulses, chanmask, l0mask, dt_hists, time_2Dhists, dt_2Dhists) 
                continue


        rowmask = straightLinePulses.row == row
        colmask = straightLinePulses.column == col
        p0 = straightLinePulses[l0mask & rowmask & colmask & chanmask]
        p1 = straightLinePulses[l1mask & rowmask & colmask & chanmask]
        p2 = straightLinePulses[l2mask & rowmask & colmask & chanmask]
        p3 = straightLinePulses[l3mask & rowmask & colmask & chanmask]
        if layer == 1:
                dt = p0.timeFit - p1.timeFit
                print(dt)
                arr_plot = dt
        if layer == 2:
                dt = p0.timeFit - p2.timeFit
                arr_plot = dt
        if layer == 3:
                dt = p0.timeFit - p3.timeFit
                arr_plot = dt
                #print(ak.to_list(dt))
        # Start with a mask of “select everything”
        #cut_mask = ak.ones_like(l0l1_dt, dtype=bool)
        # apply additional cuts
        #for (cutvar, cut_low, cut_high) in cut_specs:
        #    cut_mask = cut_mask & ((straightLinePulses[cutvar] >= cut_low) & (straightLinePulses[cutvar] <= cut_high))

        # fill histogram with masked values
        #print("File, Event, dt")
        #print(fileNumber[ak.any(straightLineMask[(l1mask & rowmask & colmask & chanmask)],axis=1)], events[ak.any(straightLineMask[(l1mask & rowmask & colmask & chanmask)],axis=1)], dt)
        h.fill(ak.ravel(arr_plot)) 

    savedir = "dt_plots_beamOn"
    outfile = f"dt_hists_run{runs[0]}_chantype{chan_type}.root"
    if(save_all): 
        draw_all(dt_hists, runs, chan_type, save_dir)
        #with uproot.recreate(outfile) as outfile:
        #        write_to_root(outfile, runs, dt_hists, 0, chan_type)
        #        write_to_root(outfile, runs, time_2Dhists, 1, chan_type)
        #        write_to_root(outfile, runs, dt_2Dhists, 2, chan_type)
        return
    # draw the histogram
    hep.style.use("CMS")
    h.plot(yerr=False)
    if do_fit:
        popt, pcov, xfit, yfit = fit(fit_func, h, fit_range)
        plt.plot(xfit,yfit)
    plt.title(f"milliQan Preliminary",loc="left")
    plt.title(f"Slab Detector",loc="right")
    plt.grid(True)
    #plt.yscale('log')
    plt.ylabel('Number of Pulses')
    if len(runs) > 1: plt.text(0.70,0.93,f"Runs {runs[0]}-{runs[-1]}", transform=plt.gca().transAxes)
    else: plt.text(0.70,0.93,f"Run {runs[0]}", transform=plt.gca().transAxes)
    #if 'chan' in locals(): plt.text(0.73,0.88,f"Channel {chan}",transform=plt.gca().transAxes)
    if do_fit:
        plt.text(0.70,0.83,rf'$\mu$={round(popt[1],1)}',transform=plt.gca().transAxes)
        plt.text(0.70,0.78,rf'$\sigma$={round(popt[2],1)}',transform=plt.gca().transAxes)
    chan = 24*layer + 2*row + 8*(2-col) + chan_type
    l0chan = 2*row + 8*(2-col) + chan_type
    plt.text(0.70,0.88,rf'Channels {l0chan}/{chan}',transform=plt.gca().transAxes)
    #if(chan_type == 0): plt.text(0.70,0.88,rf'Even Channels',transform=plt.gca().transAxes)
    #if(chan_type == 1): plt.text(0.70,0.88,rf'Odd Channels',transform=plt.gca().transAxes)
    plt.text(0.70,0.83,rf'Row {row}',transform=plt.gca().transAxes)
    plt.text(0.70,0.78,rf'Column {col}',transform=plt.gca().transAxes)
    if(layer): plt.xlabel(f"Time Layer 0 - Time Layer {layer} [ns]",loc='center')
    elif(plot_var == "height"): plt.xlabel("Pulse Height [mV]",loc='center')
    elif(plot_var == "area"): plt.xlabel("Pulse Area [pVs]",loc='center')
    else: plt.xlabel(plot_var,loc='center')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    make_plot(args.runs, args.plot_var, args.plot_range, args.layer, args.row, args.col, args.chan_type, args.cuts, args.do_fit, args.save_all, args.fit_range)

