import os,sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT as r
r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(1)
r.gErrorIgnoreLevel = r.kWarning

# fin = r.TFile("/nfs-7/userdata/bemarsh/milliqan/pmt_calib/processed/r878h_1250V_2p5V_20ns_300Hz_50000evnts.root")
fin = r.TFile("/home/rsantos/Data/MilliQan_Run505.3_default_peak_template.root")
tree = fin.Get("Events")

# if "878" in fin.GetName():
#     pmt = "r878"
# elif "7725" in fin.GetName():
#     pmt = "r7725"
# else:
#     raise Exception("unrecognized PMT type!")
pmt = "r878"
if pmt=="r878":
    plotdir = "/home/rsantos/Data/Plots"
if pmt=="r7725":
    plotdir = "/home/users/bemarsh/public_html/milliqan/pmt_calib/templates/r7725_1600V"
os.system("mkdir -p "+plotdir)

if pmt=="r878":
    vrange = (25,35)
    tstart = 280
    tend = 390
elif pmt=="r7725":
    vrange = (400,550)
    tstart = 280
    tend = 360

control_range=(-5,3)


# establish a common set of time points to sync waveforms
# Should be fast; points will be interpolated by the code
# Use 0.2 ns (5 GHz) here. May need to change upper time bound
t_baseline = np.arange(0, 500+1e-10, 0.2)

# frequencies at which to store final template
output_freqs = [0.7, 1, 2, 5]

# make area histogram for reference
print("Histogramming pulse areas")
c = r.TCanvas("c","c",604,528)
if pmt=="r878":
    h = r.TH1D("h",";pulse area [pVs];Events / 2 pVs", 125, -50, 200)
elif pmt=="r7725":
    h = r.TH1D("h",";pulse area [pVs];Events / 2 pVs", 125, -50, 750)
h.SetLineColor(r.kBlack)
h.SetLineWidth(2)
tree.Draw("area>>h","")
line = r.TLine()
line.SetLineStyle(2)
line.SetLineWidth(2)
line.SetLineColor(r.kRed)
m = h.GetMaximum()
line.DrawLine(control_range[0], 0, control_range[0], m)
line.DrawLine(control_range[1], 0, control_range[1], m)
line.SetLineColor(r.kBlue)
line.DrawLine(vrange[0], 0, vrange[0], m)
line.DrawLine(vrange[1], 0, vrange[1], m)
c.SaveAs(os.path.join(plotdir, "pulse_area.png"))
c.SaveAs(os.path.join(plotdir, "pulse_area.pdf"))


# first get average "control" waveform (i.e., 0PE, just background)
print("Loop over tree to get average 0-PE control waveform")
tree.Draw(">>entrylist", "area>{0} && area<{1}".format(*control_range), "goff")
el = r.gDirectory.Get("entrylist")
nent = el.GetN()
for i in range(nent):
    tree.GetEntry(el.GetEntry(i))

    ts = np.array(list(tree.times))
    vs = -np.array(list(tree.voltages))

    # interpolate into the predefined baseline time points
    vs = np.interp(t_baseline, ts, vs)

    if i==0:
        avg_control = np.array(vs)
    else:
        avg_control += vs

avg_control /= nent

# now get average SPE pulse, first subtracting off average control from each waveform
print("Loop over tree to get average control-subtracted SPE waveform")
tree.Draw(">>entrylist", "area>{0} && area<{1}".format(*vrange), "goff")
el = r.gDirectory.Get("entrylist")
nent = el.GetN()
waveforms = []
for i in range(nent):
    tree.GetEntry(el.GetEntry(i))

    ts = np.array(list(tree.times))
    vs = -np.array(list(tree.voltages))

    # interpolate into the predefined baseline time points
    vs = np.interp(t_baseline, ts, vs)
    
    waveforms.append(vs)

    if i==0:
        avg_spe = np.array(vs)
        avg_spe_fix = np.array(vs - avg_control)
    else:
        avg_spe += vs
        avg_spe_fix += vs - avg_control

avg_spe /= nent
avg_spe_fix /= nent

plt.figure(figsize=(6,5), dpi=100)
plt.plot(t_baseline, avg_control, 'r-', label="Average 0-PE control")
plt.plot(t_baseline, avg_spe, 'b-', label="Average SPE")
plt.plot([tstart]*2, [np.amin(avg_spe), np.amax(avg_spe)], 'k--')
plt.plot([tend]*2, [np.amin(avg_spe), np.amax(avg_spe)], 'k--')
plt.gca().set_ylim(ymax=1.3*np.amax(avg_spe))
plt.legend(loc='upper left')
plt.savefig(os.path.join(plotdir, "avg_spe_control.png"))
plt.savefig(os.path.join(plotdir, "avg_spe_control.pdf"))

plt.figure(figsize=(6,5), dpi=100)
plt.plot(t_baseline, avg_spe_fix, 'g-', label="Average corrected-SPE")
plt.plot([tstart]*2, [np.amin(avg_spe_fix), np.amax(avg_spe_fix)], 'k--')
plt.plot([tend]*2, [np.amin(avg_spe_fix), np.amax(avg_spe_fix)], 'k--')
dt = 100 if pmt=="r878" else 50
plt.gca().set_xlim(tstart-dt, tend+dt)
plt.gca().set_ylim(ymax=1.3*np.amax(avg_spe_fix))
plt.legend(loc='upper left')
plt.savefig(os.path.join(plotdir, "avg_spe_fixed.png"))
plt.savefig(os.path.join(plotdir, "avg_spe_fixed.pdf"))

def get_template(avg):
    imin = np.argmax(t_baseline > tstart)
    imax = np.argmax(t_baseline > tend)
    template = np.array(avg[imin:imax])
    template /= np.sum(template)
    return template

def make_event_display(avg, vs, draw_smoothed=False, saveAs=None):
    template = get_template(avg)
    smoothed = np.convolve(template[::-1], vs, mode='valid')
    offset = np.argmax(template)
    i_peak = offset + np.argmax(smoothed)

    fitted = template * np.trapz(vs[i_peak-offset:i_peak-offset+template.size], t_baseline[:template.size]) / np.trapz(template, t_baseline[:template.size])
    mults = np.linspace(0.5,2.0,101)
    sses = []
    for a in mults:
        test = fitted * a
        sses.append(np.sum((vs[i_peak-offset:i_peak-offset+template.size] - test)**2))    
    mult = mults[np.argmax(-np.array(sses))]
    fitted *= mult

    plt.plot(t_baseline, vs, '-', color="0.6")
    if draw_smoothed:
        plt.plot(t_baseline[offset:offset+smoothed.size], smoothed, 'r-', lw=2)
    plt.plot(t_baseline[i_peak-offset:i_peak-offset+fitted.size], fitted, 'g-', lw=2)
    plt.gca().set_xlim(tstart-70, tend+70)
    plt.gca().set_ylim(ymin = -2 if pmt=="r878" else -10, ymax=np.amax(vs)*1.2)
    plt.savefig(os.path.join(plotdir, saveAs+".png"))
    plt.savefig(os.path.join(plotdir, saveAs+".pdf"))
    
vs = waveforms[0] - avg_control
plt.figure(figsize=(6,5), dpi=100)
make_event_display(avg_spe_fix, vs, draw_smoothed=True, saveAs="samp_smoothed")

os.system("mkdir -p "+plotdir+"/events")
for i in range(10):
    plt.clf()
    vs = waveforms[i] - avg_control
    make_event_display(avg_spe_fix, vs, draw_smoothed=False, saveAs="events/event{0:03d}".format(i))


print("Loop over tree to get average time-corrected SPE waveform")
template = get_template(avg_spe_fix)
offset = np.argmax(template)
t_peaks = []
for i in range(nent):
    vs = waveforms[i] - avg_control
    
    # smooth to find the peak
    smoothed = np.convolve(template[::-1], vs, mode='valid')
    t_peak = t_baseline[offset+np.argmax(smoothed)]
    t_peaks.append(t_peak)

    # shift in time and re-interpolate into the baseline time points
    vs = np.interp(t_baseline, t_baseline-(t_peak-t_peaks[0]), vs)
    if i==0:
        avg_spe_tcorr = np.array(vs)
    else:
        avg_spe_tcorr += vs

avg_spe_tcorr /= nent

if pmt=="r878":
    # weird noise issue in 878 after averaging... smooth out here
    ones = np.ones(21)
    avg_spe_tcorr = np.convolve(ones, avg_spe_tcorr, mode='same')
    counts = np.convolve(ones, np.ones(avg_spe_tcorr.size), mode='same')
    avg_spe_tcorr /= counts

# plot histogram of peak time
plt.figure(figsize=(6,5), dpi=100)
plt.hist(t_peaks, range=[np.mean(t_peaks)-20, np.mean(t_peaks)+20], bins=40, histtype='step')
plt.savefig(os.path.join(plotdir, "time_spread.png"))
plt.savefig(os.path.join(plotdir, "time_spread.pdf"))

# plot time-corrected template
plt.figure(figsize=(6,5), dpi=100)
plt.plot(t_baseline, avg_spe_fix, 'r-', label="Non-time-shifted SPE average")
plt.plot(t_baseline, avg_spe_tcorr, 'g-', lw=2, label="Time-shifted SPE average")
plt.plot([tstart]*2, [np.amin(avg_spe_tcorr), np.amax(avg_spe_tcorr)], 'k--')
plt.plot([tend]*2, [np.amin(avg_spe_tcorr), np.amax(avg_spe_tcorr)], 'k--')
dt = 100 if pmt=="r878" else 50
plt.gca().set_xlim(tstart-dt, tend+dt)
plt.gca().set_ylim(ymax=1.35*np.amax(avg_spe_tcorr))
plt.legend(loc='upper left')
plt.savefig(os.path.join(plotdir, "avg_spe_fixed_tcorr.png"))
plt.savefig(os.path.join(plotdir, "avg_spe_fixed_tcorr.pdf"))

vs = waveforms[0] - avg_control
plt.figure(figsize=(6,5), dpi=100)
make_event_display(avg_spe_tcorr, vs, draw_smoothed=True, saveAs="samp_smoothed_tcorr")

os.system("mkdir -p "+plotdir+"/events_tcorr")
for i in range(10):
    plt.clf()
    vs = waveforms[i] - avg_control
    make_event_display(avg_spe_tcorr, vs, draw_smoothed=False, saveAs="events_tcorr/event{0:03d}".format(i))

template = get_template(avg_spe_tcorr)
offset = np.argmax(template)
# redefine t=0 to be the maximum of the template
t_template = t_baseline[:template.size] - t_baseline[offset]

print("Generating templates at frequencies", output_freqs, "GHz")
hs = {}
for freq in output_freqs:
    dt = 1.0 / freq
    np_left = int(abs(t_template[0]) / dt)+1
    np_right = int(abs(t_template[-1]) / dt)+1
    
    ts = np.append(
        np.linspace(-np_left*dt, 0, np_left+1),
        np.linspace(dt, np_right*dt, np_right),
        )

    out = np.interp(ts, t_template, template, left=0.0, right=0.0)
    out /= np.trapz(out, ts)
    hs[freq] = r.TH1D("h"+str(freq), ";time[ns]", out.size, 0, out.size*dt)
    for i in range(out.size):
        hs[freq].SetBinContent(i+1, out[i])

c = r.TCanvas("c2","c2",604, 528)
for freq in output_freqs:
    c.Clear()
    hs[freq].SetLineColor(r.kBlack)
    hs[freq].Draw("HIST")
    c.SaveAs(os.path.join(plotdir, "template_{0}_GHz.png".format(str(freq).replace(".","p"))))
    c.SaveAs(os.path.join(plotdir, "template_{0}_GHz.pdf".format(str(freq).replace(".","p"))))
    
