import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import hist
import argparse
import sys
import re
import glob
import ROOT as r
import pathlib
import re

#============================Four layers offline========================================
def offlineTrig1Check(pulse1,pulse2):    
    #Filter out the paired pulses which are not in the same layer and aren't panels. We don't care about those combinations anyway and removing them speeds things up.
    pulse1_filtered = pulse1[timemask & layermask & not_panels]
    pulse2_filtered = pulse2[timemask & layermask & not_panels]

    #Pairs of filtered pulses
    filtered_pairs = ak.zip([pulse1_filtered,pulse2_filtered])

    #Now we take the pair of pulses and form pairs of those pairs
    pairpairs = ak.combinations(filtered_pairs,2)

    pair1, pair2 = ak.unzip(pairpairs)

    pair1_pulse1,pair1_pulse2 = ak.unzip(pair1)
    pair2_pulse1,pair2_pulse2 = ak.unzip(pair2)

    #Now create masks that check whether all four pulses from each pair are within 150ns of each other. (Can't seem to take min or max between awkward arrays elementwise.)
    aaba = np.abs(pair1_pulse1.time - pair2_pulse1.time) <= window
    aabb = np.abs(pair1_pulse1.time - pair2_pulse2.time) <= window
    abba = np.abs(pair1_pulse2.time - pair2_pulse1.time) <= window    #see that girl...
    abbb = np.abs(pair1_pulse2.time - pair2_pulse2.time) <= window

    window_mask = aaba & aabb & abba & abbb     #Require that they all be true in order for all four pulses to be in a common window.

    #These are the combinations of four pulses which are our candidates for trigger 1
    trig1cand1 = pair1_pulse1[window_mask]
    trig1cand2 = pair1_pulse2[window_mask]
    trig1cand3 = pair2_pulse1[window_mask]
    trig1cand4 = pair2_pulse2[window_mask]

    #For our final offline four layers hit mask we require that none of our candidate pulses be in the same layer. Note that we already know cand1/3 and cand2/4 are not in the same layer.
    different_layers = ((trig1cand1.layer != trig1cand3.layer) & 
                        (trig1cand1.layer != trig1cand4.layer) & 
                        (trig1cand2.layer != trig1cand3.layer) & 
                        (trig1cand2.layer != trig1cand4.layer))

    fourLayersHitmask = different_layers
    offline_trig1_events = event[ak.any(fourLayersHitmask,axis=1)]

    #The following lines find the combination of four pulses which are earliest and chooses them to be put in the histogram
    min11 = ak.min(trig1cand1.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    min12 = ak.min(trig1cand2.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    min13 = ak.min(trig1cand3.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    min14 = ak.min(trig1cand4.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)

    argmin11 = ak.argmin(trig1cand1.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    argmin12 = ak.argmin(trig1cand2.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    argmin13 = ak.argmin(trig1cand3.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)
    argmin14 = ak.argmin(trig1cand4.time[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)],axis=1)

    arr1 = trig1cand1.chan[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)]
    arr2 = trig1cand2.chan[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)]
    arr3 = trig1cand3.chan[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)]
    arr4 = trig1cand4.chan[fourLayersHitmask][ak.any(fourLayersHitmask,axis=1)]

    minarray = np.minimum(min11,min12)
    minarray = np.minimum(minarray,min13)
    minarray = np.minimum(minarray,min14)

    which_min = ak.argmax(
        [min11 == minarray, min12 == minarray, min13 == minarray, min14 == minarray], 
        axis=0
    )
    
    indices = ak.zeros_like(minarray)
    indices = ak.where(which_min == 0, argmin11, indices)
    indices = ak.where(which_min == 1, argmin12, indices)
    indices = ak.where(which_min == 2, argmin13, indices)
    indices = ak.where(which_min == 3, argmin14, indices)
    indices_np = np.array(ak.to_numpy(indices), dtype=int)

    sliced_arr1 = ak.Array([row[index] for row, index in zip(arr1,indices_np)])
    sliced_arr2 = ak.Array([row[index] for row, index in zip(arr2,indices_np)])
    sliced_arr3 = ak.Array([row[index] for row, index in zip(arr3,indices_np)])
    sliced_arr4 = ak.Array([row[index] for row, index in zip(arr4,indices_np)])

    cand1_chans = sliced_arr1
    cand2_chans = sliced_arr2
    cand3_chans = sliced_arr3
    cand4_chans = sliced_arr4
    
    channels1.fill(ak.ravel(cand1_chans))
    channels1.fill(ak.ravel(cand2_chans))
    channels1.fill(ak.ravel(cand3_chans))
    channels1.fill(ak.ravel(cand4_chans))
    
    
    return offline_trig1_events, cand1_chans, cand2_chans, cand3_chans, cand4_chans

#============================Three in row offline========================================
def offlineTrig2Check(pulse1,pulse2):
    #Form a mask to require that the paired pulses are in the same column and trigger group row. We don't care about those combinations anyway and removing them speeds things up.
    rowmask = pulse1.row == pulse2.row
    columnmask = (pulse1.column == pulse2.column) | (pulse1.column == pulse2.column + 1) | (pulse1.column + 1 == pulse2.column) & ((pulse1.column != 1) | (pulse2.column != 2)) & ((pulse1.column != 2 ) | (pulse2.column != 1))

    pulse1_filtered = pulse1[timemask & layermask & not_panels & rowmask & columnmask]
    pulse2_filtered = pulse2[timemask & layermask & not_panels & rowmask & columnmask]

    #Pairs of filtered pulses
    filtered_pairs = ak.zip([pulse1_filtered,pulse2_filtered])

    #Now we take the pair of pulses and form pairs of those pairs
    pairpairs = ak.combinations(filtered_pairs,2)

    pair1, pair2 = ak.unzip(pairpairs)

    pair1_pulse1,pair1_pulse2 = ak.unzip(pair1)
    pair2_pulse1,pair2_pulse2 = ak.unzip(pair2)

    #Now create masks that check whether all four pulses from each pair are within 150ns of each other. (Can't seem to take min or max between awkward arrays elementwise.)
    #Note: we must also require the pulse time difference is not zero to remove duplicate pulses, provided that they in fact are duplicates. Check the channel to make sure they aren't just pulses from 
    #different channels arriving at the same time.
    aaba = (np.abs(pair1_pulse1.time - pair2_pulse1.time) <= window) & ( (pair1_pulse1.chan != pair2_pulse1.chan) | (np.abs(pair1_pulse1.time - pair2_pulse1.time) != 0) )
    aabb = (np.abs(pair1_pulse1.time - pair2_pulse2.time) <= window) & ( (pair1_pulse1.chan != pair2_pulse2.chan) | (np.abs(pair1_pulse1.time - pair2_pulse2.time) != 0) )
    abba = (np.abs(pair1_pulse2.time - pair2_pulse1.time) <= window) & ( (pair1_pulse2.chan != pair2_pulse1.chan) | (np.abs(pair1_pulse2.time - pair2_pulse1.time) != 0) )   #see that girl...
    abbb = (np.abs(pair1_pulse2.time - pair2_pulse2.time) <= window) & ( (pair1_pulse2.chan != pair2_pulse2.chan) | (np.abs(pair1_pulse2.time - pair2_pulse2.time) != 0) )

    #Require that at least three of the pulses be within the time window. Note that we already know pair1pulse1 and pair2pulse2 satisfy the window, ditto for pair2pulse1 and pair2pulse2
    time_123 = abba & aaba
    time_124 = aabb & abbb
    time_234 = abba & abbb
    time_134 = aaba & aabb

    #There are four different ways to have three different layers. Note that we have already required that pair1_pulse1 and pair1_pulse2 be in different layers from the first pairing, same for cand3 and cand4.
    diff_layer123 =  ((pair1_pulse1.layer != pair2_pulse1.layer) &
                    (pair1_pulse2.layer != pair2_pulse1.layer))
    diff_layer124 =  ((pair1_pulse1.layer != pair2_pulse2.layer) &
                    (pair1_pulse2.layer != pair2_pulse2.layer))
    diff_layer234 =  ((pair1_pulse2.layer != pair2_pulse1.layer) & 
                    (pair1_pulse2.layer != pair2_pulse2.layer))
    diff_layer134 =  ((pair1_pulse1.layer != pair2_pulse1.layer) & 
                    (pair1_pulse1.layer != pair2_pulse2.layer))

    #There are four ways to be in the same row. Again, cand1 and cand2 are already in the same row, as are cand3 and cand4.
    same_row123 =  ((pair1_pulse1.row == pair2_pulse1.row) &
            (pair1_pulse2.row == pair2_pulse1.row))
    same_row124 =  ((pair1_pulse1.row == pair2_pulse2.row) &
            (pair1_pulse2.row == pair2_pulse2.row))
    same_row234 =  ((pair1_pulse2.row == pair2_pulse1.row) &
            (pair1_pulse2.row == pair2_pulse2.row))
    same_row134 =  ((pair1_pulse1.row == pair2_pulse1.row) &
            (pair1_pulse1.row == pair2_pulse2.row))

    #There are four ways to be in the same trigger group (column pairs). Once again cand1 and cand2 already are as are cand3 and cand4.
    same_cols123 = (((pair1_pulse1.column == pair2_pulse1.column) | (pair1_pulse1.column == pair2_pulse1.column + 1) | (pair1_pulse1.column + 1 == pair2_pulse1.column)) & 
                    (((pair1_pulse1.column != 1) | (pair2_pulse1.column != 2)) & ((pair1_pulse1.column != 2 ) | (pair2_pulse1.column != 1))) &
                    ((pair1_pulse2.column == pair2_pulse1.column) | (pair1_pulse2.column == pair2_pulse1.column + 1) | (pair1_pulse2.column + 1 == pair2_pulse1.column)) & 
                    (((pair1_pulse2.column != 1) | (pair2_pulse1.column != 2)) & ((pair1_pulse2.column != 2 ) | (pair2_pulse1.column != 1))))
    same_cols124 = (((pair1_pulse1.column == pair2_pulse2.column) | (pair1_pulse1.column == pair2_pulse2.column + 1) | (pair1_pulse1.column + 1 == pair2_pulse2.column)) & 
                    (((pair1_pulse1.column != 1) | (pair2_pulse2.column != 2)) & ((pair1_pulse1.column != 2 ) | (pair2_pulse2.column != 1))) &
                    ((pair1_pulse2.column == pair2_pulse2.column) | (pair1_pulse2.column == pair2_pulse2.column + 1) | (pair1_pulse2.column + 1 == pair2_pulse2.column)) & 
                    (((pair1_pulse2.column != 1) | (pair2_pulse2.column != 2)) & ((pair1_pulse2.column != 2 ) | (pair2_pulse2.column != 1))))
    same_cols234 = (((pair1_pulse2.column == pair2_pulse1.column) | (pair1_pulse2.column == pair2_pulse1.column + 1) | (pair1_pulse2.column + 1 == pair2_pulse1.column)) & 
                    (((pair1_pulse2.column != 1) | (pair2_pulse1.column != 2)) & ((pair1_pulse2.column != 2 ) | (pair2_pulse1.column != 1))) &
                    ((pair1_pulse2.column == pair2_pulse2.column) | (pair1_pulse2.column == pair2_pulse2.column + 1) | (pair1_pulse2.column + 1 == pair2_pulse2.column)) & 
                    (((pair1_pulse2.column != 1) | (pair2_pulse2.column != 2)) & ((pair1_pulse2.column != 2 ) | (pair2_pulse2.column != 1))))
    same_cols134 = (((pair1_pulse1.column == pair2_pulse1.column) | (pair1_pulse1.column == pair2_pulse1.column + 1) | (pair1_pulse1.column + 1 == pair2_pulse1.column)) & 
                    (((pair1_pulse1.column != 1) | (pair2_pulse1.column != 2)) & ((pair1_pulse1.column != 2 ) | (pair2_pulse1.column != 1))) &
                    ((pair1_pulse1.column == pair2_pulse2.column) | (pair1_pulse1.column == pair2_pulse2.column + 1) | (pair1_pulse1.column + 1 == pair2_pulse2.column)) & 
                    (((pair1_pulse1.column != 1) | (pair2_pulse2.column != 2)) & ((pair1_pulse1.column != 2 ) | (pair2_pulse2.column != 1))))

    #We now have four possible candidate three in row combinations, we take the OR of them to allow any of them to be true
    trig2cand1 = time_123 & diff_layer123 & same_row123 & same_cols123
    trig2cand2 = time_124 & diff_layer124 & same_row124 & same_cols124
    trig2cand3 = time_234 & diff_layer234 & same_row234 & same_cols234
    trig2cand4 = time_134 & diff_layer134 & same_row134 & same_cols134
    threeInRow_mask = trig2cand1 | trig2cand2 | trig2cand3 | trig2cand4

    offline_trig2_events = event[ak.any(threeInRow_mask,axis=1)]
    return offline_trig2_events

#============================Two separated layers========================================
def offlineTrig3Check(pulse1,pulse2):
    nonadjacent = layerdiff > 1
    twoSeparatedLayers_mask = nonadjacent & not_panels & timemask
    offline_trig3_events = event[ak.any(twoSeparatedLayers_mask,axis=1)]
    return offline_trig3_events

#============================Two adjacent layers=========================================
def offlineTrig4Check(pulse1,pulse2):
    adjacent = layerdiff == 1
    twoAdjacentLayers_mask = adjacent & not_panels & timemask
    offline_trig4_events = event[ak.any(twoAdjacentLayers_mask,axis=1)]
    return offline_trig4_events


#============================N Layers Hit================================================
def offlineTrig5Check(pulse1,pulse2):
    pulse1_filtered = pulse1[timemask & layermask]
    pulse2_filtered = pulse2[timemask & layermask]

    #Pairs of filtered pulses
    filtered_pairs = ak.zip([pulse1_filtered,pulse2_filtered])

    #Now we take the pair of pulses and form pairs of those pairs
    pairpairs = ak.combinations(filtered_pairs,2)

    pair1, pair2 = ak.unzip(pairpairs)

    pair1_pulse1,pair1_pulse2 = ak.unzip(pair1)
    pair2_pulse1,pair2_pulse2 = ak.unzip(pair2)

    #Now create masks that check whether all four pulses from each pair are within 150ns of each other. (Can't seem to take min or max between awkward arrays elementwise.)
    #Note: we must also require the pulse time difference is not zero to remove duplicate pulses, provided that they in fact are duplicates. Check the channel to make sure they aren't just pulses from 
    #different channels arriving at the same time.
    aaba = (np.abs(pair1_pulse1.time - pair2_pulse1.time) <= window) & ( (pair1_pulse1.chan != pair2_pulse1.chan) | (np.abs(pair1_pulse1.time - pair2_pulse1.time) != 0) )
    aabb = (np.abs(pair1_pulse1.time - pair2_pulse2.time) <= window) & ( (pair1_pulse1.chan != pair2_pulse2.chan) | (np.abs(pair1_pulse1.time - pair2_pulse2.time) != 0) )
    abba = (np.abs(pair1_pulse2.time - pair2_pulse1.time) <= window) & ( (pair1_pulse2.chan != pair2_pulse1.chan) | (np.abs(pair1_pulse2.time - pair2_pulse1.time) != 0) )   #see that girl...
    abbb = (np.abs(pair1_pulse2.time - pair2_pulse2.time) <= window) & ( (pair1_pulse2.chan != pair2_pulse2.chan) | (np.abs(pair1_pulse2.time - pair2_pulse2.time) != 0) )

    #Require that at least three of the pulses be within the time window. Note that we already know pair1pulse1 and pair2pulse2 satisfy the window, ditto for pair2pulse1 and pair2pulse2
    time_123 = abba & aaba
    time_124 = aabb & abbb
    time_234 = abba & abbb
    time_134 = aaba & aabb

    #There are four different ways to have three different layers. Note that we have already required that pair1_pulse1 and pair1_pulse2 be in different layers from the first pairing, same for cand3 and cand4.
    diff_layer123 =  ((pair1_pulse1.layer != pair2_pulse1.layer) &
                    (pair1_pulse2.layer != pair2_pulse1.layer))
    diff_layer124 =  ((pair1_pulse1.layer != pair2_pulse2.layer) &
                    (pair1_pulse2.layer != pair2_pulse2.layer))
    diff_layer234 =  ((pair1_pulse2.layer != pair2_pulse1.layer) & 
                    (pair1_pulse2.layer != pair2_pulse2.layer))
    diff_layer134 =  ((pair1_pulse1.layer != pair2_pulse1.layer) & 
                    (pair1_pulse1.layer != pair2_pulse2.layer))

    #We now have four possible candidate three layers combinations, we take the OR of them to allow any of them to be true
    trig5cand1 = time_123 & diff_layer123
    trig5cand2 = time_124 & diff_layer124
    trig5cand3 = time_234 & diff_layer234
    trig5cand4 = time_134 & diff_layer134
    nLayers_mask = trig5cand1 | trig5cand2 | trig5cand3 | trig5cand4

    offline_trig5_events = event[ak.any(nLayers_mask,axis=1)]
    return offline_trig5_events

#============================Greater than N Hits=========================================
#We assume nHit = 3 so that we are always checking for at least 4 in time pulses. If nHit changes, will need to edit this code.
def offlineTrig7Check(pulse1,pulse2):
    #Filter out the paired pulses which are not in time.
    pulse1_filtered = pulse1[timemask]
    pulse2_filtered = pulse2[timemask]

    #Pairs of filtered pulses
    filtered_pairs = ak.zip([pulse1_filtered,pulse2_filtered])

    #Now we take the pair of pulses and form pairs of those pairs
    pairpairs = ak.combinations(filtered_pairs,2)

    pair1, pair2 = ak.unzip(pairpairs)

    pair1_pulse1,pair1_pulse2 = ak.unzip(pair1)
    pair2_pulse1,pair2_pulse2 = ak.unzip(pair2)

    #Now create masks that check whether all four pulses from each pair are within 150ns of each other. (Can't seem to take min or max between awkward arrays elementwise.)
    aaba = np.abs(pair1_pulse1.time - pair2_pulse1.time) <= window
    aabb = np.abs(pair1_pulse1.time - pair2_pulse2.time) <= window
    abba = np.abs(pair1_pulse2.time - pair2_pulse1.time) <= window    #see that girl...
    abbb = np.abs(pair1_pulse2.time - pair2_pulse2.time) <= window

    window_mask = aaba & aabb & abba & abbb     #Require that they all be true in order for all four pulses to be in a common window.

    gtNHits_mask = ak.any(window_mask,axis=1)
    offline_trig7_events = event[gtNHits_mask]
    return offline_trig7_events

#============================Top Panels==================================================
def offlineTrig9Check(pulses):
    topPanels_mask = pulses.row == 4
    offline_trig9_events = event[ak.any(topPanels_mask,axis=1)]
    return offline_trig9_events

#============================Top Panels + Bottom Bars====================================
def offlineTrig10Check(pulse1,pulse2):
    panel_bar_mask = ((pulse1.row == 4) & ((pulse2.row == 0) & (pulse2["type"] == 0))) | ( (pulse2.row == 4) & ((pulse1.row == 0) & (pulse1["type"] == 0)) )
    topPanelsBotBars_mask = panel_bar_mask & timemask
    #ak.any(pulses.row == 4,axis=1) & ak.any( (pulses.row == 0) & (pulses["type"] == 0),axis=1)
    offline_trig10_events = event[ak.any(topPanelsBotBars_mask,axis=1)]
    return offline_trig10_events

#============================Front/Back Panels==================================================
def offlineTrig11Check(pulse1,pulse2):
    frontBack_mask = ((pulse1.layer == -1) & (pulse2.layer == 4)) | ((pulse2.layer == -1) & (pulse1.layer == 4))
    #ak.any(pulses.layer == -1, axis=1) & ak.any(pulses.layer == 4, axis=1)
    offline_trig11_events = event[ak.any(frontBack_mask & timemask,axis=1)]
    return offline_trig11_events


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


path = sys.argv[1]  #Path to data file to run on 
quality_string = "loose"  #goodRunsList quality string
files = sorted(glob.glob("{0}.root".format(path)), key=natural_sort_key)
path_to_json = "goodRunsList.json"
goodJson_array = ak.from_json(pathlib.Path(path_to_json))
data = ak.Array(goodJson_array['data'])
goodJson = ak.zip({
            'run': data[:, 0],
            'file': data[:, 1],
            'loose': data[:, 2],
            'medium': data[:, 3],
            'tight': data[:, 4],
            'single_trigger': data[:, 5]
        }, depth_limit=1)

goodJson = ak.enforce_type(goodJson, "{run: int64, file: int64, loose: bool, medium: bool, tight: bool, single_trigger: bool}")
pattern = r"Run(\d+)\.(\d+)"

#Remove files which fail quality requirements
for file in files[:]:
    match = re.search(pattern, file)
    if match:
        run_num = int(match.group(1))
        file_num = int(match.group(2))
        matching_goodJson = goodJson[(goodJson['run'] == int(run_num)) & (goodJson['file'] == int(file_num))]
        if len(matching_goodJson) == 0 or not matching_goodJson[quality_string][0]:
            files.remove(file)
    else:
        raise TypeError("Filename does not match pattern.")
    
    


#print(files)
fig = plt.figure()
channels1 = hist.Hist(hist.axis.Regular(64,0,64,label="Channel"))
channels1_and_2 = hist.Hist(hist.axis.Regular(64,0,64,label="Channel"))
channels1_not_2 = hist.Hist(hist.axis.Regular(64,0,64,label="Channel"))
ZBpulses = hist.Hist(hist.axis.Regular(20,0,20,label="N Pulses"))
fourLayersPulses = hist.Hist(hist.axis.Regular(20,0,20,label="N Pulses"))
threeInRowPulses = hist.Hist(hist.axis.Regular(20,0,20,label="N Pulses"))

offline_trig1_events = ak.Array([])
offline_trig2_events = ak.Array([])
offline_trig3_events = ak.Array([])
offline_trig4_events = ak.Array([])
offline_trig5_events = ak.Array([])
offline_trig7_events = ak.Array([])
offline_trig9_events = ak.Array([])
offline_trig10_events = ak.Array([])
offline_trig11_events = ak.Array([])

online_trig1_events = ak.Array([])
online_trig2_events = ak.Array([])
online_trig3_events = ak.Array([])
online_trig4_events = ak.Array([])
online_trig5_events = ak.Array([])
online_trig7_events = ak.Array([])
online_trig9_events = ak.Array([])
online_trig10_events = ak.Array([])
online_trig11_events = ak.Array([])
online_trig13_events = ak.Array([])

empty_trig1_events = ak.Array([])
empty_trig2_events = ak.Array([])
empty_trig3_events = ak.Array([])
empty_trig4_events = ak.Array([])
empty_trig5_events = ak.Array([])
empty_trig7_events = ak.Array([])
empty_trig9_events = ak.Array([])
empty_trig10_events = ak.Array([])
empty_trig11_events = ak.Array([])
empty_trig13_events = ak.Array([])

n_zerobias_pulses = ak.Array([])
n_trig1_pulses = ak.Array([])
n_trig2_pulses = ak.Array([])


files_with_trees = {file_name: "t;1" for file_name in files}
#print(files_with_trees)

for branches in uproot.iterate(files_with_trees,["time","height","area","row","column","layer","chan","type","event","tTrigger","dynamicPedestal","fileNumber","runNumber"],step_size=1000):


    #Open root file and ttree
    #file = uproot.open("MilliQan_Run1415.2_v34.root")
    #tree = file["t;1"]
    #stop = 1000   #Set the number of events to run on
    #branches = tree.arrays(["time","height","area","row","column","layer","duration","chan","type"],entry_stop=stop)

    
    #Fix mislabeled channels (don't use chan!!!!)
    chan_mask78 = (branches["chan"] == 78)
    chan_mask79 = (branches["chan"] == 79)
    new_chan = ak.where(chan_mask78, 24, branches["chan"])
    new_chan = ak.where(chan_mask79, 25, new_chan)
    branches = ak.with_field(branches, new_chan, "chan")

    #tTrigger = tree["tTrigger"].array(entry_stop=stop)
    #event = tree["event"].array(entry_stop=stop)
    runNumber = branches["runNumber"]
    fileNumber = branches["fileNumber"]
    tTrigger = branches["tTrigger"]
    event = branches["event"]
    dynamicPedestal = branches["dynamicPedestal"]

    matched_mask = tTrigger != -1  #Require matched triggers
    non_empty_mask = ak.num(branches["chan"]) > 0
    is_empty_mask = ak.num(branches["chan"]) == 0
    if(ak.any(is_empty_mask)):
        empty_branches = branches[is_empty_mask]
        empty_event = event[is_empty_mask]
        empty_tTrigger = tTrigger[is_empty_mask]
    branches = branches[matched_mask & non_empty_mask]
    dynamicPedestal = dynamicPedestal[matched_mask & non_empty_mask]
    tTrigger = tTrigger[matched_mask & non_empty_mask]
    event = event[matched_mask & non_empty_mask]

    #Turn the decimal tTrigger branch into an array of bitstrings
    bin_rep_vec = np.vectorize(np.binary_repr)
    trig_np = ak.to_numpy(tTrigger).astype(int)
    triggerbits = bin_rep_vec(trig_np,width=13)
    if(ak.any(is_empty_mask)):
        empty_trig_np = ak.to_numpy(empty_tTrigger).astype(int)
        empty_triggerbits = bin_rep_vec(empty_trig_np,width=13)

    #Now zip all these pulse shaped branches together into a record called pulses
    pulses = ak.zip(
        {
            "time": branches["time"],
            "height": branches["height"],
            "area": branches["area"],
            "chan": branches["chan"],
            "row": branches["row"],
            "column": branches["column"],
            "layer": branches["layer"],
            "type": branches["type"],
        }
    )

    

    online_height = pulses.height + dynamicPedestal[pulses.chan]
    heightmask = online_height > 15
    pulses = pulses[heightmask]  #Apply dynamic pedestal correction to pulses



    #Now for timing combinations. We use ak.combinations to do N choose 2 on the existing pulses to form pairs of pulses.
    pairs = ak.combinations(pulses,2)

    #Unzip the pairs to get each half of the combinations.
    pulse1, pulse2 = ak.unzip(pairs)

    #Form a mask for pulses within 150ns of each other. This mask will cut down on the combinatorics.
    window = 150
    timemask = np.abs(pulse1.time - pulse2.time) <= window

    #Useful for later
    layerdiff = np.abs(pulse1.layer - pulse2.layer)
    layermask = layerdiff != 0   #pulse1.layer != pulse2.layer
    not_panels = (pulse1["type"] == 0) & (pulse2["type"] == 0)

    file_prefix = 1000*fileNumber[0] - 1000

    offline_trig1_chunk_events, cand1_chans, cand2_chans, cand3_chans, cand4_chans = offlineTrig1Check(pulse1,pulse2)
    offline_trig1_chunk_events = offline_trig1_chunk_events + file_prefix
    offline_trig2_chunk_events = offlineTrig2Check(pulse1,pulse2) + file_prefix
    offline_trig3_chunk_events = offlineTrig3Check(pulse1,pulse2) + file_prefix
    offline_trig4_chunk_events = offlineTrig4Check(pulse1,pulse2) + file_prefix
    offline_trig5_chunk_events = offlineTrig5Check(pulse1,pulse2) + file_prefix
    offline_trig7_chunk_events = offlineTrig7Check(pulse1,pulse2) + file_prefix
    offline_trig9_chunk_events = offlineTrig9Check(pulses) + file_prefix
    offline_trig10_chunk_events = offlineTrig10Check(pulse1,pulse2) + file_prefix
    offline_trig11_chunk_events = offlineTrig11Check(pulse1,pulse2) + file_prefix

    trig1_and_trig2 = np.isin(offline_trig1_chunk_events,offline_trig2_chunk_events)
    channels1_and_2.fill(cand1_chans[trig1_and_trig2])
    channels1_and_2.fill(cand2_chans[trig1_and_trig2])
    channels1_and_2.fill(cand3_chans[trig1_and_trig2])
    channels1_and_2.fill(cand4_chans[trig1_and_trig2])
    channels1_not_2.fill(cand1_chans[~trig1_and_trig2])
    channels1_not_2.fill(cand2_chans[~trig1_and_trig2])
    channels1_not_2.fill(cand3_chans[~trig1_and_trig2])
    channels1_not_2.fill(cand4_chans[~trig1_and_trig2])
    

    offline_trig1_events = ak.concatenate([offline_trig1_events,offline_trig1_chunk_events]) 
    offline_trig2_events = ak.concatenate([offline_trig2_events,offline_trig2_chunk_events])
    offline_trig3_events = ak.concatenate([offline_trig3_events,offline_trig3_chunk_events])
    offline_trig4_events = ak.concatenate([offline_trig4_events,offline_trig4_chunk_events])
    offline_trig5_events = ak.concatenate([offline_trig5_events,offline_trig5_chunk_events])
    offline_trig7_events = ak.concatenate([offline_trig7_events,offline_trig7_chunk_events])
    offline_trig9_events = ak.concatenate([offline_trig9_events,offline_trig9_chunk_events])
    offline_trig10_events = ak.concatenate([offline_trig10_events,offline_trig10_chunk_events])
    offline_trig11_events = ak.concatenate([offline_trig11_events,offline_trig11_chunk_events])


    #Online trigger events
    online_trig1_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-1] == '1'],dtype=int)] + file_prefix
    online_trig2_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-2] == '1'],dtype=int)] + file_prefix
    online_trig3_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-3] == '1'],dtype=int)] + file_prefix
    online_trig4_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-4] == '1'],dtype=int)] + file_prefix
    online_trig5_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-5] == '1'],dtype=int)] + file_prefix
    #online_trig6_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-6] == '1'],dtype=int)] + file_prefix
    online_trig7_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-7] == '1'],dtype=int)] + file_prefix
    #online_trig8_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-8] == '1'],dtype=int)] + file_prefix
    online_trig9_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-9] == '1'],dtype=int)] + file_prefix
    online_trig10_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-10] == '1'],dtype=int)] + file_prefix
    online_trig11_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-11] == '1'],dtype=int)] + file_prefix
    #online_trig12_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-12] == '1'],dtype=int)] + file_prefix
    online_trig13_chunk_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-13] == '1'],dtype=int)] + file_prefix

    online_trig1_events = ak.concatenate([online_trig1_events,online_trig1_chunk_events])
    online_trig2_events = ak.concatenate([online_trig2_events,online_trig2_chunk_events])
    online_trig3_events = ak.concatenate([online_trig3_events,online_trig3_chunk_events])
    online_trig4_events = ak.concatenate([online_trig4_events,online_trig4_chunk_events])
    online_trig5_events = ak.concatenate([online_trig5_events,online_trig5_chunk_events])
    online_trig7_events = ak.concatenate([online_trig7_events,online_trig7_chunk_events])
    online_trig9_events = ak.concatenate([online_trig9_events,online_trig9_chunk_events])
    online_trig10_events = ak.concatenate([online_trig10_events,online_trig10_chunk_events])
    online_trig11_events = ak.concatenate([online_trig11_events,online_trig11_chunk_events])
    online_trig13_events = ak.concatenate([online_trig13_events,online_trig13_chunk_events])


    #Empty events
    if(ak.any(is_empty_mask)):
        empty_trig1_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-1] == '1'],dtype=int)] + file_prefix
        empty_trig2_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-2] == '1'],dtype=int)] + file_prefix
        empty_trig3_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-3] == '1'],dtype=int)] + file_prefix
        empty_trig4_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-4] == '1'],dtype=int)] + file_prefix
        empty_trig5_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-5] == '1'],dtype=int)] + file_prefix
        empty_trig7_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-7] == '1'],dtype=int)] + file_prefix
        empty_trig9_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-9] == '1'],dtype=int)] + file_prefix
        empty_trig10_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-10] == '1'],dtype=int)] + file_prefix
        empty_trig11_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-11] == '1'],dtype=int)] + file_prefix
        empty_trig13_chunk_events = empty_event[np.array([i for i, bit in enumerate(empty_triggerbits) if bit[-13] == '1'],dtype=int)] + file_prefix

        empty_trig1_events = ak.concatenate([empty_trig1_events,empty_trig1_chunk_events])
        empty_trig2_events = ak.concatenate([empty_trig2_events,empty_trig2_chunk_events])
        empty_trig3_events = ak.concatenate([empty_trig3_events,empty_trig3_chunk_events])
        empty_trig4_events = ak.concatenate([empty_trig4_events,empty_trig4_chunk_events])
        empty_trig5_events = ak.concatenate([empty_trig5_events,empty_trig5_chunk_events])
        empty_trig7_events = ak.concatenate([empty_trig7_events,empty_trig7_chunk_events])
        empty_trig9_events = ak.concatenate([empty_trig9_events,empty_trig9_chunk_events])
        empty_trig10_events = ak.concatenate([empty_trig10_events,empty_trig10_chunk_events])
        empty_trig11_events = ak.concatenate([empty_trig11_events,empty_trig11_chunk_events])
        empty_trig13_events = ak.concatenate([empty_trig13_events,empty_trig13_chunk_events])
        
    n_zerobias_chunk_pulses = ak.num(pulses["chan"][np.isin(event,online_trig13_chunk_events - file_prefix)])
    n_trig1_chunk_pulses = ak.num(pulses["chan"][np.isin(event,offline_trig1_chunk_events - file_prefix)])
    n_trig2_chunk_pulses = ak.num(pulses["chan"][np.isin(event,offline_trig2_chunk_events - file_prefix)])

    n_zerobias_pulses = ak.concatenate([n_zerobias_pulses,n_zerobias_chunk_pulses])
    n_trig1_pulses = ak.concatenate([n_trig1_pulses,n_trig1_chunk_pulses])
    n_trig2_pulses = ak.concatenate([n_trig2_pulses,n_trig2_chunk_pulses])

    


def calculations(offline_events,online_events):
    if(len(offline_events) != 0):
        eff = str(round(len(offline_events[np.isin(offline_events,online_events)])/len(offline_events),6))
        if(float(eff) == 0): effunc = 0
        else: effunc = round(float(eff) * np.sqrt( (1/len(offline_events[np.isin(offline_events,online_events)])) + (1/len(offline_events)) ),6)
        onoff_frac = str(round(len(online_events)/len(offline_events),6))
    else: 
        eff = "N/A"
        effunc = "N/A"
        onoff_frac = "N/A"
    return eff, effunc, onoff_frac

eff1, eff1unc, onoff_frac1 = calculations(offline_trig1_events,online_trig1_events)
eff2, eff2unc, onoff_frac2 = calculations(offline_trig2_events,online_trig2_events)
eff3, eff3unc, onoff_frac3 = calculations(offline_trig3_events,online_trig3_events)
eff4, eff4unc, onoff_frac4 = calculations(offline_trig4_events,online_trig4_events)
eff5, eff5unc, onoff_frac5 = calculations(offline_trig5_events,online_trig5_events)
eff7, eff7unc, onoff_frac7 = calculations(offline_trig7_events,online_trig7_events)
eff9, eff9unc, onoff_frac9 = calculations(offline_trig9_events,online_trig9_events)
eff10, eff10unc, onoff_frac10 = calculations(offline_trig10_events,online_trig10_events)
eff11, eff11unc, onoff_frac11 = calculations(offline_trig11_events,online_trig11_events)


def empty_fracs(empty_events,online_events):
    if(len(empty_events) != 0):
        frac_empty = round(len(empty_events) / (len(empty_events) + len(online_events)),5)
    else: frac_empty = "0"
    return frac_empty

frac_t1_empty = empty_fracs(empty_trig1_events,online_trig1_events)
frac_t2_empty = empty_fracs(empty_trig2_events,online_trig2_events)
frac_t3_empty = empty_fracs(empty_trig3_events,online_trig3_events)
frac_t4_empty = empty_fracs(empty_trig4_events,online_trig4_events)
frac_t5_empty = empty_fracs(empty_trig5_events,online_trig5_events)
frac_t7_empty = empty_fracs(empty_trig7_events,online_trig7_events)
frac_t9_empty = empty_fracs(empty_trig9_events,online_trig9_events)
frac_t10_empty = empty_fracs(empty_trig10_events,online_trig10_events)
frac_t11_empty = empty_fracs(empty_trig11_events,online_trig11_events)
frac_t13_empty = empty_fracs(empty_trig13_events,online_trig13_events)



#print("Offline trig1 events ",ak.to_list(offline_trig1_events),"\n")
#print("Online trig1 events",ak.to_list(online_trig1_events),"\n")
#print("Offline trig2 events",ak.to_list(offline_trig2_events),"\n")
#print("Online trig2 events",ak.to_list(online_trig2_events),"\n")

print("Trigger Name".ljust(22),"nOnline".ljust(22),"nOffline".ljust(22),"Offline Efficiency".ljust(22),"Fraction Empty".ljust(22))
print("-"*105)
print("FourLayersHit".ljust(22),str(len(online_trig1_events)).ljust(22),str(len(offline_trig1_events)).ljust(22),(eff1+" +- "+str(eff1unc)).ljust(22),str(frac_t1_empty).ljust(22))
print("threeInaRow".ljust(22),str(len(online_trig2_events)).ljust(22),str(len(offline_trig2_events)).ljust(22),(eff2+" +- "+str(eff2unc)).ljust(22),str(frac_t2_empty))
print("twoSeparatedLayers".ljust(22),str(len(online_trig3_events)).ljust(22),str(len(offline_trig3_events)).ljust(22),(eff3+" +- "+str(eff3unc)).ljust(22),str(frac_t3_empty))
print("twoAdjacentLayers".ljust(22),str(len(online_trig4_events)).ljust(22),str(len(offline_trig4_events)).ljust(22),(eff4+" +- "+str(eff4unc)).ljust(22),str(frac_t4_empty))
print("NLayersHit".ljust(22),str(len(online_trig5_events)).ljust(22),str(len(offline_trig5_events)).ljust(22),(eff5+" +- "+str(eff5unc)).ljust(22),str(frac_t5_empty))
print("gtNHits".ljust(22),str(len(online_trig7_events)).ljust(22),str(len(offline_trig7_events)).ljust(22),(eff7+" +- "+str(eff7unc)).ljust(22),str(frac_t7_empty))
print("topPanels".ljust(22),str(len(online_trig9_events)).ljust(22),str(len(offline_trig9_events)).ljust(22),(eff9+" +- "+str(eff9unc)).ljust(22),str(frac_t9_empty))
print("topPanelsBotBars".ljust(22),str(len(online_trig10_events)).ljust(22),str(len(offline_trig10_events)).ljust(22),(eff10+" +- "+str(eff10unc)).ljust(22),str(frac_t10_empty))
print("Front/BackPanels".ljust(22),str(len(online_trig11_events)).ljust(22),str(len(offline_trig11_events)).ljust(22),(eff11+" +- "+str(eff11unc)).ljust(22),str(frac_t11_empty))
print("Zerobias".ljust(22),str(len(online_trig13_events)).ljust(22),"N/A".ljust(22),"N/A".ljust(22),str(frac_t13_empty))


#print("Offline trig1 events that are not found online ",ak.to_list(offline_trig1_events[np.isin(offline_trig1_events,online_trig1_events,invert=True)]))
#print("Offline trig2 events that are not found online ",ak.to_list(offline_trig2_events[np.isin(offline_trig2_events,online_trig2_events,invert=True)]))
#print("Online trig1 events that are not found offline ",ak.to_list(online_trig1_events[np.isin(online_trig1_events,offline_trig1_events,invert=True)]))
#print("Online trig2 events that are not found offline ",ak.to_list(online_trig2_events[np.isin(online_trig2_events,offline_trig2_events,invert=True)]))

'''
with open("trig1offline_py.txt","w") as outfile:
    for t1offevent in offline_trig1_events:
        outfile.write(f"{t1offevent}\n")

with open("trig1online_py.txt","w") as outfile:
    for t1onevent in online_trig1_events:
        outfile.write(f"{t1onevent}\n")

with open("trig2offline_py.txt","w") as outfile:
    for t2offevent in offline_trig2_events:
        outfile.write(f"{t2offevent}\n")

with open("trig2online_py.txt","w") as outfile:
    for t2onevent in online_trig2_events:
        outfile.write(f"{t2onevent}\n")

with open(f"efficiency_files/run{runNumber[0]}.txt","w") as outfile:
    outfile.write(f"{runNumber[0]} {eff1} {eff1unc} {onoff_frac_1} {frac_t1_empty} {eff2} {eff2unc} {onoff_frac_2} {frac_t2_empty} {frac_t13_empty} \n")
'''

def create_hists(name,title,events,run_number,nbins=5,end_bin=5):
    hist = r.TH1F(name,title,nbins,run_number,run_number + end_bin)
    hist.SetBinContent(1,len(events))
    return hist

def create_overlap_hists(name,title,offline_events,online_events,run_number,nbins=5,end_bin=5):
    overlap = np.isin(offline_events,online_events)
    hist = r.TH1F(name,title,nbins,run_number,run_number + end_bin)
    hist.SetBinContent(1,len(offline_events[overlap]))
    return hist

n_online_hists = [create_hists(f"n_online{i}",f"Number of Online Trig {i} Events in Run {runNumber[0]}",eval(f"online_trig{i}_events"),runNumber[0]) for i in [1,2,3,4,5,7,9,10,11,13]]
n_offline_hists = [create_hists(f"n_offline{i}",f"Number of Offline Trig {i} Events in Run {runNumber[0]}",eval(f"offline_trig{i}_events"),runNumber[0]) for i in [1,2,3,4,5,7,9,10,11]]
n_empty_hists = [create_hists(f"n_empty{i}",f"Number of Empty Trig {i} Events in Run {runNumber[0]}", eval(f"empty_trig{i}_events"), runNumber[0]) for i in [1,2,3,4,5,7,9,10,11,13]]
n_overlap_hists = [create_overlap_hists(f"n_overlap{i}",f"Number of Overlap Trig {i} Events in Run {runNumber[0]}",eval(f"offline_trig{i}_events"),eval(f"online_trig{i}_events"),runNumber[0]) for i in [1,2,3,4,5,7,9,10,11]]
ZBpulses.fill(n_zerobias_pulses)
fourLayersPulses.fill(n_trig1_pulses)
threeInRowPulses.fill(n_trig2_pulses)

with uproot.recreate(f"output_run_{runNumber[0]}_file_{fileNumber[0]}.root") as outfile:
    outfile["channels1"] = channels1
    outfile["channels1_and_2"] = channels1_and_2
    outfile["channels1_not_2"] = channels1_not_2
    outfile["ZBpulses"] = ZBpulses
    outfile["fourLayersPulses"] = fourLayersPulses
    outfile["threeInRowPulses"] = threeInRowPulses

    for i in [1,2,3,4,5,7,9,10,11,13]:
        if i <= 5: j = i-1
        elif i == 7: j = i-2
        elif i >= 9 and i < 13: j = i-3
        elif i == 13: j = i-4
        outfile[f"n_online{i}"] = n_online_hists[j]
        outfile[f"n_empty{i}"] = n_empty_hists[j]
        if i == 13: continue
        outfile[f"n_offline{i}"] = n_offline_hists[j]
        outfile[f"n_overlap{i}"] = n_overlap_hists[j] 





#channels1.plot()
#plt.show()




