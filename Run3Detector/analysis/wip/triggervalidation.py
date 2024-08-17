import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import vector
import hist

#Open root file and ttree
file = uproot.open("MilliQan_Run1415.1_v34.root")
tree = file["t;1"]
stop = 1000   #Set the number of files to run on



branches = tree.arrays(["time","height","area","row","column","layer","duration","chan","type"],entry_stop=stop)

#Fix mislabeled channels
chan_mask78 = (branches["chan"] == 78)
chan_mask79 = (branches["chan"] == 79)
new_chan = ak.where(chan_mask78, 24, branches["chan"])
new_chan = ak.where(chan_mask79, 25, new_chan)
branches = ak.with_field(branches, new_chan, "chan")

tTrigger = tree["tTrigger"].array(entry_stop=stop)
event = tree["event"].array(entry_stop=stop)

matched_mask = tTrigger != -1  #Require matched triggers

branches = branches[matched_mask]
tTrigger = tTrigger[matched_mask]
event = event[matched_mask]

#Turn the decimal tTrigger branch into an array of bitstrings
bin_rep_vec = np.vectorize(np.binary_repr)
trig_np = ak.to_numpy(tTrigger).astype(int)
triggerbits = bin_rep_vec(trig_np,width=13)



dynamicPedestal = tree["dynamicPedestal"].array(entry_stop=stop)
dynamicPedestal = dynamicPedestal[matched_mask]

online_height = branches["height"] + dynamicPedestal[branches["chan"]]
heightmask = online_height > 15
height_cut = branches[heightmask]  #Apply dynamic pedestal correction


#Now zip all these pulse shaped branches together into a record called pulses
pulses = ak.zip(
    {
        "time": height_cut["time"],
        "height": height_cut["height"],
        "area": height_cut["area"],
        "chan": height_cut["chan"],
        "row": height_cut["row"],
        "column": height_cut["column"],
        "layer": height_cut["layer"],
        "duration": height_cut["duration"],
        "type": height_cut["type"],
    }
)

#Now for timing combinations. We use ak.combinations to do N choose 2 on the existing pulses to form pairs of pulses.
pairs = ak.combinations(pulses,2)

#Unzip the pairs to get each half of the combinations.
pulse1, pulse2 = ak.unzip(pairs)

#Form a mask for pulses within 150ns of each other. This mask will cut down on the combinatorics.
window = 150
timemask = np.abs(pulse1.time - pulse2.time) < window

#Useful for later
layerdiff = np.abs(pulse1.layer - pulse2.layer)

#============================Four layers offline========================================
#Form a mask to require that the paired pulses are not in the same layer and aren't panels. We don't care about those combinations anyway and removing them speeds things up.
layermask = layerdiff != 0   #pulse1.layer != pulse2.layer
not_panels = (pulse1["type"] == 0) & (pulse2["type"] == 0) 

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


#============================Three in row offline========================================
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
trig2cand2 = time_124 & diff_layer124 & same_row123 & same_cols124
trig2cand3 = time_234 & diff_layer234 & same_row234 & same_cols234
trig2cand4 = time_134 & diff_layer134 & same_row134 & same_cols134
threeInRow_mask = trig2cand1 | trig2cand2 | trig2cand3 | trig2cand4

offline_trig2_events = event[ak.any(threeInRow_mask,axis=1)]

#============================Two separated layers========================================
nonadjacent = layerdiff > 1
twoSeparatedLayers_mask = nonadjacent & not_panels & timemask
offline_trig3_events = event[ak.any(twoSeparatedLayers_mask,axis=1)]

#============================Two adjacent layers=========================================
adjacent = layerdiff == 1
twoAdjacentLayers_mask = adjacent & not_panels & timemask
offline_trig4_events = event[ak.any(twoAdjacentLayers_mask,axis=1)]

#============================N Layers Hit================================================
nLayers = np.array([len(np.unique(row)) for row in ak.to_list(pulses.layer[pulses["type"] == 0 ])])
nLayers_mask = nLayers >= 3
offline_trig5_events = event[nLayers_mask]

#============================Greater than N Hits=========================================
numPulses = ak.num(pulses)
nHit = 3
gtNHits_mask = numPulses >= nHit + 1
offline_trig7_events = event[gtNHits_mask]

#============================Top Panels==================================================
topPanels_mask = ak.any(pulses.row == 4,axis=1)
offline_trig9_events = event[topPanels_mask]

#============================Top Panels + Bottom Bars====================================
topPanelsBotBars_mask = ak.any(pulses.row == 4,axis=1) & ak.any( (pulses.row == 0) & (pulses["type"] == 0),axis=1)
offline_trig10_events = event[topPanelsBotBars_mask]

#============================Front/Back Panels==================================================
frontBack_mask = ak.any(pulses.layer == -1, axis=1) & ak.any(pulses.layer == 4, axis=1)
offline_trig11_events = event[frontBack_mask]



#Online trigger events
online_trig1_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-1] == '1'])]
online_trig2_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-2] == '1'])]
online_trig3_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-3] == '1'])]
online_trig4_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-4] == '1'])]
online_trig5_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-5] == '1'])]
#online_trig6_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-6] == '1'])]
online_trig7_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-7] == '1'])]
#online_trig8_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-8] == '1'])]
online_trig9_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-9] == '1'])]
online_trig10_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-10] == '1'])]
online_trig11_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-11] == '1'])]
#online_trig12_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-12] == '1'])]
online_trig13_events = event[np.array([i for i, bit in enumerate(triggerbits) if bit[-13] == '1'])]




#print("Offline trig1 events ",ak.to_list(offline_trig1_events),"\n")
#print("Online trig1 events",ak.to_list(online_trig1_events),"\n")
#print("Offline trig2 events",ak.to_list(offline_trig2_events),"\n")
#print("Online trig2 events",ak.to_list(online_trig2_events),"\n")

print("Number of offline trig1 events ",len(offline_trig1_events))
print("Number of online trig1 events ",len(online_trig1_events))
print("Number of offline trig2 events",len(offline_trig2_events))
print("Number of online trig2 events ",len(online_trig2_events))
print("Number of offline trig3 events ",len(offline_trig3_events))
print("Number of online trig3 events ",len(online_trig3_events))
print("Number of offline trig4 events",len(offline_trig4_events))
print("Number of online trig4 events ",len(online_trig4_events))
print("Number of offline trig5 events ",len(offline_trig5_events))
print("Number of online trig5 events ",len(online_trig5_events))
print("Number of offline trig7 events",len(offline_trig7_events))
print("Number of online trig7 events ",len(online_trig7_events))
print("Number of offline trig9 events",len(offline_trig9_events))
print("Number of online trig9 events ",len(online_trig9_events))
print("Number of offline trig10 events",len(offline_trig10_events))
print("Number of online trig10 events ",len(online_trig10_events))
print("Number of offline trig11 events",len(offline_trig11_events))
print("Number of online trig11 events ",len(online_trig11_events),"\n")
print("Offline trig1 events that are not found online ",ak.to_list(offline_trig1_events[np.isin(offline_trig1_events,online_trig1_events,invert=True)]))
print("Offline trig2 events that are not found online ",ak.to_list(offline_trig2_events[np.isin(offline_trig2_events,online_trig2_events,invert=True)]))
print("Online trig1 events that are not found offline ",ak.to_list(online_trig1_events[np.isin(online_trig1_events,offline_trig1_events,invert=True)]))
print("Online trig2 events that are not found offline ",ak.to_list(online_trig2_events[np.isin(online_trig2_events,offline_trig2_events,invert=True)]))

