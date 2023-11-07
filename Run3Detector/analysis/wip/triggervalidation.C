#include <iostream>
#include <fstream>
#include <string>
#include "TSystem.h"
#include <vector>
#include <bitset>
using namespace std;

TFile* OpenRootFile(const char* filename) {
    TFile* file = new TFile(filename, "READ");

    if (file->IsZombie()) {
        std::cerr << "Error: Failed to open the ROOT file." << std::endl;
        return nullptr;
    }

    return file;

}


void triggervalidation(const char* filename)
{
	//Getting files
	//TFile *f1 = new TFile("/net/cms26/cms26r0/milliqan/Run3Offline/v32/MilliQan_Run1116.1_default_v32.root"); 
	//TFile *f1 = new TFile("/net/cms26/cms26r0/milliqan/Run3Offline/v33/MilliQan_Run1172.100_v33_firstPedestals.root");
	//TFile *f1 = new TFile("/net/cms26/cms26r0/milliqan/Run3Offline/v33/MilliQan_Run1172_combined.root");
	//TFile *f1 = new TFile("/net/cms2/cms2r0/neha/hadd_files/MilliQan_Run1173_combined.root");
	//TFile *f1 = new TFile("/net/cms26/cms26r0/milliqan/Run3Offline/v32/MilliQan_Run1172.1_default_v32.root");
	//TFile *f1 = new TFile("MilliQan_Run1114.1_matched.root");

	gROOT->Reset();
	TFile *f1 = OpenRootFile(filename);

	//Getting trees
	TTree *t1 = (TTree*)f1->Get("t");
	const Int_t nentries1 = t1->GetEntries();

	std::cout << nentries1 << std::endl;

	//Declaring tree variables
	int event; 
	Long64_t runNumber;
	Long64_t fileNumber;
	float tTrigger;
	std::vector<int> *layer = 0;	
	std::vector<int> *chan = 0;
	std::vector<int> *row = 0;
	std::vector<int> *column = 0;
	std::vector<float> *dynamicPedestal = 0;
	std::vector<float> *height = 0;
	std::vector<float> *area = 0;
	std::vector<float> *time = 0;
	std::vector<Long64_t> *v_groupTDC_g0 = 0;
	std::vector<float> *duration = 0;
	
	std::vector<std::vector<int>> overlappingHits;
    	std::vector<int> vec;
	std::vector<int> vec2;
	std::vector<int> trig1_coincidence_satisfying_comb;
	std::vector<std::vector<int>> trig2_coincidence_satisfying_comb;

	//Setting branch addresses
	t1->SetBranchAddress("runNumber",&runNumber);
	t1->SetBranchAddress("fileNumber",&fileNumber);
	t1->SetBranchAddress("event",&event);
	t1->SetBranchAddress("layer",&layer);
	t1->SetBranchAddress("chan",&chan);
	t1->SetBranchAddress("tTrigger",&tTrigger);
	t1->SetBranchAddress("row",&row);
	t1->SetBranchAddress("column",&column);
	t1->SetBranchAddress("dynamicPedestal",&dynamicPedestal);
	t1->SetBranchAddress("height",&height);
	t1->SetBranchAddress("area",&area);	
	t1->SetBranchAddress("time",&time);
	t1->SetBranchAddress("v_groupTDC_g0",&v_groupTDC_g0);
	t1->SetBranchAddress("duration",&duration);

	//Offline triggers
	bool fourLayersHit = false;
	bool threeInRow = false; 
	bool twoSeparatedLayers = false;
	bool twoAdjacentLayers = false;
	bool NLayersHit = false;
	bool gtNHits = false;
	bool topPanels = false;
	bool topPanels_plus_BottomBars = false;
	bool front_back_panels = false;


	//Online triggers
	bool trig1 = false;
	bool trig2 = false;
	bool trig3 = false;
	bool trig4 = false;
	bool trig5 = false;
	bool trig6 = false;
	bool trig7 = false;
	bool trig8 = false;
	bool trig9 = false;
	bool trig10 = false;
	bool trig11 = false;
	bool trig12 = false;
	bool trig13 = false;

	//Other variables of interest
	bool layer0hit = false;
	bool layer1hit = false;
	bool layer2hit = false;
	bool layer3hit = false;
	bool threeLayersHit = false;
	int nHits;
    
	float vMax_online = 0;
	int nLayerThreshold = 3;   //Obtained from trigger conig file: https://mcitron-public.web.cern.ch/milliqanRunLog/configs/Run1116TriggerDefault.py
	int nHitThreshold = 3;
	int layerCount = 0;
	int hitCount = 0;
	double prescale1 = 1;
	double prescale2 = 1;
	double prescale3 = 0.0001;
	double prescale4 = 0.0001;
	double prescale5 = 0.001;
	double prescale6 = 1;
	double prescale7 = 0.0001;
	double prescale8 = 1;
	double prescale9 = 0.0001;
	double prescale10 = 0.01;
	double prescale11 = 1;
	double prescale12 = 0;
	double prescale13 = 1;
	double frac1; double frac2; double frac3; double frac4; double frac5;
	double frac6; double frac7; double frac8; double frac9; double frac10;
	double frac11; double frac12; double frac13;
	double lastPulseStart = 0; double firstPulseStart = 0;
	double lastPulseStartCorr = 0; double firstPulseStartCorr = 0;
	double pulse_i_end = 0; double pulse_j_end = 0; double pulse_k_end = 0; double pulse_l_end = 0;
	double deltaT = 0; double window = 0 ;
	double firstPulseTime = 100000; double lastPulseTime = 0; 
	double deltaT_check = 0;
	bool isFirstTrig1found = true; bool isFirstTrig2found = true;

	int nTrig1_online = 0;int nTrig2_online = 0;int nTrig3_online = 0;int nTrig4_online = 0;int nTrig5_online = 0;int nTrig6_online = 0;
	int nTrig7_online = 0;int nTrig8_online = 0;int nTrig9_online = 0;int nTrig10_online = 0;int nTrig11_online = 0;int nTrig12_online = 0;int nTrig13_online = 0;

	int nTrig1_offline = 0;int nTrig2_offline = 0;int nTrig3_offline = 0;int nTrig4_offline = 0;int nTrig5_offline = 0;int nTrig6_offline = 0;
	int nTrig7_offline = 0;int nTrig8_offline = 0;int nTrig9_offline = 0;int nTrig10_offline = 0;int nTrig11_offline = 0;int nTrig12_offline = 0;int nTrig13_offline = 0;
	int nZeroBias_trig1off = 0; int nZeroBias_trig1on = 0; int nZeroBias_trig2off = 0; int nZeroBias_trig2on = 0;


	//Online/offline histograms
	TH1I *h_trigger_rates_on = new TH1I("h_trigger_rates_on","Number of Times Triggers Occur; Trigger Number; Number of Events",13,0,13);  //online
	TH1I *h_trigger_rates_off = new TH1I("h_trigger_rates_off","Number of Times Triggers Occur; Trigger Number; Number of Events",13,0,13); //offline
	TH1F *h_ctrl_deltaT = new TH1F("h_ctrl_deltaT","DeltaT of matched combinations; DeltaT [ns]; Number of Events", 120,0,500);
	TH1F *h_mismatched_deltaT = new TH1F("h_mismatched_deltaT","DeltaT of mismatched combinations; DeltaT [ns]; Number of Events", 120,0,500);
	TH1F *h_ctrl_vMax = new TH1F("h_ctrl_vMax","vMax of matched combinations; vMax [mV]; Number of Events", 80,0,1300);
	TH1F *h_mismatched_vMax = new TH1F("h_mismatched_vMax","vMax of mismatched combinations; vMax [mV]; Number of Events", 80,0,1300);
	TH1F *h_ctrl_area = new TH1F("h_ctrl_area","Pulse area of matched combinations; area [pV s]; Number of Events", 80,0,500000);
	TH1F *h_mismatched_area = new TH1F("h_mismatched_area","Pulse area of mismatched combinations; area [pV s]; Number of Events", 80,0,500000);
	TH1I *h_ctrl_chan = new TH1I("h_ctrl_chan","chan of matched combinations; channel; Number of Events",80,0,80);
	TH1I *h_mismatched_chan = new TH1I("h_mismatched_chan","chan of mismatched combinations; channel; Number of Events",80,0,80);
	TH1F *h_ctrl_startTime = new TH1F("h_ctrl_startTime","Starting time of matched combinations; time [ns]; Number of Events", 120, 0, 2600);
	TH1F *h_mismatched_startTime = new TH1F("h_mismatched_startTime","Starting time of mismatched combinations; time [ns]; Number of Events", 120, 0, 2600);
	TH2F *h_ctrl_height_v_area = new TH2F("h_ctrl_height_v_area","Height vs area for matched combinations; area [pV s]; height [mV]",80,0,30000,80,0,1300);
	TH2F *h_mismatched_height_v_area = new TH2F("h_mismatched_height_v_area","Height vs area for mismatched combinations; area [pV s]; height [mV]",80,0,30000,80,0,1300);


	h_trigger_rates_on->SetName("Online");
	h_trigger_rates_on->SetName("Offline");

    	std::ofstream onlineFile1("onlineTrig1Evts.txt");
	onlineFile1 << "Run " << runNumber << ", file " << fileNumber << ", trigger 1 events" << std::endl;
	std::ofstream offlineFile1("offlineTrig1Evts.txt");
	offlineFile1 << "Run " << runNumber << ", file " << fileNumber << " trigger 1 events" << std::endl;
	std::ofstream onlineFile2("onlineTrig2Evts.txt");
	onlineFile2 << "Run " << runNumber << ", file " << fileNumber << " trigger 2 events" << std::endl;
	std::ofstream offlineFile2("offlineTrig2Evts.txt");
	offlineFile2 << "Run " << runNumber << ", file " << fileNumber << " trigger 2 events" << std::endl;

	std::ofstream timingStudy("timingStudy.txt");
	timingStudy <<  "Run " << runNumber << ", file " << fileNumber << std::endl;
	timingStudy << "Events with trig1 failing online, trig2 passing online, both pass offline" << std::endl;
	timingStudy << "Event firstChan firstTime lastChan lastTime deltaT lastChan+dur lastTime+dur deltaT+dur" << std::endl;


    //Begin loop over number of events
    for (Int_t ievent =0;ievent<nentries1;ievent++)      
    //for(Int_t ievent = 0; ievent < 10; ievent++)
    {
	
    	t1->GetEntry(ievent);
    	if(chan == nullptr || chan->empty())
    	{
    		continue;
    	}
    	if(layer == nullptr || layer->empty())
    	{
    		continue;
    	}

    	//Setting variables for this event
    	//Offline triggers
    	fourLayersHit = false;
    	threeInRow = false; 
    	twoSeparatedLayers = false;
    	twoAdjacentLayers = false;
    	NLayersHit = false;
    	gtNHits = false;
    	topPanels = false;
    	topPanels_plus_BottomBars = false;
    	front_back_panels = false;


    	//Online triggers
    	trig1 = false;
    	trig2 = false;
    	trig3 = false;
    	trig4 = false;
    	trig5 = false;
    	trig6 = false;
    	trig7 = false;
    	trig8 = false;
    	trig9 = false;
    	trig10 = false;

    	trig11 = false;
    	trig12 = false;
    	trig13 = false;

    	//Other variables of interest
        layer0hit = false;
	layer1hit = false;
	layer2hit = false;
	layer3hit = false;
	threeLayersHit = false;
    	nHits = layer->size();
        bool goodHits[nHits];
    	for(int i = 0; i < nHits; i++) goodHits[i] = false;
    	layerCount = 0;
    	hitCount = 0;
    	lastPulseStart = 0; firstPulseStart = 0;
    	pulse_i_end = 0; pulse_j_end = 0; pulse_k_end = 0; pulse_l_end = 0;
    	deltaT = 0; window = 0;
    	overlappingHits.clear();
	trig1_coincidence_satisfying_comb.clear();
	trig2_coincidence_satisfying_comb.clear();
	
        std::bitset<13> binaryTrigger(tTrigger);
        for(int i = 0; i < 13; i++) //Getting online triggers
        {
            if(binaryTrigger[i])
            {
                h_trigger_rates_on->Fill(i+1);
                
                if(i+1 == 1) {trig1 = true; nTrig1_online++; onlineFile1  << event << std::endl;}
                if(i+1 == 2) {trig2 = true;nTrig2_online++; onlineFile2   << event << std::endl;}
                if(i+1 == 3) {trig3 = true;nTrig3_online++;}
                if(i+1 == 4) {trig4 = true;nTrig4_online++;}
                if(i+1 == 5) {trig5 = true;nTrig5_online++;}
                if(i+1 == 6) {trig6 = true;nTrig6_online++;}
                if(i+1 == 7) {trig7 = true;nTrig7_online++;}
                if(i+1 == 8) {trig8 = true;nTrig8_online++;}
                if(i+1 == 9) {trig9 = true;nTrig9_online++;}
                if(i+1 == 10) {trig10 = true;nTrig10_online++;}
                if(i+1 == 11) {trig11 = true;nTrig11_online++;}
                if(i+1 == 12) {trig12 = true;nTrig12_online++;}
                if(i+1 == 13) {trig13 = true;nTrig13_online++;}
                
            }
        }

        

    	//====Marking good channels for the entire event====
    	for(size_t j = 0; j < nHits; ++j)
	{

		//Takes channel mismapping into account. Can remove once fixed.
		if((*chan)[j] == 78) {(*chan)[j] = 24; (*layer)[j] = 1;(*row)[j] = 1;(*column)[j] = 0;}  
		if((*chan)[j] == 79) {(*chan)[j] = 25; (*layer)[j] = 1;(*row)[j] = 1;(*column)[j] = 1;}

		
		//Ensure that the online pulse height actually passed trigger requirement and it is a bar
		vMax_online = (*height)[j] + (*dynamicPedestal)[(*chan)[j]];
		//std::cout << "chan " << (*chan)[j] << " time " << (*time)[j] << "vMax_online " << vMax_online << std::endl;
		//vMax_online = (*height)[j]
		if(vMax_online >= 15 && (*chan)[j] <= 63 && (*chan)[j] >= 0) goodHits[j] = true;       
		//if((*chan)[j] <= 63 && (*chan)[j] >= 0) goodHits[j] = true; 
		
	} 

	
	
	//std::cout << "Event " << event << std::endl;
	
        //Finding pulses in coincidence window
        for(int i = 0; i < nHits; ++i)
        {
	    vec.clear();
	    if(goodHits[i] == 0) {overlappingHits.push_back(vec); continue;}
	    vec.push_back(i);
	    firstPulseStart = (*time)[i];
            firstPulseStartCorr = (*time)[i] + 5*((*v_groupTDC_g0)[int((*chan)[i]/16)]-(*v_groupTDC_g0)[0]);
            window = 348; //100ns from first pulse plus 160ns coincident time window            
	    /*
            std::cout<<" ---------New combination check -------" <<std::endl;
            std::cout << "chan " << (*chan)[i] << std::endl;
            std::cout << "time " << (*time)[i] << std::endl;
            std::cout << "height " << (*height)[i] << std::endl;
            std::cout << "duration " << (*duration)[i] << std::endl;
            std::cout << "window " << window << std::endl;
            std::cout << "goodHit: " << goodHits[i] << std::endl;
	    std::cout << " -- Possible combinations ---"<<endl;
            */
            for(int j = 0; j < nHits; ++j)
            {
                if(i==j) continue;

                //Check whether pulse satisfies timing requirement
                //pulse_i_end = (*time)[i] + (*duration)[i];
                //pulse_j_end = (*time)[j] + (*duration)[j];
                
                lastPulseStart = (*time)[j];
                lastPulseStartCorr = (*time)[j] + 5*((*v_groupTDC_g0)[int((*chan)[j]/16)]-(*v_groupTDC_g0)[0]);
                //deltaT = lastPulseStartCorr - firstPulseStartCorr;
		deltaT = lastPulseStart - firstPulseStart;
                /*
                std::cout << (*chan)[j] << std::endl;
                std::cout << (*time)[j] << std::endl;
                std::cout << (*height)[j] << std::endl;
                std::cout << (*duration)[j] << std::endl;
		std::cout << goodHits[j] << std::endl;
		std::cout << "deltaT " << deltaT << std::endl;
                */
                
		
                if(deltaT < 0 || deltaT > window) continue;
		if(goodHits[j]) vec.push_back(j);
            }
            overlappingHits.push_back(vec);
        }

	


    	
    	
    	for(int i = 0; i < overlappingHits.size(); ++i) //Looping over combinations of hits that fall in 288ns window
        {
            layer0hit = false;
            layer1hit = false;
            layer2hit = false;
            layer3hit = false;
            layerCount = 0;
            isFirstTrig1found = true;
            isFirstTrig2found = true;
	    vec2.clear();
	    //std::cout << std::endl;
            //std::cout << "time " << (*time)[i] << std::endl;
	    //std::cout << "overlappingHits " << "goodHit" << std::endl;
	    for(int j = 0; j < overlappingHits[i].size(); ++j)  //Looping over hits within a given combo that fall within 288ns window
            {

                //std::cout << (*chan)[overlappingHits[i][j]] << " " << goodHits[overlappingHits[i][j]] << std::endl;

                //==========Checking first trigger logic===========
                if((*layer)[overlappingHits[i][j]] == 0 && goodHits[overlappingHits[i][j]]) layer0hit = true;
                if((*layer)[overlappingHits[i][j]] == 1 && goodHits[overlappingHits[i][j]]) layer1hit = true;   
                if((*layer)[overlappingHits[i][j]] == 2 && goodHits[overlappingHits[i][j]]) layer2hit = true;
                if((*layer)[overlappingHits[i][j]] == 3 && goodHits[overlappingHits[i][j]]) layer3hit = true;
                if(layer0hit && layer1hit && layer2hit && layer3hit) fourLayersHit = true; 
		if(layer0hit && layer1hit && layer2hit && layer3hit && isFirstTrig1found == true) {trig1_coincidence_satisfying_comb.push_back(i); isFirstTrig1found = false;}

                //==========Checking second trigger logic===========
                if( (layer0hit && layer1hit && layer2hit) || (layer1hit && layer2hit && layer3hit) || (layer0hit && layer2hit && layer3hit) || (layer0hit && layer1hit && layer3hit)) threeLayersHit = true;
                //Looping over two other hits
                if(threeLayersHit)
                {
                    for(int k = 0; k < overlappingHits[i].size(); ++k)
                    {
                        for(int l = 0; l < overlappingHits[i].size(); ++l)
                        {
                            //Not in same layer
                            if((*layer)[overlappingHits[i][j]] == (*layer)[overlappingHits[i][k]] || (*layer)[overlappingHits[i][k]] == (*layer)[overlappingHits[i][l]] || (*layer)[overlappingHits[i][l]] == (*layer)[overlappingHits[i][j]]) continue; 

                            if(goodHits[overlappingHits[i][j]] && goodHits[overlappingHits[i][k]] && goodHits[overlappingHits[i][l]])
                            {
                                //Check to see if they are in the same trigger group
                                if( (((*column)[overlappingHits[i][j]]/2) == (*column)[overlappingHits[i][k]]/2) && ((*row)[overlappingHits[i][j]] == (*row)[overlappingHits[i][k]]))
                                {
                                    if( (((*column)[overlappingHits[i][k]]/2) == (*column)[overlappingHits[i][l]]/2) && ((*row)[overlappingHits[i][k]] == (*row)[overlappingHits[i][l]]))
                                    {
                                        threeInRow = true; //All three in same group
					if(isFirstTrig2found == true) {
					vec2.push_back(overlappingHits[i][j]);
					vec2.push_back(overlappingHits[i][j]); //Pushing back pulse number and combination
					vec2.push_back(overlappingHits[i][l]); 
					trig2_coincidence_satisfying_comb.push_back(vec2);			
					isFirstTrig2found = false;
					
					}
                                    }
                                }
                            }

                        }
                    }
                }

                //==========Checking third trigger logic===========
                if( (layer0hit && layer2hit) || (layer0hit && layer3hit) || ((layer1hit && layer3hit)) ) twoSeparatedLayers = true;

                //==========Checking fourth trigger logic==========
                if( (layer0hit && layer1hit) || (layer1hit && layer2hit) || ((layer2hit && layer3hit)) ) twoAdjacentLayers = true;

                //==========Checking fifth trigger logic===========
                layerCount = layer0hit + layer1hit + layer2hit + layer3hit;
                if(layerCount >= nLayerThreshold) NLayersHit = true;

                //==========Checking seventh trigger logic=========
                hitCount = nHits;
                if(hitCount >= (nHitThreshold + 1)) gtNHits = true;

                //==========Checking ninth trigger logic===========
                if((*chan)[overlappingHits[i][j]] == 68 || (*chan)[overlappingHits[i][j]] == 72) topPanels = true;

                //==========Checking tenth trigger logic===========             
                if((*chan)[overlappingHits[i][j]] == 68 || (*chan)[overlappingHits[i][j]] == 72) //If a top panel is hit...
                {
                    for(int k = 0; k < overlappingHits[i].size(); ++k)    //check to see whether a bar in the bottom row is hit as well
                    {
                        if((*row)[overlappingHits[i][k]] == 0) topPanels_plus_BottomBars = true;
                    }
                }

                //==========Checking eleventh trigger logic========
                if((*chan)[overlappingHits[i][j]] == 74)
                {
                    for(int k = 0; k < overlappingHits[i].size(); ++k)
                    {
                        if((*chan)[overlappingHits[i][k]] == 75) front_back_panels = true;
                    }
                }
            }


	}
            


        
	//Fill histograms
    	if(fourLayersHit) {nTrig1_offline++; h_trigger_rates_off->Fill(1); offlineFile1 << event << std::endl;}
        if(threeInRow) {nTrig2_offline++; h_trigger_rates_off->Fill(2); offlineFile2 << event << std::endl;}
        if(twoSeparatedLayers) {nTrig3_offline++;h_trigger_rates_off->Fill(3);}
        if(twoAdjacentLayers) {nTrig4_offline++;h_trigger_rates_off->Fill(4);}
        if(NLayersHit) {nTrig5_offline++; h_trigger_rates_off->Fill(5);}
        if(gtNHits) {nTrig7_offline++; h_trigger_rates_off->Fill(7);}
        if(topPanels) {nTrig9_offline++;h_trigger_rates_off->Fill(9);}
        if(topPanels_plus_BottomBars) {nTrig10_offline++;h_trigger_rates_off->Fill(10);}
        if(front_back_panels) {nTrig11_offline++;h_trigger_rates_off->Fill(11);}
	
	
	//============Begin comparison of matched vs mismatched events
	//
	if(trig1 && !trig2 && fourLayersHit && !threeInRow) //Control group: matched trig1 only
	{
		//std::cout << "Event " << event << std::endl;
		for(int i = 0; i < trig1_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig1
		{
			lastPulseTime = 0; firstPulseTime = 100000;
			//std::cout << " " << std::endl;
			for(int j = 0; j < overlappingHits[trig1_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
			{
				//std::cout << "Chan: " << (*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << std::endl;
				if((*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] < firstPulseTime) 
				{
					firstPulseTime = (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]];
				}
				if((*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] > lastPulseTime) 
                                {
                                        lastPulseTime = (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]];
                                }
				h_ctrl_vMax->Fill((*height)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
                        	h_ctrl_area->Fill((*area)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
				h_ctrl_chan->Fill((*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
				h_ctrl_height_v_area->Fill((*area)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]], (*height)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
			}
			deltaT_check = lastPulseTime - firstPulseTime;
			//std::cout << "lastPulseTime: " << lastPulseTime << " firstPulseTime: " << firstPulseTime << " deltaT: " << deltaT_check << std::endl;
			h_ctrl_deltaT->Fill(deltaT_check);
			h_ctrl_startTime->Fill(firstPulseTime);
		}
	}
	

	if(!trig1 && trig2 && !fourLayersHit && threeInRow) //Control group: matched trig2 only
        {       
                //std::cout << "Event " << event << std::endl;
                for(int i = 0; i < trig2_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig2
                {       
                        lastPulseTime = 0; firstPulseTime = 100000;
                        //std::cout << " " << std::endl;
                        
			for(int j = 0; j < trig2_coincidence_satisfying_comb[i].size(); j++) //Looping over j pulses in ith coincident window
                        {       
                                //std::cout << "Chan: " << (*chan)[trig2_coincidence_satisfying_comb[i][j]] << " Time: " << (*time)[trig2_coincidence_satisfying_comb[i][j]] << std::endl;
                                if((*time)[trig2_coincidence_satisfying_comb[i][j]] < firstPulseTime)
                                {       
                                        firstPulseTime = (*time)[trig2_coincidence_satisfying_comb[i][j]];
                                }
                                if((*time)[trig2_coincidence_satisfying_comb[i][j]] > lastPulseTime)
                                {       
                                        lastPulseTime = (*time)[trig2_coincidence_satisfying_comb[i][j]];
                                }
				h_ctrl_vMax->Fill((*height)[trig2_coincidence_satisfying_comb[i][j]]);
                        	h_ctrl_area->Fill((*area)[trig2_coincidence_satisfying_comb[i][j]]);
				h_ctrl_chan->Fill((*chan)[trig2_coincidence_satisfying_comb[i][j]]);
                        	h_ctrl_height_v_area->Fill((*area)[trig2_coincidence_satisfying_comb[i][j]], (*height)[trig2_coincidence_satisfying_comb[i][j]]);
			} 
                        deltaT_check = lastPulseTime - firstPulseTime;
                        
			//std::cout << "lastPulseTime: " << lastPulseTime << " firstPulseTime: " << firstPulseTime << " deltaT: " << deltaT_check << std::endl;
                        h_ctrl_deltaT->Fill(deltaT_check);
			h_ctrl_startTime->Fill(firstPulseTime);
                }
        }
	

	if(!trig1 && trig2 && fourLayersHit && threeInRow) //Mismatched group
        {
		for(int i = 0; i < trig1_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig1
                {       
                        lastPulseTime = 0; firstPulseTime = 100000;
                        //std::cout << " " << std::endl;
                        for(int j = 0; j < overlappingHits[trig1_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
                        {       
                                //std::cout << "Chan: " << (*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << std::endl;
                                if((*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] < firstPulseTime)
                                {       
                                        firstPulseTime = (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]];
                                }
                                if((*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] > lastPulseTime)
                                {       
                                        lastPulseTime = (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]];
                                }
                        	h_mismatched_vMax->Fill((*height)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
                        	h_mismatched_area->Fill((*area)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
				h_mismatched_chan->Fill((*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
                        	h_mismatched_height_v_area->Fill((*area)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]], (*height)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]]);
			}
                        deltaT_check = lastPulseTime - firstPulseTime;
                        //std::cout << "lastPulseTime: " << lastPulseTime << " firstPulseTime: " << firstPulseTime << " deltaT: " << deltaT_check << std::endl;
                        h_mismatched_deltaT->Fill(deltaT_check);
			h_mismatched_startTime->Fill(firstPulseTime);
                }
		
                //std::cout << "Event " << event << std::endl;
                for(int i = 0; i < trig2_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig2
                {
                        lastPulseTime = 0; firstPulseTime = 100000;
                        //std::cout << " " << std::endl;
                        for(int j = 0; j < trig2_coincidence_satisfying_comb[i].size(); j++) //Looping over pulses in coincident window
                        {
                                //std::cout << "Chan: " << (*chan)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << std::endl;
                                if((*time)[trig2_coincidence_satisfying_comb[i][j]] < firstPulseTime)
                                {
                                        firstPulseTime = (*time)[trig2_coincidence_satisfying_comb[i][j]];
                                }
                                if((*time)[trig2_coincidence_satisfying_comb[i][j]] > lastPulseTime)
                                {
                                        lastPulseTime = (*time)[trig2_coincidence_satisfying_comb[i][j]];
                                }
				h_mismatched_vMax->Fill((*height)[trig2_coincidence_satisfying_comb[i][j]]);
                        	h_mismatched_area->Fill((*area)[trig2_coincidence_satisfying_comb[i][j]]);
				h_mismatched_chan->Fill((*chan)[trig2_coincidence_satisfying_comb[i][j]]);
                        	h_mismatched_height_v_area->Fill((*area)[trig2_coincidence_satisfying_comb[i][j]], (*height)[trig2_coincidence_satisfying_comb[i][j]]);
			}
                        deltaT_check = lastPulseTime - firstPulseTime;
                        //std::cout << "lastPulseTime: " << lastPulseTime << " firstPulseTime: " << firstPulseTime << " deltaT: " << deltaT_check << std::endl;
                        h_mismatched_deltaT->Fill(deltaT_check);
			h_mismatched_startTime->Fill(firstPulseTime);
                }
        }

	//Zerobias study
	if(trig13 && fourLayersHit){nZeroBias_trig1off++;}
	if(trig13 && trig1){nZeroBias_trig1on++;}
	if(trig13 && threeInRow){nZeroBias_trig2off++;}
	if(trig13 && trig2){nZeroBias_trig2on++;}



	/*	
	if(trig1 && !trig2 && fourLayersHit && threeInRow)
	{
		std::cout << " " << std::endl;
		std::cout << "Event with trig1 passing online and trig2 failing online, with both found offline: " << event << std::endl;
		
		for(int i = 0; i < trig1_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig1
		{
			std::cout << "Combination " << i << " satisfying fourLayersHit logic." << std::endl;
			
			std::cout << "Time of pulses in trig1 combo: " << std::endl;
			for(int j = 0; j < overlappingHits[trig1_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
			{	
				std::cout << "Chan: " << (*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << std::endl;
				
			}
		}
		for(int i = 0; i < trig2_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig2
                {
                        std::cout << "Combination " << i << " satisfying threeInRow logic." << std::endl;
                	
			std::cout << "Time of pulses in trig2 combo: " << std::endl;
			for(int j = 0; j < overlappingHits[trig2_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
                        {
				std::cout << "Chan: " << (*chan)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << std::endl;
                        }
		}
	}
	*/

	/*
	if(trig1 && !trig2 && fourLayersHit && threeInRow)
        {
                std::cout << " " << std::endl;
		std::cout << "Event with trig1 passing online and trig2 failing online, with both found offline: " << event << std::endl;
        	for(int i = 0; i < trig1_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig1
                {
                        std::cout << "Combination " << i << " satisfying fourLayersHit logic." << std::endl;
                        //std::cout << "First pulse channel: " << (*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]].front()] << std::endl;
                        //std::cout << "First pulse start time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]].front()] << std::endl;
                        //std::cout << "Last pulse start time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]].back()] << std::endl;
                        //std::cout << "Last pulse end time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]].back()] + (*duration)[overlappingHits[trig1_coincidence_satisfying_comb[i]].back()] << std::endl;
                
			std::cout << "Time of pulses in trig1 combo: " << std::endl;
                	for(int j = 0; j < overlappingHits[trig1_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
                	{
				std::cout << "Chan: " << (*chan)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig1_coincidence_satisfying_comb[i]][j]] << std::endl;       
                	}
                }
                //std::cout << " " << std::endl;
                for(int i = 0; i < trig2_coincidence_satisfying_comb.size(); i++) //Looping over indices of combs that satisfy trig1
                {
                        std::cout << "Combination " << i << " satisfying threeInRow logic." << std::endl;
                        //std::cout << "First pulse channel: " << (*chan)[overlappingHits[trig2_coincidence_satisfying_comb[i]].front()] << std::endl;
                        //std::cout << "First pulse start time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]].front()] << std::endl;
                        //std::cout << "Last pulse start time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]].back()] << std::endl;
                        //std::cout << "Last pulse end time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]].back()] + (*duration)[overlappingHits[trig2_coincidence_satisfying_comb[i]].back()] << std::endl;
                	
			std::cout << "Time of pulses in trig2 combo: " << std::endl;
			for(int j = 0; j < overlappingHits[trig2_coincidence_satisfying_comb[i]].size(); j++) //Looping over pulses in coincident window
                        {
				std::cout << "Chan: " << (*chan)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << " Time: " << (*time)[overlappingHits[trig2_coincidence_satisfying_comb[i]][j]] << std::endl;
                        }
		}
	}
	*/


	//Check online/offline mismatch
    	/*
    	if(fourLayersHit != trig1)
    	{
    		std::cout << " " << std::endl;
    		std::cout << "Potential trigger1 mistake at event "<< event << std::endl;
    		std::cout << "trig1 = " << trig1 << " fourLayersHit = " << fourLayersHit << std::endl;
    		//std::cout << "layer0hit " << layer0hit << " layer1hit " << layer1hit << " layer2hit " << layer2hit << " layer3hit " << layer3hit << std::endl;

		
    		for(size_t j = 0; j < nHits; ++j)
		{
			std::cout << "chan: " << (*chan)[j] << std::endl;
        	        std::cout << "vMax_online = vMax_offline + dynamicPedestal = " << (*height)[j] << " + " << (*dynamicPedestal)[(*chan)[j]] << " = " << (*height)[j] + (*dynamicPedestal)[(*chan)[j]] << std::endl;
			std::cout << "time: " << (*time)[j] << std::endl;
		}
		
    	}*/
    	
		//=================================================




    	/*
    	if(threeInRow != trig2)
        {
                std::cout << "Potential trigger2 mistake at event "<< event << std::endl;
                std::cout << "trig2 = " << trig2 << " threeInRow = " << threeInRow << std::endl;
                std::cout << "threeLayersHit " << threeLayersHit << std::endl;
                //std::cout << "layer0hit " << layer0hit << " layer1hit " << layer1hit << " layer2hit " << layer2hit << " layer3hit " << layer3hit << std::endl;
		for(size_t j = 0; j < nHits; ++j)
                {
                        std::cout << "chan: " << (*chan)[j] << std::endl;
                        std::cout << "vMax_online = vMax_offline + dynamicPedestal = " << (*height)[j] << " + " << (*dynamicPedestal)[(*chan)[j]] << " = " << (*height)[j] + (*dynamicPedestal)[(*chan)[j]] << std::endl;
                        std::cout << "time: " << (*time)[j] << std::endl;
                }
        }
        */
    	//=================================================




    	
        /*
        if(trig3 && !twoSeparatedLayers)
        {
        	std::cout << "Potential trigger3 mistake at event "<< event << std::endl;
            std::cout << "trig3 = " << trig3 << " twoSeparatedLayers = " << twoSeparatedLayers << std::endl;
            std::cout << "layer0hit " << layer0hit << " layer1hit " << layer1hit << " layer2hit " << layer2hit << " layer3hit " << layer3hit << std::endl;
        }
        */
    	//=================================================




        /*
        if(trig4 && !twoAdjacentLayers)
        {
        	std::cout << "Potential trigger4 mistake at event "<< event << std::endl;
            std::cout << "trig4 = " << trig4 << " twoAdjacentLayers = " << twoAdjacentLayers << std::endl;
            std::cout << "layer0hit " << layer0hit << " layer1hit " << layer1hit << " layer2hit " << layer2hit << " layer3hit " << layer3hit << std::endl;
        }
        */
    	//=================================================




    	
        /*
        if(trig5 && !NLayersHit)
        {
        	std::cout << "Potential trigger5 mistake at event "<< event << std::endl;
            std::cout << "trig5 = " << trig5 << " NLayersHit = " << NLayersHit << std::endl;
            std::cout << "layer0hit " << layer0hit << " layer1hit " << layer1hit << " layer2hit " << layer2hit << " layer3hit " << layer3hit << std::endl;
        }
		*/
    	//=================================================




    	//==========Checking sixth trigger logic===========

    	//=================================================




        /*
        if(trig7 && !gtNHits)
        {
        	std::cout << "Potential trigger7 mistake at event "<< event << std::endl;
            std::cout << "trig7 = " << trig7 << " gtNHits = " << gtNHits << std::endl;
        }
		*/
    	//=================================================


    	//==========Checking eighth trigger logic==========

    	//=================================================


        
        /*
        if(trig9 && !topPanels)
        {
        	std::cout << "Potential trigger9 mistake at event "<< event << std::endl;
            std::cout << "trig9 = " << trig9 << " topPanels = " << topPanels << std::endl;
        }
		*/
    	//=================================================

	
        
        /*
        if(trig10 && !topPanels_plus_BottomBars)
        {
        	std::cout << "Potential trigger10 mistake at event "<< event << std::endl;
            std::cout << "trig10 = " << trig10 << " topPanels_plus_BottomBars = " << topPanels_plus_BottomBars << std::endl;
        }
        */
    	//=================================================


        
        /*
        if(front_back_panels != trig11)
        {
                std::cout << "Potential trigger11 mistake at event "<< event << std::endl;
                std::cout << "trig11 = " << trig11 << " front_back_panels = " << front_back_panels << std::endl;
        }
        */
    	//=================================================


	




    }  //End of event loop


//Calculate online/offline fractions
if(nTrig1_online > 0) frac1 = double(nTrig1_online)/double(nTrig1_offline);
if(nTrig2_online > 0) frac2 = double(nTrig2_online)/double(nTrig2_offline);
if(nTrig3_online > 0) frac3 = double(nTrig3_online)/double(nTrig3_offline);
if(nTrig4_online > 0) frac4 = double(nTrig4_online)/double(nTrig4_offline);
if(nTrig5_online > 0) frac5 = double(nTrig5_online)/double(nTrig5_offline);
if(nTrig6_online > 0) frac6 = double(nTrig6_online)/double(nTrig6_offline);
if(nTrig7_online > 0) frac7 = double(nTrig7_online)/double(nTrig7_offline);
if(nTrig8_online > 0) frac8 = double(nTrig8_online)/double(nTrig8_offline);
if(nTrig9_online > 0) frac9 = double(nTrig9_online)/double(nTrig9_offline);
if(nTrig10_online > 0) frac10 = double(nTrig10_online)/double(nTrig10_offline);
if(nTrig11_online > 0) frac11 = double(nTrig11_online)/double(nTrig11_offline);
if(nTrig12_online > 0) frac12 = double(nTrig12_online)/double(nTrig12_offline);
if(nTrig13_online > 0) frac13 = double(nTrig13_online)/double(nTrig13_offline);


std::cout << "Zerobias and fourLayersHit offline: " << nZeroBias_trig1off++ << std::endl;
std::cout << "Zerobias and fourLayersHit online: " << nZeroBias_trig1on++ << std::endl;
std::cout << "Zerobias and threeInRow offline: " << nZeroBias_trig2off++ << std::endl;
std::cout << "Zerobias and threeInRow online: " << nZeroBias_trig2on++ << std::endl;


std::cout << std::endl;
int width = 17;
std::cout << "Trigger Name-----nTrig_online-----nTrig_offline----prescale---------online/offline fraction" << std::endl;
std::cout << std::setw(width) << std::left << "4LayersHit" << std::setw(width) << nTrig1_online << std::setw(width) << nTrig1_offline << std::setw(width) << prescale1 << std::setw(width) << frac1 << std::endl;
std::cout << std::setw(width) << "threeInRow" << std::setw(width) << nTrig2_online << std::setw(width) << nTrig2_offline << std::setw(width) << prescale2 << std::setw(width) << frac2 << std::endl;
std::cout << std::setw(width) << "twoSeparatedLayer" << std::setw(width) << nTrig3_online << std::setw(width) << nTrig3_offline << std::setw(width) << prescale3 << std::setw(width) << frac3 << std::endl;
std::cout << std::setw(width) << "twoAdjacentLayer" << std::setw(width) << nTrig4_online << std::setw(width) << nTrig4_offline << std::setw(width) << prescale4 << std::setw(width) << frac4 << std::endl;
std::cout << std::setw(width) << "NLayersHit" << std::setw(width) << nTrig5_online << std::setw(width) << nTrig5_offline << std::setw(width) << prescale5 << std::setw(width) << frac5 << std::endl;
std::cout << std::setw(width) << "External" << std::setw(width) << nTrig6_online << std::setw(width) << nTrig6_offline << std::setw(width) << prescale6 << std::setw(width) << frac6 << std::endl;
std::cout << std::setw(width) << "gtNHits" << std::setw(width) << nTrig7_online << std::setw(width) << nTrig7_offline << std::setw(width) << prescale7 << std::setw(width) << frac7 << std::endl;
std::cout << std::setw(width) << "internal" << std::setw(width) << nTrig8_online << std::setw(width) << nTrig8_offline << std::setw(width) << prescale8 << std::setw(width) << frac8 << std::endl;
std::cout << std::setw(width) << "topPanels" << std::setw(width) << nTrig9_online << std::setw(width) << nTrig9_offline << std::setw(width) << prescale9 << std::setw(width) << frac9 << std::endl;
std::cout << std::setw(width) << "topPanel+BotBar" << std::setw(width) << nTrig10_online << std::setw(width) << nTrig10_offline << std::setw(width) << prescale10 << std::setw(width) << frac10 << std::endl;
std::cout << std::setw(width) << "front_back_panels" << std::setw(width) << nTrig11_online << std::setw(width) << nTrig11_offline << std::setw(width) << prescale11 << std::setw(width) << frac11 << std::endl;






TCanvas *c1 = new TCanvas();
h_trigger_rates_on->SetStats(0);
h_trigger_rates_off->SetStats(0);
h_trigger_rates_on->SetLineColor(kRed);
h_trigger_rates_off->SetLineColor(kGreen);
auto rp = new TRatioPlot(h_trigger_rates_on, h_trigger_rates_off);
rp->Draw();
//rp->GetUpperRefYaxis()->SetRangeUser(0., 1000.);
rp->GetUpperPad()->cd();
TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
leg->AddEntry(h_trigger_rates_on,"Online","l");
leg->AddEntry(h_trigger_rates_off,"Offline","l");
leg->Draw();
c1->SaveAs("trigger_hist.png");


TCanvas *c2 = new TCanvas();
h_ctrl_deltaT->Draw();
c2->SetLogy();
c2->SaveAs("h_ctrl_deltaT.png");


TCanvas *c3 = new TCanvas();
h_mismatched_deltaT->Draw("hist");
c3->SetLogy();
c3->SaveAs("h_mismatched_deltaT.png");

TCanvas *c4 = new TCanvas();
h_ctrl_vMax->Draw("hist");
c4->SetLogy();
c4->SaveAs("h_ctrl_vMax.png");

TCanvas *c5 = new TCanvas();
h_mismatched_vMax->Draw("hist");
c5->SetLogy();
c5->SaveAs("h_mismatched_vMax.png");

TCanvas *c6 = new TCanvas();
h_ctrl_area->Draw("hist");
c6->SetLogy();
c6->SaveAs("h_ctrl_area.png");

TCanvas *c7 = new TCanvas();
h_mismatched_area->Draw("hist");
c7->SetLogy();
c7->SaveAs("h_mismatched_area.png");

TCanvas *c8 = new TCanvas();
h_ctrl_chan->Draw("hist");
c8->SaveAs("h_ctrl_chan.png");

TCanvas *c9 = new TCanvas();
h_mismatched_chan->Draw("hist");
c9->SaveAs("h_mismatched_chan.png");

TCanvas *c10 = new TCanvas();
h_ctrl_startTime->Draw("hist");
c10->SetLogy();
c10->SaveAs("h_ctrl_startTime.png");

TCanvas *c11 = new TCanvas();
h_mismatched_startTime->Draw("hist");
c11->SetLogy();
c11->SaveAs("h_mismatched_startTime.png");


TCanvas *c12 = new TCanvas();
h_ctrl_height_v_area->Draw("colz");
c12->SaveAs("h_ctrl_height_v_area.png");

TCanvas *c13 = new TCanvas();
h_mismatched_height_v_area->Draw("colz");
c13->SaveAs("h_mismatched_height_v_area.png");



}










