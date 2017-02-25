#if defined(__ROOTCLING__) || defined(__ACLIC__)

#ifndef TreeAnalyzer_cc
#define TreeAnalyzer_cc

#include "Analyzers/ttH_analyzer/interface/TreeAnalyzer.h"

TreeAnalyzer::TreeAnalyzer(TTree* tree, Analysis_types analysis,
						   Selection_types selection, bool isdata)
{
	_tree = tree;
	_AnaType = analysis;
	_SelType = selection;
	_isdata = isdata;

	// SFHelper for updating scale factors
	_sfhelper = new SFHelper(_AnaType, _SelType, _isdata);
	
	// set branch address
	_ntuple.set_branch_address(_tree);

	_fourvectorsbuilt = false;
}

TreeAnalyzer::~TreeAnalyzer()
{
	delete _sfhelper;
}

void TreeAnalyzer::fill_Datacards_MC(std::map<TString,TH1D*>& hists, TTree* tree)
{
	assert(hists.size());

	int nEntries = tree->GetEntries();

	// loop over events in the tree
	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);
		
		if (not passTriggers()) continue;
		if (not passFilters()) continue;

		buildFourVectors();

		if (not passAdditionalSelection()) continue;
		
		// update bin index
		int ib = _ntuple.ibin;
		
		// update weights here
		updateWeights();
		
		// Fill the histograms
		// Assume keys of the histogram map are in the following format:
		// (_HiggsDecayMode)+(_gentau/_faketau)+(_systSuffix)
		// empty string "" for the central inclusive one
		for (auto & hist : hists) {
			TString key = hist.first;
			// filters
			if (key.Contains("_htt") and abs(_ntuple.HiggsDecayType) != 15)
				continue;
			if (key.Contains("_hzz") and abs(_ntuple.HiggsDecayType) != 23)
				continue;
			if (key.Contains("_hww") and abs(_ntuple.HiggsDecayType) != 24)
				continue;
			if (key.Contains("_gentau") and not _ntuple.isGenMatched) continue;
			if (key.Contains("_faketau") and _ntuple.isGenMatched) continue;

			// decide which weight to use
			float w_sf = 1;
			if (key.EndsWith("_LFUp"))
				w_sf = _btagSF_weight_LFUp / _bTagSF_weight;
			else if (key.EndsWith("_LFDown"))
				w_sf = _btagSF_weight_LFDown / _bTagSF_weight;
			else if (key.EndsWith("_HFUp"))
				w_sf = _btagSF_weight_HFUp / _bTagSF_weight;
			else if (key.EndsWith("_HFDown"))
				w_sf = _btagSF_weight_HFDown / _bTagSF_weight;
			else if (key.EndsWith("_HFStats1Up"))
				w_sf = _btagSF_weight_HFStats1Up / _bTagSF_weight;
			else if (key.EndsWith("_HFStats1Down"))
				w_sf = _btagSF_weight_HFStats1Down / _bTagSF_weight;
			else if (key.EndsWith("_HFStats2Up"))
				w_sf = _btagSF_weight_HFStats2Up / _bTagSF_weight;
			else if (key.EndsWith("_HFStats2Down"))
				w_sf = _btagSF_weight_HFStats2Down / _bTagSF_weight;
			else if (key.EndsWith("_LFStats1Up"))
				w_sf = _btagSF_weight_LFStats1Up / _bTagSF_weight;
			else if (key.EndsWith("_LFStats1Down"))
				w_sf = _btagSF_weight_LFStats1Down / _bTagSF_weight;
			else if (key.EndsWith("_LFStats2Up"))
				w_sf = _btagSF_weight_LFStats2Up / _bTagSF_weight;
			else if (key.EndsWith("_LFStats2Down"))
				w_sf = _btagSF_weight_LFStats2Down / _bTagSF_weight;
			else if (key.EndsWith("_cErr1Up"))
				w_sf = _btagSF_weight_cErr1Up / _bTagSF_weight;
			else if (key.EndsWith("_cErr1Down"))
				w_sf = _btagSF_weight_cErr1Down / _bTagSF_weight;
			else if (key.EndsWith("_cErr2Up"))
				w_sf = _btagSF_weight_cErr2Up / _bTagSF_weight;
			else if (key.EndsWith("_cErr2Down"))
				w_sf = _btagSF_weight_cErr2Down / _bTagSF_weight;
			else if (key.EndsWith("_x1Up"))
				w_sf = _MC_weight_scale_muF2 / _MC_weight;
			else if (key.EndsWith("_x1Down"))
				w_sf = _MC_weight_scale_muF0p5 / _MC_weight;
			else if (key.EndsWith("_y1Up"))
				w_sf = _MC_weight_scale_muR2 / _MC_weight;
			else if (key.EndsWith("_y2Down"))
				w_sf = _MC_weight_scale_muR0p5 / _MC_weight;

			hist.second -> Fill(ib, _event_weight * w_sf);
		} // end of keys loop
	} // end of event loop

	return;
}

void TreeAnalyzer::fill_Datacards_Data(TH1D* h, TTree* tree, vector<vector<unsigned long long>>& eventList)
{
	int nEntries = tree->GetEntries();

	// loop over events in the tree
	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);

		if (not passTriggers()) continue;
		if (not passFilters()) continue;

		vector<unsigned long long> eventid =
			{static_cast<unsigned long long>(_ntuple.run),
			 static_cast<unsigned long long>(_ntuple.ls), _ntuple.nEvent};
		
		bool alreadyIncluded =
			find(eventList.begin(), eventList.end(), eventid) != eventList.end();

		if (alreadyIncluded) continue;

		eventList.push_back(eventid);

		buildFourVectors();
		
		if (not passAdditionalSelection()) continue;
		
		// update bin index
		int ib = _ntuple.ibin;
		
		// update weights here
		updateWeights();

		// fill histogram
		h->Fill(ib, _event_weight);
	}
}

void TreeAnalyzer::dump_Events(TTree* tree, TString channel)
{
	std::ofstream eventDump;

	TString fname = "EventDump_"+channel+".py";
	
	eventDump.open(fname);
	cout << "File " << fname <<" is created." << endl;
	
	eventDump << "\"\"\"Event list\n";
	eventDump << "Per event :\n";
	eventDump << " - leptons\n";
	eventDump << " - tau\n";
	eventDump << " - MET\n";
	eventDump << " - MET Cov. matrix\n";
	eventDump << " - 2 b-jets\n";
	eventDump << " - Untagged jets\n";
	eventDump << "\"\"\"\n\n";
	eventDump << "Events = [\n";

	int nEntries = tree->GetEntries();

	// loop over events in the tree
	for (int i = 0; i < nEntries; ++i) {
		tree->GetEntry(i);

		if (not passTriggers()) continue;
		if (not passFilters()) continue;

		buildFourVectors();

		if (not passAdditionalSelection()) continue;

		cout << "Event " << i << "passed the selection. ";
		cout << "Dumping its contents..." << endl;
		
		eventDump << "{\n";
		// leptons
		eventDump << "\'_recolep_sel_px\': ["
				  << _lep0.Px() << "," << _lep1.Px() << "],\n";
		eventDump << "\'_recolep_sel_py\': ["
				  << _lep0.Py() << "," << _lep1.Py() << "],\n";
		eventDump << "\'_recolep_sel_pz\': ["
				  << _lep0.Pz() << "," << _lep1.Pz() << "],\n";
		eventDump << "\'_recolep_sel_e\': ["
				  << _lep0.E() << "," << _lep1.E() << "],\n";
		eventDump << "\'_recolep_sel_pdg\': ["
				  << _leps_id[0] << "," << _leps_id[1] << "],\n";
		// tau
		eventDump << "\'_recotauh_sel_px\': [" << _tau.Px() << "],\n";
		eventDump << "\'_recotauh_sel_py\': [" << _tau.Py() << "],\n";
		eventDump << "\'_recotauh_sel_pz\': [" << _tau.Pz() << "],\n";
		eventDump << "\'_recotauh_sel_e\': [" << _tau.E() << "],\n";
		eventDump << "\'_recotauh_sel_decayMode\': ["
				  << _ntuple.tau0_decayMode << "],\n";
		// MET
		eventDump << "\'_PFMET\': " << _ntuple.PFMET
				  << ", \'_PFMET_phi\': " << _ntuple.PFMETphi << ",\n";
		eventDump << "\'_PFMETCov00\': " << _ntuple.METCov00 << ", "
				  << "\'_PFMETCov11\': " << _ntuple.METCov11 << ",\n";
		eventDump << "\'_PFMETCov01\': " << _ntuple.METCov01 << ", "
				  << "\'_PFMETCov10\': " << _ntuple.METCov10 << ",\n";
		// btagged jets
		eventDump << "\'_recoPFJet_btag_px\': ["
				  << _bjet0.Px() << "," << _bjet1.Px() << "],\n";
		eventDump << "\'_recoPFJet_btag_py\': ["
				  << _bjet0.Py() << "," << _bjet1.Py() << "],\n";
		eventDump << "\'_recoPFJet_btag_pz\': ["
				  << _bjet0.Pz() << "," << _bjet1.Pz() << "],\n";
		eventDump << "\'_recoPFJet_btag_e\': ["
				  << _bjet0.E() << "," << _bjet1.E() << "],\n";
		// untagged jets
		eventDump << "\'_n_recoPFJet_untag\': " << _untag_jets.size() << ",\n";
		eventDump << "\'_recoPFJet_untag_px\': [";
		size_t ipx = 0;
		for (const auto & jet : _untag_jets) {
			eventDump << jet.Px();
			++ipx;
			if (ipx < _untag_jets.size())
				eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_py\': [";
		size_t ipy = 0;
		for (const auto & jet : _untag_jets) {
			eventDump << jet.Py();
			++ipy;
			if (ipy < _untag_jets.size())
				eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_pz\': [";
		size_t ipz = 0;
		for (const auto & jet : _untag_jets) {
			eventDump << jet.Pz();
			++ipz;
			if (ipz < _untag_jets.size())
				eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_e\': [";
		size_t ie = 0;
		for (const auto & jet : _untag_jets) {
			eventDump << jet.E();
			++ie;
			if (ie < _untag_jets.size())
				eventDump << ",";
		}
		eventDump << "]\n";
		
		eventDump << "},\n";
	}

	eventDump << "]";
	eventDump.close();
}

vector<TH1D*> TreeAnalyzer::makeHistograms(TTree* tree, bool isdata, vector<vector<unsigned long long>>& eventList)
{
	vector<TH1D*> out_hists;

	return out_hists;
}

void TreeAnalyzer::buildFourVectors()
{
	// leptons
	if (_ntuple.lepCategory == 0) {// mumu
		_lep0.SetPtEtaPhiE(_ntuple.mu0_pt, _ntuple.mu0_eta, 
						   _ntuple.mu0_phi, _ntuple.mu0_E);
		_lep1.SetPtEtaPhiE(_ntuple.mu1_pt, _ntuple.mu1_eta, 
						   _ntuple.mu1_phi, _ntuple.mu1_E);
		_leps_istight[0] = _ntuple.mu0_ismvasel;
		_leps_istight[1] = _ntuple.mu1_ismvasel;
		_leps_charge[0] = _ntuple.mu0_charge;
		_leps_charge[1] = _ntuple.mu1_charge;
		_leps_id[0] = -13*_leps_charge[0]; _leps_id[1] = -13*_leps_charge[1];
	}
	else if (_ntuple.lepCategory == 1) {// ee
		_lep0.SetPtEtaPhiE(_ntuple.ele0_pt, _ntuple.ele0_eta, 
						   _ntuple.ele0_phi, _ntuple.ele0_E);
		_lep1.SetPtEtaPhiE(_ntuple.ele1_pt, _ntuple.ele1_eta, 
						   _ntuple.ele1_phi, _ntuple.ele1_E);
		_leps_istight[0] = _ntuple.ele0_ismvasel;
		_leps_istight[1] = _ntuple.ele1_ismvasel;
		_leps_charge[0] = _ntuple.ele0_charge;
		_leps_charge[1] = _ntuple.ele1_charge;
		_leps_id[0] = -11*_leps_charge[0]; _leps_id[1] = -11*_leps_charge[1];
	}
	else if (_ntuple.lepCategory == 2) {// emu
		if (_ntuple.ele0_pt > _ntuple.mu0_pt) {
			_lep0.SetPtEtaPhiE(_ntuple.ele0_pt, _ntuple.ele0_eta, 
							   _ntuple.ele0_phi, _ntuple.ele0_E);
			_lep1.SetPtEtaPhiE(_ntuple.mu0_pt, _ntuple.mu0_eta, 
							   _ntuple.mu0_phi, _ntuple.mu0_E);
			_leps_istight[0] = _ntuple.ele0_ismvasel;
			_leps_istight[1] = _ntuple.mu0_ismvasel;
			_leps_charge[0] = _ntuple.ele0_charge;
			_leps_charge[1] = _ntuple.mu0_charge;
			_leps_id[0] = -11*_leps_charge[0]; _leps_id[1] = -13*_leps_charge[1];
		}
		else {
			_lep0.SetPtEtaPhiE(_ntuple.mu0_pt, _ntuple.mu0_eta, 
							   _ntuple.mu0_phi, _ntuple.mu0_E);
			_lep1.SetPtEtaPhiE(_ntuple.ele0_pt, _ntuple.ele0_eta, 
							   _ntuple.ele0_phi, _ntuple.ele0_E);
			_leps_istight[0] = _ntuple.mu0_ismvasel;
			_leps_istight[1] = _ntuple.ele0_ismvasel;
			_leps_charge[0] = _ntuple.mu0_charge;
			_leps_charge[1] = _ntuple.ele0_charge;
			_leps_id[0] = -13*_leps_charge[0]; _leps_id[1] = -11*_leps_charge[1];
		}
	}

	_leps_conept[0] = _ntuple.lep0_conept;
	_leps_conept[1] = _ntuple.lep1_conept;

	// tau
	_tau.SetPtEtaPhiE(_ntuple.tau0_pt, _ntuple.tau0_eta, _ntuple.tau0_phi, _ntuple.tau0_E);

	// jets
	int njets = _ntuple.jets_pt.size();
	
	int icsv0 = -1; int icsv1 = -1;
	float csv0 = -10.; float csv1 = -10.;
	
	for (int i = 0; i < njets; ++i) {
		if (_ntuple.jets_csv[i] > csv0) {
			csv1 = csv0;
			icsv1 = icsv0;
			csv0 = _ntuple.jets_csv[i];
			icsv0 = i;
		}
		else if (_ntuple.jets_csv[i] > csv1) {
			csv1 = _ntuple.jets_csv[i];
			icsv1 = i;
		}
	}

	assert(icsv0 > -1 and icsv1 > -1);
	_bjet0.SetPtEtaPhiE(_ntuple.jets_pt[icsv0],_ntuple.jets_eta[icsv0],
						_ntuple.jets_phi[icsv0],_ntuple.jets_E[icsv0]);
	_bjet1.SetPtEtaPhiE(_ntuple.jets_pt[icsv1],_ntuple.jets_eta[icsv1],
						_ntuple.jets_phi[icsv1],_ntuple.jets_E[icsv1]);

	for (int i = 0; i < njets; ++i) {
		if (i == icsv0 or i == icsv1) continue;
		TLorentzVector jet;
		jet.SetPtEtaPhiE(_ntuple.jets_pt[i],_ntuple.jets_eta[i],
						 _ntuple.jets_phi[i],_ntuple.jets_E[i]);
		_untag_jets.push_back(jet);
	}
	
	_fourvectorsbuilt = true;
}

void TreeAnalyzer::updateWeights()
{
	if (not _isdata) {
		// pileup
		_PU_weight = _sfhelper->Get_PUWeight(_ntuple.npuTrue);

		// trigger SF
		_triggerSF_weight = _sfhelper->Get_HLTSF(_ntuple.lepCategory);
		
		// tau scale factor
		_tauSF_weight = _sfhelper->Get_TauIDSF(_ntuple.tau0_pt,_ntuple.tau0_eta,
											   _ntuple.isGenMatched);

		// lepton ID scale factors
		_leptonSF_weight =
			_sfhelper->Get_LeptonIDSF(_lep0.Pt(), _lep0.Eta(),
									  abs(_leps_id[0])==11, abs(_leps_id[0])==13,
									  _leps_istight[0]);
		_leptonSF_weight *=
			_sfhelper->Get_LeptonIDSF(_lep1.Pt(), _lep1.Eta(),
									  abs(_leps_id[1])==11, abs(_leps_id[1])==13,
									  _leps_istight[1]);

		// btag
		// TODO
		// for now
		_bTagSF_weight = _ntuple.bTagSF_weight;
		_btagSF_weight_LFUp = _ntuple.btagSF_weight_LFUp;
		_btagSF_weight_LFDown = _ntuple.btagSF_weight_LFDown;
		_btagSF_weight_HFUp = _ntuple.btagSF_weight_HFUp;
		_btagSF_weight_HFDown = _ntuple.btagSF_weight_HFDown;
		_btagSF_weight_HFStats1Up = _ntuple.btagSF_weight_HFStats1Up;
		_btagSF_weight_HFStats1Down = _ntuple.btagSF_weight_HFStats1Down;
		_btagSF_weight_HFStats2Up = _ntuple.btagSF_weight_HFStats2Up;
		_btagSF_weight_HFStats2Down = _ntuple.btagSF_weight_HFStats2Down;
		_btagSF_weight_LFStats1Up = _ntuple.btagSF_weight_LFStats1Up;
		_btagSF_weight_LFStats1Down = _ntuple.btagSF_weight_LFStats1Down;
		_btagSF_weight_LFStats2Up = _ntuple.btagSF_weight_LFStats2Up;
		_btagSF_weight_LFStats2Down = _ntuple.btagSF_weight_LFStats2Down;
		_btagSF_weight_cErr1Up = _ntuple.btagSF_weight_cErr1Up;
		_btagSF_weight_cErr1Down = _ntuple.btagSF_weight_cErr1Down;
		_btagSF_weight_cErr2Up = _ntuple.btagSF_weight_cErr2Up;
		_btagSF_weight_cErr2Down = _ntuple.btagSF_weight_cErr2Down;

		// MC weight
		_MC_weight = _ntuple.MC_weight;
		_MC_weight_scale_muF0p5 = _ntuple.MC_weight_scale_muF0p5;
		_MC_weight_scale_muF2 = _ntuple.MC_weight_scale_muF2;
		_MC_weight_scale_muR0p5 = _ntuple.MC_weight_scale_muR0p5;
		_MC_weight_scale_muR2 = _ntuple.MC_weight_scale_muR2;

		// total event weight
		_event_weight = _PU_weight * _triggerSF_weight * _tauSF_weight *
			_leptonSF_weight * _bTagSF_weight * _MC_weight;
	}
	else {
		_event_weight = _ntuple.event_weight;

		if (_SelType==Control_1lfakeable) {  // fakes
			_event_weight =
				_sfhelper->Get_FR_weight(_leps_conept[0], _lep0.Eta(),
										 abs(_leps_id[0])==11,
										 abs(_leps_id[0])==13,
										 static_cast<bool>(_leps_istight[0]),
										 _leps_conept[1], _lep1.Eta(),
										 abs(_leps_id[1])==11,
										 abs(_leps_id[1])==13,
										 static_cast<bool>(_leps_istight[1]));
		}
		else if (_SelType==Control_2los1tau) {  // flips
			float p1 = _sfhelper->Get_EleChargeMisIDProb(
			    _ntuple.ele0_pt,_ntuple.ele0_eta,_ntuple.ele0_charge,_ntuple.tau0_charge);
			float p2 = _sfhelper->Get_EleChargeMisIDProb(
				_ntuple.ele1_pt,_ntuple.ele1_eta,_ntuple.ele1_charge,_ntuple.tau0_charge);
			if (_ntuple.lepCategory == 0) // mumu
				_event_weight = 0.;
			else if (_ntuple.lepCategory == 1) // ee
				_event_weight = p1 + p2;
			else if (_ntuple.lepCategory == 2) // emu
				_event_weight = p1;
		}
	}
}

bool TreeAnalyzer::passTriggers()
{
	// check trigger bits if necessary
	// _ntuple.triggerBits

	return _ntuple.matchHLTPath;
}

bool TreeAnalyzer::passFilters()
{
	// MET filter suggestions: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
	// list of filter names contained in the bits can be found in
	// python/ttHtaus_cfi.py

	// event tagged as bad muons are corrected rather than removing
	
	if (_isdata) {
		return _ntuple.filterBits & 63;
	}
	else {
		return _ntuple.filterBits & 47;
	}

	return false;
}

bool TreeAnalyzer::passAdditionalSelection(bool controlRegion)
{
	assert(_fourvectorsbuilt);

	bool passSel = false;
	
	if (_SelType != Control_2los1tau) { 
		// lepton charge and tau charge are oppostie sign
		assert(_leps_charge[0]==_leps_charge[1]);
		passSel = (_leps_charge[0] + _ntuple.tau0_charge == 0);
	}

	if (controlRegion)
		passSel = !passSel;

	return passSel;
}

#endif
#endif
