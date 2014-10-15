/* Simple xAna analysis example. */

#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <stdio.h>
#include <fstream>
#include <TProfile.h>


#include "../HiddenFile/untuplizer.h"
#include "../HiddenFile/ElectronSelections.h"
#include "../HiddenFile/MuonSelections.h"
#include "../HiddenFile/PhotonSelections.h"
#include <iostream> 
#include <vector>
#include "rochcor2012jan22.h"
#include "../HiddenFile/puweicalc.h"
#include "/afs/cern.ch/user/c/ctu/work/ele_veto/Zgammagamma/SF_Modify/SF_Modify_Photon.C"
#include "/afs/cern.ch/user/c/ctu/work/ele_veto/Zgammagamma/SF_Modify/SF_Modify_SingleMu.C"
#include "/afs/cern.ch/user/c/ctu/work/ele_veto/Zgammagamma/SF_Modify/SF_Modify_diMu.C"
#include "/afs/cern.ch/user/c/ctu/CMSSW_5_3_7/src/MyPhosphorDir/PhosphorCorrectorFunctor.hh"

using namespace std;

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}

void xAna_simplify(Int_t MC_N, Int_t MC, Int_t WP=1) {
  // if MC == 0, it's S10
  // if MC == 1, it's RD1
  rochcor2012 *rmcor = new rochcor2012(229);
  Char_t fname_File[30];
  sprintf(fname_File,"Result/Zrr_%d_short.root",MC_N);
  cout<<fname_File<<endl;
  
  TFile *fout_ = new TFile(fname_File,"recreate");
  TTree *outtree_ = new TTree("t", "mini tree");
  
  Float_t MuSel_1_Pt;
  Float_t MuSel_1_Eta;
  Float_t MuSel_1_Phi;

  Float_t MuSel_2_Pt;
  Float_t MuSel_2_Eta;
  Float_t MuSel_2_Phi;

  Float_t Phot_1_Pt;
  Float_t Phot_1_Eta;
  Float_t Phot_1_Phi;

  Float_t Phot_2_Pt;
  Float_t Phot_2_Eta;
  Float_t Phot_2_Phi;
  
  Float_t Zmass_;
  Float_t Zgmass_;
  Float_t Zggmass_;
  Float_t tpmass_;

  Float_t SFMu1;
  Float_t SFMu2;
  Float_t SFdiMu;
  Float_t SFPhot1;
  Float_t SFPhot2;
  Float_t SFResult;
  Float_t PUwei;

  Int_t   EventNumber;
  Int_t   run_;
  Long64_t   event_;

  Float_t phoSCRChIso1_;
  Float_t phoSCRChIso2_;
  Float_t phoRandConeChIso_;
  Float_t rho2012_;
  Float_t SCEta_;
 
  Float_t pho1eleveto_;
  Float_t pho2eleveto_;

  Float_t PIDcheck1_;
  Float_t PIDcheck2_;

  Float_t mcParentage1_;
  Float_t mcParentage2_;

  Float_t deltaR1_;
  Float_t deltaR2_;

  outtree_->Branch("MuSel_1_Pt",  &MuSel_1_Pt,  "MuSel_1_Pt/F");
  outtree_->Branch("MuSel_1_Eta", &MuSel_1_Eta, "MuSel_1_Eta/F");
  outtree_->Branch("MuSel_1_Phi", &MuSel_1_Phi, "MuSel_1_Phi/F");
  outtree_->Branch("MuSel_2_Pt",  &MuSel_2_Pt,  "MuSel_2_Pt/F");
  outtree_->Branch("MuSel_2_Eta", &MuSel_2_Eta, "MuSel_2_Eta/F");
  outtree_->Branch("MuSel_2_Phi", &MuSel_2_Phi, "MuSel_2_Phi/F");

  outtree_->Branch("Phot_1_Pt",  &Phot_1_Pt,  "Phot_1_Pt/F");
  outtree_->Branch("Phot_1_Eta", &Phot_1_Eta, "Phot_1_Eta/F");
  outtree_->Branch("Phot_1_Phi", &Phot_1_Phi, "Phot_1_Phi/F");
  outtree_->Branch("Phot_2_Pt",  &Phot_2_Pt,  "Phot_2_Pt/F");
  outtree_->Branch("Phot_2_Eta", &Phot_2_Eta, "Phot_2_Eta/F");
  outtree_->Branch("Phot_2_Phi", &Phot_2_Phi, "Phot_2_Phi/F");

  outtree_->Branch("Zmass_", &Zmass_, "Zmass_/F");
  outtree_->Branch("Zgmass_", &Zgmass_, "Zgmass_/F");
  outtree_->Branch("Zggmass_", &Zggmass_, "Zggmass_/F");
  outtree_->Branch("tpmass_", &tpmass_, "tpmass_/F");

  outtree_->Branch("SFResult", &SFResult, "SFResult/F");
  outtree_->Branch("SFMu1", &SFMu1, "SFMu1/F");
  outtree_->Branch("SFMu2", &SFMu2, "SFMu2/F");
  outtree_->Branch("SFdiMu", &SFdiMu, "SFdiMu/F");
  outtree_->Branch("SFPhot1", &SFPhot1, "SFPhot1/F");
  outtree_->Branch("SFPhot2", &SFPhot2, "SFPhot2/F");

  outtree_->Branch("PUwei", &PUwei, "PUwei/F");
  outtree_->Branch("EventNumber", &EventNumber, "EventNumber/I");
  outtree_->Branch("run_", &run_, "run_/I");
  outtree_->Branch("event_", &event_, "event_/L");
  outtree_->Branch("phoSCRChIso1_", &phoSCRChIso1_, "phoSCRChIso1_/F");
  outtree_->Branch("phoSCRChIso2_", &phoSCRChIso2_, "phoSCRChIso2_/F");
  outtree_->Branch("rho2012_", &rho2012_, "rho2012_/F");
  outtree_->Branch("SCEta_", &SCEta_, "SCEta_/F");
  outtree_->Branch("phoRandConeChIso_", &phoRandConeChIso_, "phoRandConeChIso_/F");
  outtree_->Branch("pho1eleveto_", &pho1eleveto_, "pho1eleveto_/F");
  outtree_->Branch("pho2eleveto_", &pho2eleveto_, "pho2eleveto_/F");
  outtree_->Branch("PIDcheck1_", &PIDcheck1_, "PIDcheck1_/F");
  outtree_->Branch("PIDcheck2_", &PIDcheck2_, "PIDcheck2_/F");
  outtree_->Branch("mcParentage1_", &mcParentage1_, "mcParentage1_/F");
  outtree_->Branch("mcParentage2_", &mcParentage2_, "mcParentage2_/F");
  outtree_->Branch("deltaR1_", &deltaR1_, "deltaR1_/F");
  outtree_->Branch("deltaR2_", &deltaR2_, "deltaR2_/F");

  vector <string> pathes;

  if(MC_N==1) pathes.push_back("/data5/ggNtuples/V05-03-12-03/job_summer12_DYJetsToLL_s10.root");
  if(MC_N==2) pathes.push_back("/data5/ggNtuples/V05-03-12-06/job_summer12_Zg_s10.root");
  if(MC_N==7) pathes.push_back("/data5/ggNtuples/V05-03-12-05/job_summer12_Zgg.root");
  
  TreeReader data(pathes);
  cout<<data.GetEntriesFast()<<endl;

  vector<int> rvector;
  vector<int> lvector;
  vector<Long64_t> evector;
  vector<float> puvector;
  vector<float> pho1EtOvector;
  vector<float> pho2EtOvector;
  vector<float> pho1EtAvector;
  vector<float> pho2EtAvector;
  
  // pileup reweighting for MC
  PUWeightCalculator puCalcS10;
  PUWeightCalculator puCalcRD1;
  if (data.HasMC()) {
    puCalcS10.Init("../mcweights/mcwei_S10.root");
    puCalcRD1.Init("../mcweights/mcwei_PU_RD1_muo.root");
  }

  vector<int> ev_s;
  vector<double> check_Eta;
  vector<double> check_Et;
  vector<double> check_deltaR1;
  vector<double> check_deltaR2;

  Int_t DYTest = 0;
  
  /*
  Double_t L = 19.7;
  Double_t Ne[8] = {1,30458871,6588161,1497445,954911,1936727,800047,199982};
  Double_t XS[8] = {1,3503.71,156.2,0.1767,0.32,1.28,0.0,0.125};
  */

  // event loop
  //for (Long64_t ev = 0; ev < 3000000; ev++) {
  Float_t normalY = 1;
  Double_t L = 19.7;

  if(MC_N==1) normalY = 3503.7*L*1000/data.GetEntriesFast();
  if(MC_N==2) normalY = 159.12*L*1000/data.GetEntriesFast();
  if(MC_N==7) normalY = 0.125*L*1000/data.GetEntriesFast();

  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    // print progress
    if (ev % 1000000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);
    // IMPORTANT: branches with counters must be loaded BEFORE dependent branches
    // accessing HLT information
    Int_t nHLT = data.GetInt("nHLT");
    Int_t* HLT = data.GetPtrInt("HLT");
    Int_t* HLTIndex = data.GetPtrInt("HLTIndex");
    if(HLT[HLTIndex[14]]!=1) continue;
    //for (Int_t i = 0; i < nHLT; ++i) printf(" %i", HLTIndex[i]);
    //printf("\n");
    
    // muon selection
    Int_t nMu        = data.GetInt("nMu");
    Float_t* muPt    = data.GetPtrFloat("muPt");
    Float_t* muEta   = data.GetPtrFloat("muEta");
    Float_t* muPhi   = data.GetPtrFloat("muPhi");
    Int_t* muCharge  = data.GetPtrInt("muCharge");

    Int_t    run     = data.GetInt("run");
    Long64_t event    = data.GetLong64("event");
    Int_t    lumis      = data.GetInt("lumis");
    if(lumis!=7246) continue;
    if(event!=2173243) continue;
    cout<<"bitch"<<endl;

    Int_t    nMC         = 0;
    Int_t*   mcPID       = NULL;
    Int_t*   mcMomPID    = NULL;
    Int_t*   mcGMomPID   = NULL;
    Int_t*   mcParentage = NULL;
    Float_t* mcPt        = NULL;
    Float_t* mcEta       = NULL;
    Float_t* mcPhi       = NULL;
    Float_t* mcE         = NULL;

    if (data.HasMC()) {
      nMC         = data.GetInt("nMC");
      mcPID       = data.GetPtrInt("mcPID");
      mcMomPID    = data.GetPtrInt("mcMomPID");
      mcGMomPID   = data.GetPtrInt("mcGMomPID");
      mcParentage = data.GetPtrInt("mcParentage");
      mcPt        = data.GetPtrFloat("mcPt");
      mcEta       = data.GetPtrFloat("mcEta");
      mcPhi       = data.GetPtrFloat("mcPhi");
      mcE         = data.GetPtrFloat("mcE");
    }
    vector<int> acc_mu;
    select_muons(data, acc_mu);
    
    TLorentzVector mu1, mu2, Mmumu, Z, mu1_wa, mu2_wa;
    mu1_wa.SetPtEtaPhiM(0, 0, 0, 0.1057);
    mu2_wa.SetPtEtaPhiM(0, 0, 0, 0.1057);
    
    Double_t M_wa = 0;
    Int_t check_mu=0;
	
    Int_t runopt = 0;
    /*
    if(run >= 203894){
      runopt = 1;
    }
    */

    TLorentzVector Z_test;
    cout<<"Muons : "<<acc_mu.size()<<endl;
    
    for (size_t i = 0; i < acc_mu.size(); i++) {  
      for (size_t j = i + 1; j < acc_mu.size(); j++) {
	mu1.SetPtEtaPhiM(muPt[acc_mu[i]], muEta[acc_mu[i]], muPhi[acc_mu[i]], 0.1057);
	mu2.SetPtEtaPhiM(muPt[acc_mu[j]], muEta[acc_mu[j]], muPhi[acc_mu[j]], 0.1057);
	
	float qter = 1.0;
	if(data.HasMC()) {
	  rmcor->momcor_mc(mu1, muCharge[acc_mu[i]], runopt, qter); 
	  rmcor->momcor_mc(mu2, muCharge[acc_mu[j]], runopt, qter);
	}
	else{
	  rmcor->momcor_data(mu1, muCharge[acc_mu[i]], runopt, qter);
          rmcor->momcor_data(mu2, muCharge[acc_mu[j]], runopt, qter);
	}

	Mmumu = mu1 + mu2;
	cout<<"Mmumu : "<<Mmumu.M()<<endl;

	if (muCharge[acc_mu[i]]==muCharge[acc_mu[j]]) continue;
	if (mu1.Pt()<20||mu2.Pt()<10) continue;
	if (fabs(mu1.Eta())>2.4||fabs(mu2.Eta())>2.4) continue;
	if (Mmumu.M()<50) continue;
	
	if(check_mu==0||fabs(91.1876-M_wa)>fabs(91.1876-Mmumu.M())){
	  check_mu++;
	  mu1_wa=mu1;
	  mu2_wa=mu2;
	  M_wa = (mu1_wa+mu2_wa).M();
	}
      }
    }
    if (mu1_wa.Pt()==0||mu2_wa.Pt()==0) continue;
    Z      = mu1_wa + mu2_wa;
    if (Z.M()<50) continue;

    // photon selection
    Int_t nPho        = data.GetInt("nPho");
    Float_t* PhoEt    = data.GetPtrFloat("phoEt");
    Float_t* PhoE     = data.GetPtrFloat("phoE");
    Float_t* PhoEta   = data.GetPtrFloat("phoEta");
    Float_t* PhoPhi   = data.GetPtrFloat("phoPhi");
    Int_t* phoEleVeto = data.GetPtrInt("phoEleVeto");
    Float_t  rho2012     = data.GetFloat("rho2012");
    Float_t* phoR9     = data.GetPtrFloat("phoR9");
      
    TLorentzVector Phot, Phot_DY, Phot_wa1, Phot_wa2, Phot_two, Phot_MC;
    Phot_wa1.SetPtEtaPhiM(0,0,0,0);
    Phot_wa2.SetPtEtaPhiM(0,0,0,0);
    vector<int> acc_photon;
    select_photon(WP, data, acc_photon);
    
    Float_t* phoSCEta   = data.GetPtrFloat("phoSCEta");    
   	
    Int_t count_pho = 0;
    Int_t Check_SecPhoton = 0;

    Float_t SCEta_wa1 = 0;
    Float_t SCEta_wa2 = 0;
    Float_t phoSCRChIsoR1 = 0;
    Float_t phoSCRChIsoR2 = 0;
    Float_t phoRandConeChIsoR = 0;
    Float_t rho2012R = 0;
    Float_t pho1EtOR = 0;
    Float_t pho2EtOR = 0;
    Float_t pho1elevetoR = 0;
    Float_t pho2elevetoR = 0;

    Float_t deltaR1_R = 0;
    Float_t deltaR2_R = 0;
    
    Float_t puwei = 1.;
    if (data.HasMC()) {
      Float_t* puTrue = data.GetPtrFloat("puTrue");
      if (MC == 0) puwei = (Float_t) puCalcS10.GetWeight(run, puTrue[1]);
      if (MC == 1) puwei = (Float_t) puCalcRD1.GetWeight(run, puTrue[1]);
    }

    for (size_t g=0; g < acc_photon.size(); g++) {
      if (count_pho==2) continue;
      Phot.SetPtEtaPhiM(PhoEt[acc_photon[g]], PhoEta[acc_photon[g]], PhoPhi[acc_photon[g]],0);
      
      Int_t check_Zg = 0;
      if(MC_N==1) {
	int phoInd = -1;
	for(int imc = 0; imc < nMC; ++imc){
	  if(mcPID[imc] != 22) continue;
	  bool match_gen = deltaR(mcEta[imc], mcPhi[imc], PhoEta[acc_photon[g]], PhoPhi[acc_photon[g]]) < 0.2 &&(fabs(PhoEt[acc_photon[g]] - mcPt[imc]) / mcPt[imc] < 1.0);
	  if(match_gen && phoInd < 0) phoInd = imc;
	}
	if(phoInd >= 0){
	  if((mcParentage[phoInd]& 4) == 0) check_Zg = 1;
	  else
            check_Zg = 2;
	} else {
          check_Zg = 3;
	}
      }
      if(check_Zg==1) continue;
            
      if (fabs(phoSCEta[acc_photon[g]])>2.5) continue;
      if (fabs(phoSCEta[acc_photon[g]])<1.566&&fabs(phoSCEta[acc_photon[g]])>1.4442) continue;
      if (Phot.Pt()<15) continue;
      Float_t deltaR1 = deltaR(mu1_wa.Eta(),mu1_wa.Phi(),Phot.Eta(),Phot.Phi());
      Float_t deltaR2 = deltaR(mu2_wa.Eta(),mu2_wa.Phi(),Phot.Eta(),Phot.Phi());
      Float_t deltamin = deltaR1;
      if (deltaR2<deltaR1) deltamin = deltaR2;
      if (deltamin < 0.5) continue;
      if (phoEleVeto[acc_photon[g]]!=0) continue;
      
      if(count_pho==1){
	Phot_wa2.SetPtEtaPhiM(PhoEt[acc_photon[g]], PhoEta[acc_photon[g]], PhoPhi[acc_photon[g]],0);
	SCEta_wa2 = phoSCEta[acc_photon[g]];
	pho2EtOR = PhoEt[acc_photon[g]];
	pho2elevetoR = phoEleVeto[acc_photon[g]];
	rho2012R = rho2012;

	rvector.push_back(run);
        lvector.push_back(lumis);
        evector.push_back(event);
	puvector.push_back(puwei);	

	pho1EtOvector.push_back(pho1EtOR);
	pho2EtOvector.push_back(pho2EtOR);
	pho1EtAvector.push_back(Phot_wa1.Et());
	pho2EtAvector.push_back(Phot_wa2.Et());
	
	deltaR2_R = deltamin;
	count_pho++;
      }

      if(count_pho==0){
	Phot_wa1.SetPtEtaPhiM(PhoEt[acc_photon[g]], PhoEta[acc_photon[g]], PhoPhi[acc_photon[g]],0);
        SCEta_wa1 = phoSCEta[acc_photon[g]];
	pho1EtOR = PhoEt[acc_photon[g]];
	pho1elevetoR = phoEleVeto[acc_photon[g]];
	deltaR1_R = deltamin;
	count_pho++;
      }
    }
    if (Phot_wa1.Pt()==0||Phot_wa2.Pt()==0) continue;
    //cout<<"hi"<<endl;
    
    PUwei = puwei;
    MuSel_1_Pt  = mu1_wa.Pt();
    MuSel_1_Eta = mu1_wa.Eta();
    MuSel_1_Phi = mu1_wa.Phi();
    MuSel_2_Pt  = mu2_wa.Pt();
    MuSel_2_Eta = mu2_wa.Eta();
    MuSel_2_Phi = mu2_wa.Phi();
    Phot_1_Pt   = Phot_wa1.Pt();
    Phot_1_Eta  = Phot_wa1.Eta();
    Phot_1_Phi  = Phot_wa1.Phi();
    Phot_2_Pt   = Phot_wa2.Pt();
    Phot_2_Eta  = Phot_wa2.Eta();
    Phot_2_Phi  = Phot_wa2.Phi();
    EventNumber = data.GetEntriesFast();
    run_ = run;
    event_ = event;
    SCEta_ = SCEta_wa1;
    rho2012_ = rho2012R;
    phoSCRChIso1_ = phoSCRChIsoR1;
    phoSCRChIso2_ = phoSCRChIsoR2;
    pho1eleveto_ = pho1elevetoR;
    pho2eleveto_ = pho2elevetoR;
    deltaR1_ = deltaR1_R;
    deltaR2_ = deltaR2_R;

    outtree_->Fill();

  } // event loop

  Char_t fname[11];
  sprintf(fname,"evMuon_Zgg_%d.txt",MC_N);
  ofstream fout(fname);
  for(Int_t k = 0;k < rvector.size();k++){
    fout<<rvector[k]<<" "<<lvector[k]<<" "<<evector[k]<<" "<<normalY<<" "<<puvector[k]<<endl;
  }
  fout.close();

  fprintf(stderr, "Processed all events\n");
  
  fout_->cd();
  outtree_->Write();
  fout_->Close();
  
}

