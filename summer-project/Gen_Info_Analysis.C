#define baseAna_cxx

#include "baseAna.h"
#include <TH2.h>
#include <TStyle.h>

Bool_t baseAna::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.

  fReader.SetLocalEntry(entry);
  if(_data == 0){
    fReader_MC.SetLocalEntry(entry);
  }

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;

  nEvtRan++;                             //Total number of events containing everything (including the trash events).

  GoodEvt = true;
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents


    GenE.clear();
    GenMu.clear();
    for(int i=0; i<(int)*nGenPart; i++){
      if(GenPart_status[i] == 1){
	Lepton temp;
	temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
	temp.GenMomId = GenPart_genPartIdxMother[i];
	temp.PdgId = GenPart_pdgId[i];

	bool passCheck = temp.v.Pt() > 5.0 && fabs(temp.v.Eta()) < 2.4;
	
	if((GenPart_pdgId[i] == 13 || GenPart_pdgId[i] == -13) && passCheck){
	  GenMu.push_back(temp);
	}	
        if((GenPart_pdgId[i] == 11 || GenPart_pdgId[i] == -11) && passCheck){
	  GenE.push_back(temp);
	}
      }
    }
    SortPtGen(0);
    SortPtGen(1);

    Lepton motherMu;
    if((int)GenMu.size() > 0){
      for(int i=0; i < (int)GenMu.size(); i++){
	motherMu.GenMomId = GenPart_genPartIdxMother[GenMu.at(i).GenMomId];
	motherMu.PdgId = GenPart_pdgId[GenMu.at(i).GenMomId];
	
	while((GenMu.at(i).PdgId == motherMu.PdgId) || (GenMu.at(i).PdgId == (-motherMu.PdgId))){
	  motherMu.PdgId = GenPart_pdgId[motherMu.GenMomId];
	  motherMu.GenMomId = GenPart_genPartIdxMother[motherMu.GenMomId];
	}

	h.GenMomPdgId[0]->Fill(motherMu.PdgId);	
      }
    }

    Lepton motherE;
    if((int)GenE.size() > 0){
      for(int i=0; i < (int)GenE.size(); i++){
	motherE.GenMomId = GenPart_genPartIdxMother[GenE.at(i).GenMomId];
	motherE.PdgId = GenPart_pdgId[GenE.at(i).GenMomId];

	while((GenE.at(i).PdgId == motherE.PdgId) || (GenE.at(i).PdgId == (-motherE.PdgId))){
	  motherE.PdgId = GenPart_pdgId[motherE.GenMomId];
	  motherE.GenMomId = GenPart_genPartIdxMother[motherE.GenMomId];
	}

	h.GenMomPdgId[1]->Fill(motherE.PdgId);
      }
    }
    
    auto[dr_min_mu1, min_index_mu1] = FindMatchGentoGood(0,0); //Find match for gen_mu1
    auto[dr_min_mu2, min_index_mu2] = FindMatchGentoGood(0,1); //Find match for gen_mu2
    auto[dr_min_e1, min_index_e1] = FindMatchGentoGood(1,0); //Find match for gen_e1
    auto[dr_min_e2, min_index_e2] = FindMatchGentoGood(1,1); //Find match for gen_e2

    
    Lepton momMu1;
    for(int m=0; m < (int)goodMu.size(); m++){
      auto[dr_min_good_mu1, min_index_good_mu1] = FindMatchGoodtoGen(0,m);
      if(((int)GenMu.size() >= min_index_good_mu1) && (min_index_good_mu1 > 0)){
	momMu1.GenMomId = GenPart_genPartIdxMother[GenMu.at(min_index_good_mu1).GenMomId];
	momMu1.PdgId = GenPart_pdgId[GenMu.at(min_index_good_mu1).GenMomId];
	
	while((GenMu.at(min_index_good_mu1).PdgId == momMu1.PdgId) || (GenMu.at(min_index_good_mu1).PdgId == (-momMu1.PdgId))){
	  momMu1.PdgId = GenPart_pdgId[momMu1.GenMomId];
	  momMu1.GenMomId = GenPart_genPartIdxMother[momMu1.GenMomId];
	}	
      }
      
      if((int)goodMu.size() > 0){
	if(fabs(momMu1.PdgId) < 100){
	  h.Mu_IP[0]->Fill(Muon_dxy[m]);
	  h.Mu_IP[1]->Fill(Muon_dz[m]);
	  h.Mu_Iso[0]->Fill(Muon_pfRelIso04_all[m]); //Iso can be for R= 0.4 or 0.3, we are doing 0.4. It can also be for charged (just tracks in R=0.4) or all, we are doing for all.
	}
	if(fabs(momMu1.PdgId) > 100){
	  h.Mu_IP[2]->Fill(Muon_dxy[m]);
	  h.Mu_IP[3]->Fill(Muon_dz[m]);
	  h.Mu_Iso[1]->Fill(Muon_pfRelIso04_all[m]);
	}
      } 
    }
    
    Lepton momE1;
    for(int e=0; e < (int)goodE.size(); e++){
      auto[dr_min_good_e1, min_index_good_e1] = FindMatchGoodtoGen(1,e);
      if(((int)GenE.size() >= min_index_good_e1) && (min_index_good_e1 > 0)){
	momE1.GenMomId = GenPart_genPartIdxMother[GenE.at(min_index_good_e1).GenMomId];
	momE1.PdgId = GenPart_pdgId[GenE.at(min_index_good_e1).GenMomId];
	
	while((GenE.at(min_index_good_e1).PdgId == momE1.PdgId) || (GenE.at(min_index_good_e1).PdgId == (-momE1.PdgId))){
	  momE1.PdgId = GenPart_pdgId[momE1.GenMomId];
	  momE1.GenMomId = GenPart_genPartIdxMother[momE1.GenMomId];
	}	
      }

      if((int)goodE.size() > 0){
	if(fabs(momE1.PdgId) < 100){	
	  h.E_IP[0]->Fill(Electron_dxy[e]);
	  h.E_IP[1]->Fill(Electron_dz[e]);
	  h.E_Iso[0]->Fill(Electron_pfRelIso03_all[e]);
	}
	if(fabs(momE1.PdgId) > 100){	
	  h.E_IP[2]->Fill(Electron_dxy[e]);
	  h.E_IP[3]->Fill(Electron_dz[e]);
	  h.E_Iso[1]->Fill(Electron_pfRelIso03_all[e]);
	}
      } 
    }

    if((int)goodMu.size() > 0){
      float mu1_mT = transv_mass(goodMu.at(0).v.Pt(), *PuppiMET_pt, *PuppiMET_phi);
      h.mu_mT[0]->Fill(mu1_mT);
      if((int)goodMu.size() > 1){
	float mu2_mT = transv_mass(goodMu.at(1).v.Pt(), *PuppiMET_pt, *PuppiMET_phi);
	h.mu_mT[1]->Fill(mu2_mT);
      }
    }

    if((int)goodE.size() > 0){
      float e1_mT = transv_mass(goodE.at(0).v.Pt(), *PuppiMET_pt, *PuppiMET_phi);
      h.e_mT[0]->Fill(e1_mT);
      if((int)goodE.size() > 1){
	float e2_mT = transv_mass(goodE.at(1).v.Pt(), *PuppiMET_pt, *PuppiMET_phi);
	h.e_mT[1]->Fill(e2_mT);
      }
    }

      

void baseAna::SortPtGen(int opt)
{
  //This functions sorts an array in the decreasing order of pT for Gen arrays
  if(opt==0){
    for(int i=0; i<(int)GenMu.size()-1; i++){
      for(int j=i+1; j<(int)GenMu.size(); j++){
	if( GenMu[i].v.Pt() < GenMu[j].v.Pt() ) swap(GenMu.at(i),GenMu.at(j));
      }
    }
  }
  
  if(opt==1){
    for(int i=0; i<(int)GenE.size()-1; i++){
      for(int j=i+1; j<(int)GenE.size(); j++){
	if( GenE[i].v.Pt() < GenE[j].v.Pt() ) swap(GenE.at(i),GenE.at(j));
      }
    }
  }
  //Repeat this for the other arrays here.
}

std::tuple<float, int> baseAna::FindMatchGentoGood(int opt1, int opt2)
{
  float dr_min= 999;
  int min_index= 0;
  int temp_min = 0;
  
  if(opt1==0){
    if((int)GenMu.size()>opt2 && (int)goodMu.size()>0){
      for(unsigned int j=0; j < goodMu.size(); j++){
	float temp_dr = GenMu.at(opt2).v.DeltaR(goodMu.at(j).v);
	if(temp_dr < dr_min){
	  dr_min = temp_dr;
	  temp_min = j;
	}	
      }
    }
  }

  if(opt1==1){
    if((int)GenE.size()>opt2 && (int)goodE.size()>0){
      for(unsigned int j=0; j < goodE.size(); j++){
	float temp_dr = GenE.at(opt2).v.DeltaR(goodE.at(j).v);
	if(temp_dr < dr_min){
	  dr_min = temp_dr;
	  temp_min = j;
	}	
      }
    }
  }

  if(dr_min < 0.05){
    min_index = temp_min;
  }
  
  if(dr_min > 0.05){
    min_index = -1;
  }	
  
  return {dr_min, min_index};
}

std::tuple<float, int> baseAna::FindMatchGoodtoGen(int opt1, int opt2)
{
  float dr_min= 999;
  int min_index= 0;
  int temp_min = 0;
  
  if(opt1==0){
    if((int)goodMu.size()>opt2 && (int)GenMu.size()>0){
      for(unsigned int j=0; j < GenMu.size(); j++){
	float temp_dr = goodMu.at(opt2).v.DeltaR(GenMu.at(j).v);
	if(temp_dr < dr_min){
	  dr_min = temp_dr;
	  temp_min = j;
	}	
      }
    }
  }
  
  if(opt1==1){
    if((int)goodE.size()>opt2 && (int)GenE.size()>0){
      for(unsigned int j=0; j < GenE.size(); j++){
	float temp_dr = goodE.at(opt2).v.DeltaR(GenE.at(j).v);
	if(temp_dr < dr_min){
	  dr_min = temp_dr;
	  temp_min = j;
	}	
      }
    }
  }
  
  if(dr_min < 0.05){
    min_index = temp_min;
  }
  
  if(dr_min > 0.05){
    min_index = -1;
  }	
  
  return {dr_min, min_index};
}

   
int baseAna::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}

float baseAna::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float baseAna::transv_mass(float pT_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* pT_lep * MET *(1-cos(dphi)));
  return mT;
}

