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

  GoodEvt = true;                        //This is because I am working with a generated sample and not actual data, so every event is considered GoodEvt
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
    
    //This is the muon selection code chunk
    goodMu.clear();
    for(int i=0; i<(int)*nMuon; i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105);
      temp.id = -13*Muon_charge[i];
      temp.ind = i;

      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;

      if(passCuts){
	goodMu.push_back(temp);
      }     
    }
    SortPt(0);

    //This is the electron selection code
    goodE.clear();
    for(int i=0; i<(int)*nElectron; i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.05);
      temp.id = -11*Electron_charge[i];
      temp.ind = i;
      
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      passCuts = passCuts && Electron_pfRelIso03_all[i]<0.15;
      passCuts = passCuts && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1;
      passCuts = passCuts && Electron_cutBased[i] > 2;
      bool is_par_far = true;
      for(int j=0; j < (int)goodMu.size(); j++){
	float dR = temp.v.DeltaR(goodMu.at(j).v);
	if(dR < 0.4){
	  is_par_far = false;
	}
      }

      passCuts = passCuts && is_par_far;
      
      if(passCuts){
	goodE.push_back(temp);
      }      
    }
    SortPt(1);

    //This is the tau selection code (This was not necessary for my analysis, I wrote a code to test whether tau selection works, using the logic I used for mu and e) 
    goodTau.clear();
    for(int i=0; i<(int)*nTau; i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Tau_pt[i],Tau_eta[i],Tau_phi[i], 0.0000001);
      temp.id = 15*Tau_charge[i];
      temp.ind = i;
      
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      
      bool is_lep_far = true;
      for(unsigned int j=0; j < goodMu.size(); j++){
	if( temp.v.DeltaR(goodMu.at(j).v) < 0.4){
	  is_lep_far = false;
	  break;
	}
      }      
      if(is_lep_far){
	for(unsigned int j=0; j< goodE.size(); j++){;
	  if(temp.v.DeltaR(goodE.at(j).v) < 0.4){
	    is_lep_far = false;
	    break;
	  }
	}
      }
      
      passCuts = passCuts && is_lep_far;	
      
      if(passCuts){
	goodTau.push_back(temp);
      }      
    }
    SortPt(2);

   //This is the photon selection code
    goodPh.clear();
    for(int i=0; i<(int)*nPhoton; i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],Photon_phi[i], 0.0000001);
      temp.id = 22;
      temp.ind = i;
      
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      passCuts = passCuts && Photon_mvaID_WP80[i]==true;   
      
      bool is_lep_far = true;
      for(unsigned int j=0; j < goodMu.size(); j++){
	if( temp.v.DeltaR(goodMu.at(j).v) < 0.4){
	  is_lep_far = false;
	  break;
	}
      }      
      if(is_lep_far){
	for(unsigned int j=0; j< goodE.size(); j++){;
	  if(temp.v.DeltaR(goodE.at(j).v) < 0.4){
	    is_lep_far = false;
	    break;
	  }
	}
      }
      
      passCuts = passCuts && is_lep_far;	
      
      if(passCuts){
	goodPh.push_back(temp);
      }      
    }
    SortPt(3);

    //This is the jet selection code. I made separate arrays for B-tagged jets to help with my analysis. 
    goodJet.clear();
    goodBJet.clear();
    NoBJets.clear();
    for(int i=0; i<(int)*nJet; i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], 0.0000001);
      temp.ind = i;
      
      bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4;
      
      bool is_lep_far = true;
      for(unsigned int j=0; j < goodMu.size(); j++){
	if( temp.v.DeltaR(goodMu.at(j).v) < 0.4){
	  is_lep_far = false;
	  break;
	}
      }      
      if(is_lep_far){
	for(unsigned int j=0; j< goodE.size(); j++){;
	  if(temp.v.DeltaR(goodE.at(j).v) < 0.4){
	    is_lep_far = false;
	    break;
	  }
	}
      }
      if(is_lep_far){
	for(unsigned int j=0; j< goodPh.size(); j++){;
	  if(temp.v.DeltaR(goodPh.at(j).v) < 0.4){
	    is_lep_far = false;
	    break;
	  }
	}
      }
      
      passCuts = passCuts && is_lep_far;	
      
      if(passCuts){
	goodJet.push_back(temp);
      }

      bool passBTag = passCuts && Jet_btagDeepFlavB[i] > 0.7;
      bool notpassBTag = passCuts && Jet_btagDeepFlavB[i] < 0.7;

      if(passBTag){
	goodBJet.push_back(temp);
      }
      if(notpassBTag){
	NoBJets.push_back(temp);
      }
    }
    SortPt(4);
    SortPt(5);
    SortPt(6);

    //Invariant mass declaration code. I made this to quickly access the IVMs I needed frequently during the analysis.
    float M_yy = 0;
    if((int)goodPh.size() > 1){
       M_yy = (goodPh.at(0).v + goodPh.at(1).v).M();
    }

    float M_mumu = 0;
    if((int)goodMu.size() > 1){
      M_mumu = (goodMu.at(0).v + goodMu.at(1).v).M();
    }
    float M_ee = 0;
    if((int)goodE.size() > 1){
      M_ee = (goodE.at(0).v + goodE.at(1).v).M();
    }
    float M_jet = 0;
    if((int)goodJet.size() > 1){
      M_jet = (goodJet.at(0).v + goodJet.at(1).v).M();
    }
    float M_nobjet = 0;
    if((int)NoBJets.size() > 1){
      M_nobjet = (NoBJets.at(0).v + NoBJets.at(1).v).M();
    }

    float mT_mumu = 0;
    if((int)goodMu.size() > 0){
      mT_mumu = transv_mass(goodMu.at(0).v.E(), *PuppiMET_pt, (goodMu.at(0).v.Phi() - *PuppiMET_phi));
    }
    float mT_ee = 0;
    if((int)goodE.size() > 0){
      mT_ee = transv_mass(goodE.at(0).v.E(), *PuppiMET_pt, (goodE.at(0).v.Phi() - *PuppiMET_phi));
    }

//========================================================== THIS IS WHERE THE PLOTS START =============================================================================================================================
//**********************************************************************************************************************************************************************************************************************  
//========================================================== This is where I made the various plots to analyse the simulation information ==============================================================================

//The plots filled in the following order: i) Plots of Mu properties, E properties, Tau properties, Y properties, jet properties
//                                         ii) Plots I made to analyse Sample A, B, C, D 

    //===============================================================================================================================
    //                                                M U O N S
    //===============================================================================================================================

    h.nmuon->Fill((int)goodMu.size());
    if((int)goodMu.size()>0){
      h.muprop[0]->Fill(goodMu.at(0).v.Pt());
      h.muprop[1]->Fill(goodMu.at(0).v.Eta());
      h.muprop[2]->Fill(goodMu.at(0).v.Phi());
      h.muprop[3]->Fill(goodMu.at(0).id);
      if ((int)goodMu.size()>1){
	h.muprop2[0]->Fill(goodMu.at(1).v.Pt());
	h.muprop2[1]->Fill(goodMu.at(1).v.Eta());
	h.muprop2[2]->Fill(goodMu.at(1).v.Phi());
	h.muprop2[3]->Fill(goodMu.at(1).id);

	h.mudPhi->Fill((goodMu.at(0).v.DeltaPhi(goodMu.at(1).v)));
	h.mudR->Fill((goodMu.at(0).v.DeltaR(goodMu.at(1).v)));

	h.muInvM[0]->Fill((goodMu.at(0).v + goodMu.at(1).v).M());
	if(goodMu.at(0).v.Pt() > 30 && goodMu.at(1).v.Pt() > 10){
	  h.muInvM[1]->Fill((goodMu.at(0).v + goodMu.at(1).v).M());
	}
	if(goodMu.at(0).v.Pt() > 60 && goodMu.at(1).v.Pt() > 10){
	  h.muInvM[2]->Fill((goodMu.at(0).v + goodMu.at(1).v).M());
	}
	if(goodMu.at(0).v.Pt() >35 && goodMu.at(1).v.Pt() > 25){
	  h.muInvM[3]->Fill((goodMu.at(0).v + goodMu.at(1).v).M());
	}
      }
    }

    //===============================================================================================================================
    //                                                E L E C T R O N S
    //===============================================================================================================================

    h.nelec->Fill((int)goodE.size());
    if((int)goodE.size()>0){
      h.eprop[0]->Fill(goodE.at(0).v.Pt());
      h.eprop[1]->Fill(goodE.at(0).v.Eta());
      h.eprop[2]->Fill(goodE.at(0).v.Phi());
      h.eprop[3]->Fill(goodE.at(0).id);
      if ((int)goodE.size()>1){
	h.eprop2[0]->Fill(goodE.at(1).v.Pt());
	h.eprop2[1]->Fill(goodE.at(1).v.Eta());
	h.eprop2[2]->Fill(goodE.at(1).v.Phi());
	h.eprop2[3]->Fill(goodE.at(1).id);

	h.edPhi->Fill((goodE.at(0).v.DeltaPhi(goodE.at(1).v)));
	h.edR->Fill((goodE.at(0).v.DeltaR(goodE.at(1).v)));

	h.eInvM[0]->Fill((goodE.at(0).v + goodE.at(1).v).M());
	if(goodE.at(0).v.Pt() > 30 && goodE.at(1).v.Pt() > 10){
	  h.eInvM[1]->Fill((goodE.at(0).v + goodE.at(1).v).M());
	}
	if(goodE.at(0).v.Pt() > 60 && goodE.at(1).v.Pt() > 10){
	  h.eInvM[2]->Fill((goodE.at(0).v + goodE.at(1).v).M());
	}
	if(goodE.at(0).v.Pt() >35 && goodE.at(1).v.Pt() > 25){
	  h.eInvM[3]->Fill((goodE.at(0).v + goodE.at(1).v).M());
	}

	h.dielePt->Fill((goodE.at(0).v + goodE.at(1).v).Pt());
      }
    }
    

    /*   if((int)goodMu.size() > 0 && (int)goodE.size() > 0){
      h.emuIVM->Fill((goodMu.at(0).v + goodE.at(0).v).M());
      }*/

    //===============================================================================================================================
    //                                                 T A U S
    //===============================================================================================================================

    h.ntau->Fill((int)goodTau.size());
    if((int)goodTau.size() > 0){
      h.tauprop[0]->Fill(goodTau.at(0).v.Eta());
      h.tauprop[1]->Fill(goodTau.at(0).v.Phi());
      h.tauprop[2]->Fill(goodTau.at(0).v.Pt());
      if ((int)goodTau.size()>1){
	h.tauInvM->Fill((goodTau.at(0).v + goodTau.at(1).v).M());
	h.taudPhi->Fill((goodTau.at(0).v.DeltaPhi(goodTau.at(1).v)));
	h.taudR->Fill((goodTau.at(0).v.DeltaR(goodTau.at(1).v)));
      }
    }
    

    //===============================================================================================================================
    //                                            P H O T O N S
    //===============================================================================================================================

    h.nphoton->Fill((int)goodPh.size());
    if((int)goodPh.size() > 0){
      h.phprop[0]->Fill(goodPh.at(0).v.Eta());
      h.phprop[1]->Fill(goodPh.at(0).v.Phi());
      h.phprop[2]->Fill(goodPh.at(0).v.Pt());
      if ((int)goodPh.size()>1){
	h.phInvM->Fill((goodPh.at(0).v + goodPh.at(1).v).M());
	h.phdPhi->Fill((goodPh.at(0).v.DeltaPhi(goodPh.at(1).v)));
	h.phdR->Fill((goodPh.at(0).v.DeltaR(goodPh.at(1).v)));
      }
    }
    
    //===============================================================================================================================
    //                                               J E T S
    //===============================================================================================================================

    h.njets->Fill((int)goodJet.size());
    if((int)goodJet.size() > 0){
      h.jetprop[0]->Fill(goodJet.at(0).v.Eta());
      h.jetprop[1]->Fill(goodJet.at(0).v.Phi());
      h.jetprop[2]->Fill(goodJet.at(0).v.Pt());
      if((int)goodJet.size() > 1){
	h.jetInvM[0]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
	h.jetdPhi->Fill((goodJet.at(0).v.DeltaPhi(goodJet.at(1).v)));
	h.jetdR->Fill((goodJet.at(0).v.DeltaR(goodJet.at(1).v)));
	h.jetprop[3]->Fill(goodJet.at(1).v.Pt());

	if(goodJet.at(0).v.Pt() > 45 && goodJet.at(1).v.Pt() > 35){
	  h.jetInvM[1]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
	}
	if(goodJet.at(0).v.Pt() > 60 && goodJet.at(1).v.Pt() > 35){
	  h.jetInvM[2]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
	}
	if(goodJet.at(0).v.Pt() > 55 && goodJet.at(1).v.Pt() > 45){
	  h.jetInvM[3]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
	}

	h.dijetPt->Fill((goodJet.at(0).v + goodJet.at(1).v).Pt());

	if((int)goodJet.size() > 2){
	  h.jetprop[4]->Fill(goodJet.at(2).v.Pt());
	  if((int)goodJet.size() > 3){
	    h.jetprop[5]->Fill(goodJet.at(3).v.Pt());
	  }
	}
      }
    }

    /*  for(int j=0; j < (int)goodJet.size(); j++){
      h.btagVal->Fill(Jet_btagDeepFlavB[goodJet.at(j).ind]);
      }*/

    h.nbjets->Fill((int)goodBJet.size());
    /*   if(goodBJet.size() > 1){
      h.bjetInvM->Fill((goodBJet.at(0).v + goodBJet.at(1).v).M());
      }*/

    //===============================================================================================================================
    //                                             M E T
    //===============================================================================================================================

    h.MET->Fill(*PuppiMET_pt);
    //  h.MET[1]->Fill(*PuppiMET_phi);

    //==================================================================================================================================================
    //                                                            SAMPLE C
    //==================================================================================================================================================

    if(100 < M_yy && M_yy < 140){
      h.nphoton_test->Fill((int)goodPh.size());
      h.njets_test->Fill((int)goodJet.size());
      if((int)goodJet.size() > 1){
	h.jetInvM_test->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
	h.jetdR_test->Fill((goodJet.at(0).v.DeltaR(goodJet.at(1).v)));
	h.jetdPhi_test->Fill((goodJet.at(0).v.DeltaPhi(goodJet.at(1).v)));
	h.jetdEta_test->Fill(goodJet.at(0).v.Eta() - goodJet.at(1).v.Eta());
	
	if((int)goodJet.size() > 2){
	  h.jetdR_j2j3_test->Fill((goodJet.at(1).v.DeltaR(goodJet.at(2).v)));
	  h.jetdPhi_j2j3_test->Fill((goodJet.at(1).v.DeltaPhi(goodJet.at(2).v)));
	  
	  h.higgs_check->Fill((goodJet.at(0).v + goodJet.at(1).v + goodJet.at(2).v).M());
	}
      }
      
      h.MET_test->Fill(*PuppiMET_pt);
    }
    
    //===================================================================================================================================================
    //                                                             SAMPLE D
    //===================================================================================================================================================
    
    if((70 < M_mumu && M_mumu < 106) || (70 < M_ee && M_ee < 106)){
      h.nmuon_test->Fill((int)goodMu.size());
      h.nelec_test->Fill((int)goodE.size());
      h.njets_test->Fill((int)goodJet.size());
      h.nbjets_test->Fill((int)goodBJet.size());
      h.nphoton_test->Fill((int)goodPh.size());
      
      if((int)goodE.size() > 1){
	h.eInvM_test->Fill((goodE.at(0).v + goodE.at(1).v).M());
      }
      
      if((int)goodE.size() > 0 && (int)goodMu.size() > 0){
	h.emuIVM_test->Fill((goodMu.at(0).v + goodE.at(0).v).M());
      }
      
      h.MET_test->Fill(*PuppiMET_pt);
      
      if((int)goodBJet.size() > 1){
	h.bjetdR_test->Fill((goodBJet.at(0).v.DeltaR(goodBJet.at(1).v)));
	h.bjetdPhi_test->Fill((goodBJet.at(0).v.DeltaPhi(goodBJet.at(1).v)));
      }
      
      if(0 < (int)goodBJet.size() && (int)goodBJet.size() < 3){
	h.nphoton_test_bjet[0]->Fill((int)goodPh.size());
      }
      
      h.nnobjets_test->Fill((int)NoBJets.size());
      
      if((int)NoBJets.size() > 1){
	h.only_jetInvM_1_mu[0]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
	if(NoBJets.at(0).v.Pt() > 60 && NoBJets.at(1).v.Pt() > 40){
	  h.only_jetInvM_1_mu[1]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
	}
	if(NoBJets.at(0).v.Pt() > 80 && NoBJets.at(1).v.Pt() > 40){
	  h.only_jetInvM_1_mu[3]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
	}	    
	if(NoBJets.at(0).v.Pt() >55 && NoBJets.at(1).v.Pt() > 45){
	  h.only_jetInvM_1_mu[2]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
	}

	h.di_nobjetPt_test->Fill((NoBJets.at(0).v + NoBJets.at(1).v).Pt());
      }
      
      if((int)goodBJet.size() > 0 && (int)NoBJets.size() > 0){
	h.test_IVM[0]->Fill((goodBJet.at(0).v + NoBJets.at(0).v).M());
      }
      if((int)NoBJets.size() > 1){
	h.test_IVM[1]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M() + *PuppiMET_pt);
      }
      if((int)goodBJet.size() > 1){
	h.test_IVM[2]->Fill((goodBJet.at(0).v + goodBJet.at(1).v).M() + *PuppiMET_pt);
      }
      if((int)goodBJet.size() > 0 && (int)NoBJets.size() > 0){
	h.test_IVM[3]->Fill((goodBJet.at(0).v + NoBJets.at(0).v).M() + *PuppiMET_pt);
      }
      if((int)goodBJet.size() > 0 && (int)NoBJets.size() > 0 && (70 < M_mumu && M_mumu < 106)){
	h.test_IVM[4]->Fill((goodBJet.at(0).v + NoBJets.at(0).v + goodMu.at(0).v + goodMu.at(1).v).M());
      }
      if((int)goodBJet.size() > 0 && (int)NoBJets.size() > 0 && (70 < M_ee && M_ee < 106)){
	h.test_IVM[5]->Fill((goodBJet.at(0).v + NoBJets.at(0).v + goodE.at(0).v + goodE.at(1).v).M());
      }
      if((int)goodBJet.size() > 1 && (int)NoBJets.size() > 1){
	h.test_IVM[6]->Fill((goodBJet.at(1).v + NoBJets.at(1).v).M());
      }
    }
    
    if(0 < (int)goodBJet.size() && (int)goodBJet.size() < 4){
      h.nphoton_test_bjet[1]->Fill((int)goodPh.size());
    }
    
    if((int)NoBJets.size() > 1){
      h.only_jetInvM_1[0]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
      if(NoBJets.at(0).v.Pt() > 60 && NoBJets.at(1).v.Pt() > 40){
	h.only_jetInvM_1[1]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
      }
      if(NoBJets.at(0).v.Pt() > 80 && NoBJets.at(1).v.Pt() > 40){
	h.only_jetInvM_1[3]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
      }	    
      if(NoBJets.at(0).v.Pt() >55 && NoBJets.at(1).v.Pt() > 45){
	h.only_jetInvM_1[2]->Fill((NoBJets.at(0).v + NoBJets.at(1).v).M());
      }	
      
      if((int)NoBJets.size() > 2){
	h.only_jetInvM_2->Fill((NoBJets.at(1).v + NoBJets.at(2).v).M());
      }
    }
    
    if((int)goodMu.size() > 1 && (int)goodBJet.size() > 1){
      h.BjetLepInvM[0]->Fill((goodMu.at(0).v + goodMu.at(1).v + goodBJet.at(0).v + goodBJet.at(0).v).M());
    }
    if((int)goodE.size() > 1 && (int)goodBJet.size() > 1){
      h.BjetLepInvM[1]->Fill((goodE.at(0).v + goodE.at(1).v + goodBJet.at(0).v + goodBJet.at(0).v).M());
    }
    
    h.tot_lep->Fill((int)goodMu.size() + (int)goodE.size());

    if((int)goodJet.size() == 4){
      h.nphoton_test_j4->Fill((int)goodPh.size());
    }

    if((int)goodMu.size() == 0){
      h.nelec_test_mu[0]->Fill((int)goodE.size());
    }
    if((int)goodMu.size() == 1){
      h.nelec_test_mu[1]->Fill((int)goodE.size());
    }
    if((int)goodMu.size() == 2){
      h.nelec_test_mu[2]->Fill((int)goodE.size());
    }

    if((int)goodE.size() == 0){
      h.nmuon_test_e[0]->Fill((int)goodMu.size());
    }
    if((int)goodE.size() == 1){
      h.nmuon_test_e[1]->Fill((int)goodMu.size());
    }
    if((int)goodE.size() == 2){
      h.nmuon_test_e[2]->Fill((int)goodMu.size());
    }

    if(66 < M_ee && M_ee < 116){
      h.dielePt_test->Fill((goodE.at(0).v + goodE.at(1).v).Pt());
    }

    if((66 < M_ee && M_ee < 116) || (66 < M_mumu && M_mumu < 116)){
      if((int)NoBJets.size() > 1 && (int)goodBJet.size() > 0){
	h.trijetIVM[0]-> Fill((NoBJets.at(0).v + NoBJets.at(1).v + goodBJet.at(0).v).M());
      }
      if((int)goodJet.size() > 2){
	h.trijetIVM[1]-> Fill((goodJet.at(0).v + goodJet.at(1).v + goodJet.at(2).v).M());
      }

      if(65 < M_nobjet && M_nobjet < 95){
	if((int)goodBJet.size() > 0){
	  h.trijetIVM[2]-> Fill((NoBJets.at(0).v + NoBJets.at(1).v + goodBJet.at(0).v).M());
	}
      }
    }
   
//=============================================================================================================================
//                                              SAMPLE B
//============================================================================================================================= 

    if((int)goodJet.size() > 1 && goodMu.size() > 0 && goodPh.size() > 0) {
      h.test_IVM_B->Fill((goodJet.at(0).v + goodJet.at(1).v + goodPh.at(0).v + goodMu.at(0).v).M());
    }
    
    if(60 < M_jet && M_jet < 100){
      if((int)goodMu.size() == 0){
	h.nelec_test[0]->Fill((int)goodE.size());
      }
      if((int)goodMu.size() == 1){
	h.nelec_test[1]->Fill((int)goodE.size());
      }
      if((int)goodMu.size() == 2){
	h.nelec_test[2]->Fill((int)goodE.size());
      }
      
      h.nphoton_test[0]->Fill((int)goodPh.size());
      
      if((int)goodMu.size() > 0){
	h.MuTransMass[0]->Fill(transv_mass(goodMu.at(0).v.E(), *PuppiMET_pt, (goodMu.at(0).v.Phi() - *PuppiMET_phi)));
	if((int)goodE.size() > 0){
	  h.eTransMass[0]->Fill(transv_mass(goodE.at(0).v.E(), *PuppiMET_pt, (goodE.at(0).v.Phi() - *PuppiMET_phi)));
	}
      }
      if((int)goodE.size() > 0){
	h.eTransMass[1]->Fill(transv_mass(goodE.at(0).v.E(), *PuppiMET_pt, (goodE.at(0).v.Phi() - *PuppiMET_phi)));
	if((int)goodMu.size() > 0){
	  h.MuTransMass[1]->Fill(transv_mass(goodMu.at(0).v.E(), *PuppiMET_pt, (goodMu.at(0).v.Phi() - *PuppiMET_phi)));
	}
      }
    }
    
    if((int)goodMu.size() == 0){
      h.njet_test[0]->Fill((int)goodJet.size());
      h.nphoton_test[1]->Fill((int)goodPh.size());
    }
    
    if((int)goodMu.size() == 1){
      h.njet_test[1]->Fill((int)goodJet.size());
      h.nphoton_test[3]->Fill((int)goodPh.size());
    }
    
    if((int)goodE.size() == 0){
      h.njet_test[2]->Fill((int)goodJet.size());
      h.nphoton_test[2]->Fill((int)goodPh.size());
    }
    
    if((int)goodE.size() == 1){
      h.njet_test[3]->Fill((int)goodJet.size());
      h.nphoton_test[4]->Fill((int)goodPh.size());
    }
    
    if((int)goodJet.size() == 1){
      h.nphoton_test[5]->Fill((int)goodPh.size());
    }
    if((int)goodJet.size() == 2){
      h.nphoton_test[6]->Fill((int)goodPh.size());
    }

    if((int)goodMu.size() > 0){
      h.MuTransMass[2]->Fill(transv_mass(goodMu.at(0).v.E(), *PuppiMET_pt, (goodMu.at(0).v.Phi() - *PuppiMET_phi)));
    }
    if((int)goodE.size() > 0){
      h.eTransMass[2]->Fill(transv_mass(goodE.at(0).v.E(), *PuppiMET_pt, (goodE.at(0).v.Phi() - *PuppiMET_phi)));
    }

    if(60 < mT_ee && mT_ee < 100){
      if((int)goodJet.size() > 1){
	h.mjj_test[0]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
      }
    }
    if(60 < mT_mumu && mT_mumu < 100){
      if((int)goodJet.size() > 1){
	h.mjj_test[1]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
      }
    }

    //==================================================================================================================================
    //                                                SAMPLE A
    //==================================================================================================================================

    if(84 < M_ee && M_ee < 105){
      h.nelec_test[0]->Fill((int)goodE.size());
      h.nmuon_test[0]->Fill((int)goodMu.size());
      h.nphoton_test[0]->Fill((int)goodPh.size());
      h.njet_test[0]->Fill((int)goodJet.size());
      
      h.edPhi_test[0]->Fill((goodE.at(0).v.DeltaPhi(goodE.at(1).v)));
      h.edR_test[0]->Fill((goodE.at(0).v.DeltaR(goodE.at(1).v)));
    }
    
    if(84 < M_mumu && M_mumu < 105){
      h.nelec_test[2]->Fill((int)goodE.size());
      h.nmuon_test[2]->Fill((int)goodMu.size());
      h.nphoton_test[2]->Fill((int)goodPh.size());
      h.njet_test[2]->Fill((int)goodJet.size());
      
      h.mudPhi_test[0]->Fill((goodMu.at(0).v.DeltaPhi(goodMu.at(1).v)));
      h.mudR_test[0]->Fill((goodMu.at(0).v.DeltaR(goodMu.at(1).v)));
    }
    
    if(50 < M_ee && M_ee < 85){
      h.nelec_test[1]->Fill((int)goodE.size());
      h.nmuon_test[1]->Fill((int)goodMu.size());
      h.nphoton_test[1]->Fill((int)goodPh.size());
      h.njet_test[1]->Fill((int)goodJet.size());
      
      h.edPhi_test[1]->Fill((goodE.at(0).v.DeltaPhi(goodE.at(1).v)));
      h.edR_test[1]->Fill((goodE.at(0).v.DeltaR(goodE.at(1).v)));
    }
    
    if(50 < M_mumu && M_mumu < 85){
      h.nelec_test[3]->Fill((int)goodE.size());
      h.nmuon_test[3]->Fill((int)goodMu.size());
      h.nphoton_test[3]->Fill((int)goodPh.size());
      h.njet_test[3]->Fill((int)goodJet.size());
      
      h.mudPhi_test[1]->Fill((goodMu.at(0).v.DeltaPhi(goodMu.at(1).v)));
      h.mudR_test[1]->Fill((goodMu.at(0).v.DeltaR(goodMu.at(1).v)));
    }
    
    if((84 < M_mumu && M_mumu < 105) || (84 < M_ee && M_ee < 105)){
      if((int)goodJet.size() == 0){
	h.nphoton_test[4]->Fill((int)goodPh.size());
      }
      if((int)goodJet.size() == 1){
	h.nphoton_test[5]->Fill((int)goodPh.size());
      }
      if((int)goodJet.size() == 2){
	h.nphoton_test[6]->Fill((int)goodPh.size());
      }

      if(goodJet.size() > 1){
	h.jetInvM_test[0]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
      }
    }

    if((50 < M_mumu && M_mumu < 85) || (50 < M_ee && M_ee < 85)){
      if(goodJet.size() > 1){
	h.jetInvM_test[1]->Fill((goodJet.at(0).v + goodJet.at(1).v).M());
      }
    }

    if((int)goodPh.size() > 0){
      if((int)goodE.size() > 1){
	h.m_eey[0]->Fill((goodE.at(0).v + goodE.at(1).v + goodPh.at(0).v).M());
      }
      if(50 < M_ee && M_ee < 85){
	h.m_eey[1]->Fill((goodE.at(0).v + goodE.at(1).v + goodPh.at(0).v).M());
      }
      if(84 < M_ee && M_ee < 105){
	h.m_eey[2]->Fill((goodE.at(0).v + goodE.at(1).v + goodPh.at(0).v).M());
      }

      if((int)goodMu.size() > 1){
	h.m_mumuy[0]->Fill((goodMu.at(0).v + goodMu.at(1).v + goodPh.at(0).v).M());
      }
      if(50 < M_mumu && M_mumu < 85){
	h.m_mumuy[1]->Fill((goodMu.at(0).v + goodMu.at(1).v + goodPh.at(0).v).M());
      }
      if(84 < M_mumu && M_mumu < 105){
	h.m_mumuy[2]->Fill((goodMu.at(0).v + goodMu.at(1).v + goodPh.at(0).v).M());
      }
    }
    
  }//GoodEvt  
  return kTRUE;
}

//============================================================================================================================================================
// ============================================ END OF GOOD EVENTS FUNCTIONS =================================================================================
//============================================================================================================================================================

//######################################
//        USER DEFINED FUNCTIONS
//######################################
void baseAna::SortPt(int opt)
{
  //This functions sorts an array in the decreasing order of pT.
  //For recoMu:
  if(opt==0){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j));
      }
    }
  }

  if(opt==1){
    for(int i=0; i<(int)goodE.size()-1; i++){
      for(int j=i+1; j<(int)goodE.size(); j++){
	if( goodE[i].v.Pt() < goodE[j].v.Pt() ) swap(goodE.at(i),goodE.at(j));
      }
    }
  }

  if(opt==2){
    for(int i=0; i<(int)goodTau.size()-1; i++){
      for(int j=i+1; j<(int)goodTau.size(); j++){
	if( goodTau[i].v.Pt() < goodTau[j].v.Pt() ) swap(goodTau.at(i),goodTau.at(j));
      }
    }
  }

  if(opt==3){
    for(int i=0; i<(int)goodPh.size()-1; i++){
      for(int j=i+1; j<(int)goodPh.size(); j++){
	if( goodPh[i].v.Pt() < goodPh[j].v.Pt() ) swap(goodPh.at(i),goodPh.at(j));
      }
    }
  }

  if(opt==4){
    for(int i=0; i<(int)goodJet.size()-1; i++){
      for(int j=i+1; j<(int)goodJet.size(); j++){
	if( goodJet[i].v.Pt() < goodJet[j].v.Pt() ) swap(goodJet.at(i),goodJet.at(j));
      }
    }
  }

  if(opt==5){
    for(int i=0; i<(int)goodBJet.size()-1; i++){
      for(int j=i+1; j<(int)goodBJet.size(); j++){
	if( goodBJet[i].v.Pt() < goodBJet[j].v.Pt() ) swap(goodBJet.at(i),goodBJet.at(j));
      }
    }
  }

  if(opt==6){
    for(int i=0; i<(int)NoBJets.size()-1; i++){
      for(int j=i+1; j<(int)NoBJets.size(); j++){
	if( NoBJets[i].v.Pt() < NoBJets[j].v.Pt() ) swap(NoBJets.at(i),NoBJets.at(j));
      }
    }
  }

  //Repeat this for the other arrays here.
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

float baseAna::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}

void baseAna::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.
  //The histrograms are booked in the following order: i) Hist. of Mu properties, E properties, Tau properties, Y properties, jet properties
  //                                                   ii) Hist. I made to analyse Sample A, B, C, D   
  //Example : new TH1F ("hst_name", "hst title", total bins, initial x, final x);

  //===============================================================================================================================
  //                                   M U O N S
  //===============================================================================================================================
  h.nmuon = new TH1F("nmuon","NMuons",10,0,10);
  
  h.muprop[0] = new TH1F("mu1_pt","Leading muon pT",200,0,200);
  h.muprop[1] = new TH1F("mu1_eta","Leading muon Eta",120,-3,3);
  h.muprop[2] = new TH1F("mu1_phi", "Leading muon Phi", 120, -3, 3);
  h.muprop[3] = new TH1F("mu1_charge", "Leading muon charge", 40, -20, 20);
  h.muprop2[0] = new TH1F("mu2_pt","SubLeading muon pT",200,0,200);
  h.muprop2[1] = new TH1F("mu2_eta","SubLeading muon Eta",120,-3,3);
  h.muprop2[2] = new TH1F("mu2_phi", "SubLeading muon Phi", 120, -3, 3);
  h.muprop2[3] = new TH1F("mu2_charge", "SubLeading muon charge", 40, 20, 20);
  
  h.muInvM[0] = new TH1F("M_mumu", "M_mumu", 200, 0, 200);
  h.muInvM[1] = new TH1F("M_mumu_con1", "M_mumu (Mu1>30, Mu2>10)", 200, 0, 200);
  h.muInvM[2] = new TH1F("M_mumu_con2", "M_mumu (Mu1>60, Mu2>10)", 200, 0, 200);
  h.muInvM[3] = new TH1F("M_mumu_con3", "M_mumu (Mu1>35, Mu2>25)", 200, 0, 200);

  h.mudPhi = new TH1F("mudPhi", "mu1 mu2 dPhi", 120, -3, 3);
  h.mudR = new TH1F("mudR", "mu1 mu2 dR", 120, 0, 10);

  //===============================================================================================================================
  //                                   E L E C T R O N
  //===============================================================================================================================

  h.nelec = new TH1F("nelec","NELectron",10,0,10);
  
  h.eprop[0] = new TH1F("e1_pt","Leading elec pT",200,0,200);
  h.eprop[1] = new TH1F("e1_eta","Leading elec Eta",120,-3,3);
  h.eprop[2] = new TH1F("e1_phi", "Leading elec Phi", 120, -3, 3);
  h.eprop[3] = new TH1F("e1_charge", "Leading elec charge", 40, -20, 20);

  h.eprop2[0] = new TH1F("e2_pt","SubLeading elec pT",200,0,200);
  h.eprop2[1] = new TH1F("e2_eta","SubLeading elec Eta",120,-3,3);
  h.eprop2[2] = new TH1F("e2_phi", "SubLeading elec Phi", 120, -3, 3);
  h.eprop2[3] = new TH1F("e2_charge", "SubLeading elec charge", 40, -20, 20);

  h.eInvM[0] = new TH1F("M_ee", "M_ee", 200, 0, 200);
  h.eInvM[1] = new TH1F("M_ee_con1", "M_ee, e1>30 e2>10", 200, 0, 200);
  h.eInvM[2] = new TH1F("M_ee_con2", "M_ee; e1>60 e2>10", 200, 0, 200);
  h.eInvM[3] = new TH1F("M_ee_con3", "M_ee, e1>35 e2>25", 200, 0, 200);

  h.dielePt = new TH1F("Pt_ee", "Pt_ee", 400, 0, 400);

  h.edPhi = new TH1F("edPhi", "e1 e2 dPhi", 120, -3, 3);
  h.edR = new TH1F("edR", "e1 e2 dR", 120, 0, 10);

  // h.emuIVM = new TH1F("M_mue", "M_emu", 200, 0, 200);

  //===============================================================================================================================
  //                                   T A U S
  //===============================================================================================================================

  h.ntau = new TH1F("ntau","NTau",10,0,10);
  
  h.tauInvM = new TH1F("M_tt", "M_tt", 200, 0, 200);
  
  h.tauprop[0] = new TH1F("tau1_eta","Leading tau Eta",120,-3,3);
  h.tauprop[1] = new TH1F("tau1_phi", "Leading tau  Phi", 120, -3, 3);
  h.tauprop[2] = new TH1F("tau1_pt", "Leading tau Pt", 200, 0, 200);
  
  h.taudPhi = new TH1F("taudPhi", "t1 t2 dPhi", 120, -3, 3);
  h.taudR = new TH1F("taudR", "t1 t2 dR", 120, 0, 10);


  //===============================================================================================================================
  //                                   P H O T O N S
  //===============================================================================================================================
  h.nphoton = new TH1F("ngamma","NPhoton",10,0,10);
  
  h.phInvM = new TH1F("M_gamgam", "M_yy", 200, 0, 200);
  
  h.phprop[0] = new TH1F("ph1_eta","Leading photon Eta",120,-3,3);
  h.phprop[1] = new TH1F("ph1_phi", "Leading photon Phi", 120, -3, 3);
  h.phprop[2] = new TH1F("ph1_pt", "Leading photon Pt", 200, 0, 200);
  
  h.phdPhi = new TH1F("phdPhi", "Ph1 Ph2 dPhi", 120, -3, 3);
  h.phdR = new TH1F("phdR", "ph1 ph2 dR", 120, 0, 10);

  //===============================================================================================================================
  //                                   J E T S
  //===============================================================================================================================

  h.njets =  new TH1F("njets","NJets",10,0,10);
  h.nbjets =  new TH1F("nbjets","NBJets",10,0,10);

  h.jetprop[0] = new TH1F("j1_eta","Leading jet Eta",120,-3,3);
  h.jetprop[1] = new TH1F("j1_phi", "Leading jet Phi", 120, -3, 3);
  h.jetprop[2] = new TH1F("j1_pt","Leading jet pT",200,0,200);
  h.jetprop[3] = new TH1F("j2_pt","SubLeading jet pT",200,0,200);
  h.jetprop[4] = new TH1F("j3_pt","j3  pT",200,0,200);
  h.jetprop[5] = new TH1F("j4_pt","j4  pT",200,0,200);
  
  h.jetdPhi = new TH1F("jetdPhi", "j1 j2 dPhi", 120, -3, 3);
  h.jetdR = new TH1F("jetdR", "j1 j2 dR", 120, 0, 10);
  
  h.jetInvM[0] = new TH1F("M_jet", "M_jet", 200, 0, 200);
  h.jetInvM[1] = new TH1F("M_jet_1", "M_jet j1>45 j2>35", 200, 0, 200);
  h.jetInvM[2] = new TH1F("M_jet_2", "M_jet j1>60 j2>35", 200, 0, 200);
  h.jetInvM[3] = new TH1F("M_jet_3", "M_jet j1>55 j2>45", 200, 0, 200);
  h.dijetPt = new TH1F("Pt_jet", "Pt_jet", 400, 0, 400);
  h.bjetInvM = new TH1F("M_bjet", "M_bjet", 200, 0, 200);
  
  // h.btagVal = new TH1F("btagVal", "B-tag Value", 300, -1, 2);

  //===============================================================================================================================
  //                                   M E T
  //===============================================================================================================================

  h.MET = new TH1F("MET_pt", "MET Pt", 200, 0, 200);
  //  h.MET[1] = new TH1F("MET_phi", "MET Phi", 64,-3, 3);

  //===============================================================================================================================
  //                                             TEST PLOTS 
  //===============================================================================================================================

  //===============================================================================================================================
  //                                              SAMPLE C
  //===============================================================================================================================

  h.nphoton_test=  new TH1F("ngamma_test","NPhoton (100<Myy<140)",10,0,10); 
  h.jetInvM_test = new TH1F("M_jet_test", "M_jet (100<Myy<140)", 200, 0, 200);
  h.njets_test =  new TH1F("njets_test","NJets (100<Myy<140)",10,0,10);
  h.jetdPhi_test = new TH1F("jetdPhi_test", "j1 j2 dPhi (100<Myy<140)", 120, -3, 3);
  h.jetdEta_test = new TH1F("jetdEta_test", "j1 j2 dEta (100<Myy<140)", 120, -3, 3);
  h.jetdR_test = new TH1F("jetdR_test", "j1 j2 dR (100<Myy<140)", 120, 0, 10); 
  
  h.MET_test = new TH1F("MET_pt_test", "MET Pt (100<Myy<140)", 200, 0, 200);
  
  h.jetdPhi_j2j3_test = new TH1F("jetdPhi_test_j2j3", "j3 j2 dPhi (100<Myy<140)", 120, -3, 3);
  h.jetdR_j2j3_test = new TH1F("jetdR_test_j2j3", "j3 j2 dR (100<Myy<140)", 120, 0, 10);
  
  h.higgs_check = new TH1F("higgs_check", "IVM of j1+j2+j3", 200, 0, 200); 

  //===============================================================================================================================
  //                                          SAMPLE D
  //===============================================================================================================================

  h.nmuon_test = new TH1F("nmuon_test","NMuon (70<M_mumu<110)",10,0,10);
  h.nelec_test = new TH1F("nelec_test","NELectron (70<M_mumu<110)",10,0,10);
  h.eInvM_test = new TH1F("M_ee_test", "M_ee (70<M_mumu<110)", 200, 0, 200);
  h.dielePt_test = new TH1F("Pt_ee_test", "Pt_ee_test (66 < M_ee < 116)", 400, 0, 400);
  h.emuIVM_test = new TH1F("M_mue_test", "M_emu (70<M_mumu<110)", 200, 0, 200);
  h.njets_test =  new TH1F("njets_test","NJets (70<Mmumu<110)",10,0,10);
  h.nbjets_test =  new TH1F("nbjets_test","NBJets (70<Mmumu<110)",10,0,10);
  h.nnobjets_test =  new TH1F("nnobjets_test","Only No B-Jets (70<Mmumu<110)",10,0,10);
  h.bjetdPhi_test = new TH1F("bjetdPhi_test", "j1 j2 dPhi B (70<Mmm<110)", 120, -3, 3);
  h.bjetdR_test = new TH1F("bjetdR_test", "j1 j2 dR B (70<Mmm<110)", 120, 0, 10);
  h.nphoton_test = new TH1F("ngamma_test","NPhoton (70<Mmm<110)",10,0,10);
  h.nphoton_test_bjet[0] = new TH1F("ngamma_test_bjet1","NPhoton 0<NBJets<3  (70<Mmm<110)",10,0,10);
  h.nphoton_test_bjet[1] = new TH1F("ngamma_test_bjet2","NPhoton 0<NBJets<4",10,0,10);
  h.MET_test = new TH1F("MET_pt_test", "MET Pt (70<Mmll<106)", 200, 0, 200);

  h.only_jetInvM_1[0] = new TH1F("M_only_jet1", "M_only_jet j1j2", 200, 0, 200);
  h.only_jetInvM_1[1] = new TH1F("M_only_jet1_con1", "M_only_jet j1j2 j1>60,j2>40", 200, 0, 200);
  h.only_jetInvM_1[2] = new TH1F("M_only_jet1_con2", "M_only_jet j1j2 j1>55,j2>45", 200, 0, 200);
  h.only_jetInvM_1[3] = new TH1F("M_only_jet1_con3", "M_only_jet j1j2 j1>80,j2>40", 200, 0, 200);

  h.only_jetInvM_1_mu[0] = new TH1F("M_only_jet1_mu", "M_only_jet j1j2 (70<Mmm<110)", 200, 0, 200);
  h.only_jetInvM_1_mu[1] = new TH1F("M_only_jet1_mu_con1", "M_only_jet j1j2 j1>60,j2>40  (70<Mmm<110)", 200, 0, 200);
  h.only_jetInvM_1_mu[2] = new TH1F("M_only_jet1_mu_con2", "M_only_jet j1j2 j1>55,j2>45  (70<Mmm<110)", 200, 0, 200);
  h.only_jetInvM_1_mu[3] = new TH1F("M_only_jet1_mu_con3", "M_only_jet j1j2 j1>80,j2>40 (70<Mmm<110)", 200, 0, 200);
  h.di_nobjetPt_test = new TH1F("Pt_nobjet", "Pt_nobjet (70<Mmm<110)", 400, 0, 400);

  h.only_jetInvM_2 = new TH1F("M_only_jet2", "M_only_jet j2j3", 200, 0, 200);

  h.BjetLepInvM[0] = new TH1F("M_bjet_mu", "M_bjet_mu", 200, 0, 200);
  h.BjetLepInvM[1] = new TH1F("M_bjet_e", "M_bjet_e", 200, 0, 200);

  h.tot_lep =  new TH1F("nlep","NLepton",10,0,10);

  h.nphoton_test_j4 = new TH1F("ngamma_test_j4","NPhoton: NJets = 4",10,0,10);
  h.nmuon_test_e[0] =  new TH1F("nmuon_test_e1","NMuon: NElec = 0",10,0,10);
  h.nmuon_test_e[1] =  new TH1F("nmuon_test_e2","NMuon: NElec = 1",10,0,10);
  h.nmuon_test_e[2] =  new TH1F("nmuon_test_e3","NMuon: NElec = 2",10,0,10);
  h.nelec_test_mu[0] =  new TH1F("nelec_test_mu1","NElec: NMuon = 0",10,0,10);
  h.nelec_test_mu[1] =  new TH1F("nelec_test_mu2","NElec: NMuon = 1",10,0,10);
  h.nelec_test_mu[2] =  new TH1F("nelec_test_mu3","NElec: NMuon = 2",10,0,10);

  h.test_IVM[0] =  new TH1F("b1 + j1", "b1 j1 IVM  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[1] =  new TH1F("j1+j2+MET", "j1+j2+MET IVM  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[2] =  new TH1F("b1+b2+MET", "b1 b2 MET IVM  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[3] =  new TH1F("b1+j2+MET", "b1 j1 MET IVM  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[4] =  new TH1F("M_all_mu", "M_all (b+j+mu+mu)  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[5] =  new TH1F("M_all_e", "M_all (b+j+e+e)  (70<Mmm<110)", 200, 0, 200);
  h.test_IVM[6] =  new TH1F("b1j1b2j2", "b1j1b2j2 IVM  (70<Mmm<110)", 200, 0, 200);

  h.trijetIVM[0] =  new TH1F("M_jjj_1", "M_bjj (66 < Mll < 116)", 300, 0, 300);
  h.trijetIVM[1] =  new TH1F("M_jjj_2", "M_jjj (66 < Mll < 116)", 300, 0, 300);
  h.trijetIVM[2] =  new TH1F("M_jjj_3", "M_bjj (65 < Mjj < 95) (66 < Mll < 116)", 300, 0, 300);
  
  //=======================================================================================================================
  //                                           SAMLPLE B
  //=======================================================================================================================
 
  h.test_IVM_B = new TH1F("b1 + j1", "b1 j1 IVM  (70<Mmm<110)", 200, 0, 200);
  h.nelec_test[0] =  new TH1F("nelec_test1","NELectron NMuon = 0 (60 < Mjj < 120)",10,0,10);
  h.nelec_test[1] =  new TH1F("nelec_test2","NELectron NMuon = 1 (60 < Mjj < 120)",10,0,10);
  h.nelec_test[2] =  new TH1F("nelec_test3","NELectron NMuon = 2 (60 < Mjj < 120)",10,0,10);
  
  h.nphoton_test[0] = new TH1F("nphoton_test1","NPhoton (60 < Mjj < 120)",10,0,10);
  h.nphoton_test[1] =  new TH1F("nphoton_test2","NPhoton NMuon = 0",10,0,10);
  h.nphoton_test[2] =  new TH1F("nphoton_test3","NPhoton NElectron = 0",10,0,10);
  h.nphoton_test[3] =  new TH1F("nphoton_test4","NPhoton NMuon = 1",10,0,10);
  h.nphoton_test[4] =  new TH1F("nphoton_test5","NPhoton NElectron = 1",10,0,10);
  h.nphoton_test[5] =  new TH1F("nphoton_test6","NPhoton NJet = 1",10,0,10);
  h.nphoton_test[6] =  new TH1F("nphoton_test7","NPhoton NJet = 2",10,0,10);
  
  h.njet_test[0] =  new TH1F("njet_test1","NJet NMuon=0",10,0,10);
  h.njet_test[1] =  new TH1F("njet_test2","NJet NMuon=1",10,0,10);
  h.njet_test[2] =  new TH1F("njet_test4","NJet NElectron=0",10,0,10);
  h.njet_test[3] =  new TH1F("njet_test5","NJet NElectron=1",10,0,10);

  h.MuTransMass[0] = new TH1F("mu_trans_mass1", "Muon Transverse Mass (60 < Mjj < 100)", 200, 0, 200);
  h.eTransMass[0] = new TH1F("e_trans_mass1", "Electron Transverse Mass when NMuon > 0 (60 < Mjj < 100)", 200, 0, 200);
  h.MuTransMass[1] = new TH1F("mu_trans_mass2", "Muon Transverse Mass when NElec > 0 (60 < Mjj < 100)", 200, 0, 200);
  h.eTransMass[1] = new TH1F("e_trans_mass2", "Electron Transverse Mass (60 < Mjj < 100)", 200, 0, 200);
  h.MuTransMass[2] = new TH1F("mu_trans_mass3", "Muon Transverse Mass", 200, 0, 200);
  h.eTransMass[2] = new TH1F("e_trans_mass3", "Electron Transverse Mass", 200, 0, 200);

  h.mjj_test[0] = new TH1F ("m_jj_test1", "Mjj (60 < e_mT < 100)", 200, 0, 200);
  h.mjj_test[1] = new TH1F ("m_jj_test2", "Mjj (60 < mu_mT < 100)", 200, 0, 200);

  //======================================================================================================================
  //                                      SAMPLE A
  //======================================================================================================================

  h.nelec_test[0] =  new TH1F("nelec_test1","NELectron (84 < Mee < 105)",10,0,10);
  h.nphoton_test[0] = new TH1F("nphoton_test1","NPhoton (84 < Mee < 105)",10,0,10);
  h.nmuon_test[0] =  new TH1F("nmuon_test1","NMuon (84 < Mll < 105)",10,0,10);
  h.njet_test[0] = new TH1F("njet_test1","NJet (84 < Mll < 105)",10,0,10);

  h.nelec_test[1] =  new TH1F("nelec_test2","NELectron (50 < Mee < 85)",10,0,10);
  h.nphoton_test[1] = new TH1F("nphoton_test2","NPhoton (50 < Mee < 85)",10,0,10);
  h.nmuon_test[1] =  new TH1F("nmuon_test2","NMuon (50 < Mee < 85)",10,0,10);
  h.njet_test[1] = new TH1F("njet_test2","NJet (50 < Mee < 85)",10,0,10);

  h.nelec_test[2] =  new TH1F("nelec_test3","NELectron (84 < Mmm < 105)",10,0,10);
  h.nphoton_test[2] = new TH1F("nphoton_test3","NPhoton (84 < Mmm < 105)",10,0,10);
  h.nmuon_test[2] =  new TH1F("nmuon_test3","NMuon (84 < Mmm < 105)",10,0,10);
  h.njet_test[2] = new TH1F("njet_test3","NJet (84 < Mmm < 105)",10,0,10);

  h.nelec_test[3] =  new TH1F("nelec_test4","NELectron (50 < Mmm < 85)",10,0,10);
  h.nphoton_test[3] = new TH1F("nphoton_test4","NPhoton (50 < Mmm < 85)",10,0,10);
  h.nmuon_test[3] =  new TH1F("nmuon_test4","NMuon (50 < Mmm < 85)",10,0,10);
  h.njet_test[3] = new TH1F("njet_test4","NJet (50 < Mmm < 85)",10,0,10);

  h.edPhi_test[0] = new TH1F("edPhi_test1", "e1 e2 dPhi (84 < Mee < 105)", 120, -3, 3);
  h.edR_test[0] = new TH1F("edR_test1", "e1 e2 dR (84 < Mee < 105)", 120, 0, 10);
  h.edPhi_test[1] = new TH1F("edPhi_test2", "e1 e2 dPhi (50 < Mee < 85)", 120, -3, 3);
  h.edR_test[1] = new TH1F("edR_test2", "e1 e2 dR (50 < Mee < 85)", 120, 0, 10);

  h.mudPhi_test[0] = new TH1F("mudPhi_test1", "mu1 mu2 dPhi  (84 < Mee < 105)", 120, -3, 3);
  h.mudR_test[0] = new TH1F("mudR_test1", "mu1 mu2 dR  (84 < Mee < 105)", 120, 0, 10);
  h.mudPhi_test[1] = new TH1F("mudPhi_test2", "mu1 mu2 dPhi  (50 < Mee < 85)", 120, -3, 3);
  h.mudR_test[1] = new TH1F("mudR_test2", "mu1 mu2 dR  (50 < Mee < 85)", 120, 0, 10);

  h.nphoton_test[4] = new TH1F("nphoton_test5","NPhoton NJet=0 (84 < Mll < 105)",10,0,10);
  h.nphoton_test[5] = new TH1F("nphoton_test6","NPhoton NJet=1 (84 < Mll < 105)",10,0,10);
  h.nphoton_test[6] = new TH1F("nphoton_test7","NPhoton NJet=2 (84 < Mll < 105)",10,0,10);

  h.m_eey[0] =  new TH1F("M_eey", "M_eey", 200, 0, 200);
  h.m_eey[1] =  new TH1F("M_eey1", "M_eey (50 < Mee < 85)", 200, 0, 200);
  h.m_eey[2] =  new TH1F("M_eey2", "M_eey (84 < Mee < 105)", 200, 0, 200);

  h.m_mumuy[0] =  new TH1F("M_mumuy", "M_mumuy", 200, 0, 200);
  h.m_mumuy[1] =  new TH1F("M_mumuy1", "M_mumuy (50 < Mmumu < 85)", 200, 0, 200);
  h.m_mumuy[2] =  new TH1F("M_mumuy2", "M_mumuy (84 < Mmumu < 105)", 200, 0, 200);

  h.jetInvM_test[0] = new TH1F("M_jet_test1", "M_jet (84 < Mll < 105)", 200, 0, 200);
  h.jetInvM_test[1] = new TH1F("M_jet_test2", "M_jet (50 < Mll < 85)", 200, 0, 200);
}


