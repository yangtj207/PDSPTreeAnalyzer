#include "PionQE.h"
#include "TLorentzVector.h"

void PionQE::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionQE::executeEvent(){

  pi_type = GetPi2ParType();
  pi_type_str = Form("%d", pi_type);
  FillHist("beam_cut_flow", 0.5, 1., 20, 0., 20.);

  pi_truetype = GetPiTrueType();
  CalTrueQEVars();
  if (pi_truetype == pitrue::kQE) FillQEPlots("QE_All");
  
  // -- 1. Beam instruments
  if(!PassBeamScraperCut()) return;
  FillHist("beam_cut_flow", 1.5, 1., 20, 0., 20.);
  if(!PassBeamMomentumWindowCut()) return;
  FillHist("beam_cut_flow", 2.5, 1., 20, 0., 20.);
  
  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  KE_beam_inst = map_BB[211] -> MomentumtoKE(P_beam_inst);
  exp_trk_len_beam_inst = map_BB[211] -> RangeFromKESpline(KE_beam_inst);
  trk_len_ratio = evt.reco_beam_alt_len / exp_trk_len_beam_inst;

  mass_beam = 139.57;
  //P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  //KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_ff_reco = KE_beam_inst; // FIXME : study energy measurement bias due to energy loss and spectrometer
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;
  double P_reweight = 1.;
  if(!IsData) P_reweight = MCCorr -> MomentumReweight_SF("Pion_Calosize_0p5", P_beam_inst, 0.);

  if(!Pass_Beam_PID(211)) return;
  FillHist("beam_cut_flow", 3.5, 1., 20, 0., 20.);
  
  // -- 2. TPC info 
  if(!PassPandoraSliceCut()) return;
  FillHist("beam_cut_flow", 4.5, 1., 20, 0., 20.);

  if(!PassCaloSizeCut()) return;
  FillHist("beam_cut_flow", 5.5, 1., 20, 0., 20.);
  FillBeamPlots("Beam_CaloSize", P_reweight);
  if (pi_type == pi2::kPiQE) FillQEPlots("QE_CaloSize");

  if(!PassAPA3Cut()) return;
  FillHist("beam_cut_flow", 6.5, 1., 20, 0., 20.);
  FillBeamPlots("Beam_APA3", P_reweight);
  if (pi_type == pi2::kPiQE) FillQEPlots("QE_APA3");

  if(!PassMichelScoreCut()) return;
  FillHist("beam_cut_flow", 7.5, 1., 20, 0., 20.);
  FillBeamPlots("Beam_MichelScore", P_reweight);
  if (pi_type == pi2::kPiQE) FillQEPlots("QE_MichelScore");

  if(KE_end_reco < 100.) return;
  FillHist("beam_cut_flow", 8.5, 1., 20, 0., 20.);
  FillBeamPlots("Beam_KE_end", P_reweight);
  if (pi_type == pi2::kPiQE) FillQEPlots("QE_KE_end");

  // == Functions to study beam
  //Run_beam_dEdx_vector();
 
  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> loose_pions = SelectLoosePions(daughters_all);
  Run_Daughter("AllRecoDaughters", daughters_all);
  Run_Daughter("LoosePions", loose_pions);

  if (loose_pions.size() == 1){ //QE candidates
    FillHist("beam_cut_flow", 9.5, 1., 20, 0., 20.);
    FillBeamPlots("Beam_QE", P_reweight);
    if (pi_type == pi2::kPiQE){
      FillQEPlots("QE_QE");
    }
    CalRecoQEVars(loose_pions[0]);
    if (recoQE_KEPi1>0){
      FillQERecoPlots("QE_QE", P_reweight);
    }
  }
  
  // == pion selection cutflow
  if(!IsData){
    int N_pion_nocut = 0;
    int N_pion_trkscore = 0;
    int N_pion_chi2proton = 0;
    int N_pion_trklen_upper = 0;
    int N_pion_trklen_lower = 0;

    int N_proton_nocut = 0;
    int N_proton_trkscore = 0;
    int N_proton_chi2proton = 0;
    int N_proton_trklen_upper = 0;
    int N_proton_trklen_lower = 0;
    
    vector<Daughter> pion_trkscore = SelectPions_trkscore(daughters_all, 0.5);
    vector<Daughter> pion_chi2proton = SelectPions_chi2(pion_trkscore, 60.);
    vector<Daughter> pion_trklen_upper = SelectPions_trklen_upper(pion_chi2proton, 180.);
    vector<Daughter> pion_trklen_lower = SelectPions_trklen_lower(pion_trklen_upper, 10.);
    for(unsigned int i = 0; i < daughters_all.size(); i++){
      if(daughters_all.at(i).PFP_true_byHits_PDG() == 211) N_pion_nocut++;
      if(daughters_all.at(i).PFP_true_byHits_PDG() == 2212) N_proton_nocut++;
    }
    for(unsigned int i = 0; i < pion_trkscore.size(); i++){
      if(pion_trkscore.at(i).PFP_true_byHits_PDG() == 211) N_pion_trkscore++;
      if(pion_trkscore.at(i).PFP_true_byHits_PDG() == 2212) N_proton_trkscore++;
    }
    for(unsigned int i = 0; i < pion_chi2proton.size(); i++){
      if(pion_chi2proton.at(i).PFP_true_byHits_PDG() == 211) N_pion_chi2proton++;
      if(pion_chi2proton.at(i).PFP_true_byHits_PDG() == 2212) N_proton_chi2proton++;
    }
    for(unsigned int i = 0; i < pion_trklen_upper.size(); i++){
      if(pion_trklen_upper.at(i).PFP_true_byHits_PDG() == 211) N_pion_trklen_upper++;
      if(pion_trklen_upper.at(i).PFP_true_byHits_PDG() == 2212) N_proton_trklen_upper++;
    }
    for(unsigned int i = 0; i < pion_trklen_lower.size(); i++){
      if(pion_trklen_lower.at(i).PFP_true_byHits_PDG() == 211) N_pion_trklen_lower++;
      if(pion_trklen_lower.at(i).PFP_true_byHits_PDG() == 2212) N_proton_trklen_lower++;
    }
    
    FillHist("pion_sel_cutflow", 0., N_pion_nocut, 5., -0.5, 4.5);
    FillHist("pion_sel_cutflow", 1., N_pion_trkscore, 5., -0.5, 4.5);
    FillHist("pion_sel_cutflow", 2., N_pion_chi2proton, 5., -0.5, 4.5);
    FillHist("pion_sel_cutflow", 3., N_pion_trklen_upper, 5., -0.5, 4.5);
    FillHist("pion_sel_cutflow", 4., N_pion_trklen_lower, 5., -0.5, 4.5);

    FillHist("proton_sel_cutflow", 0., N_proton_nocut, 5., -0.5, 4.5);
    FillHist("proton_sel_cutflow", 1., N_proton_trkscore, 5., -0.5, 4.5);
    FillHist("proton_sel_cutflow", 2., N_proton_chi2proton, 5., -0.5, 4.5);
    FillHist("proton_sel_cutflow", 3., N_proton_trklen_upper, 5., -0.5, 4.5);
    FillHist("proton_sel_cutflow", 4., N_proton_trklen_lower, 5., -0.5, 4.5);
  }
}

void PionQE::FillBeamPlots(TString beam_selec_str, double weight){

  // == Fit results after calo-size cut
  /*
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData){
    Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;
    Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;
    Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;
  }
  */
     
  // == Comparison between beam spectrometer track and TPC reco track
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX", evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY", evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ", evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ", evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst", P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh", beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton", chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX_" + pi_type_str, evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY_" + pi_type_str, evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_" + pi_type_str, evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_" + pi_type_str, delta_X_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_" + pi_type_str, delta_Y_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_cos_delta_spec_TPC_" + pi_type_str, cos_delta_spec_TPC, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh_" + pi_type_str, beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_theta_" + pi_type_str, beam_TPC_theta, weight, 5000., -1., 4.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_phi_" + pi_type_str, beam_TPC_phi, weight, 8000., -4., 4.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ_" + pi_type_str, evt.reco_beam_calo_endZ, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_alt_len_" + pi_type_str, evt.reco_beam_alt_len, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst_" + pi_type_str, P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton_" + pi_type_str, chi2_proton, weight, 10000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_trk_len_ratio_" + pi_type_str, trk_len_ratio, weight, 1000., 0., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_ff_" + pi_type_str, KE_ff_reco, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_end_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);

  /*
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma", Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma_" + pi_type_str, Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_X_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_Y_over_sigma, weight, 2000., -10., 10.);
  */
}

void PionQE::FillQEPlots(TString beam_selec_str){

  JSFillHist(beam_selec_str, beam_selec_str + "_QE_Q2", QE_Q2, 1.0, 20, 0., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_KEPi0", QE_KEPi0, 1.0, 20, 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_KEPi1", QE_KEPi1, 1.0, 20, 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_AngPi", QE_AngPi, 1.0, 20, 0., 180.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_nu", QE_nu, 1.0, 20, 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_EQE", QE_EQE, 1.0, 20, 0., 700.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QE_nuQ2", QE_Q2, QE_nu, 1.0, 100, 0., 1., 100, 0., 500.);
}

void PionQE::FillQERecoPlots(TString beam_selec_str, double weight){

  if (pi_type == pi2::kPiQE){
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_Q2", QE_Q2, recoQE_Q2, weight, 20, 0., 0.6,20,0.,0.6);
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_KEPi0", QE_KEPi0, recoQE_KEPi0, weight, 20, 0., 500., 20, 0., 500.);
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_KEPi1", QE_KEPi1, recoQE_KEPi1, weight, 20, 0., 500., 20, 0., 500.);
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_AngPi", QE_AngPi, recoQE_AngPi, weight, 20, 0., 180., 20, 0., 180.);
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_nu", QE_nu, recoQE_nu, weight, 20, 0., 500., 20, 0., 500.);
    JSFillHist(beam_selec_str, beam_selec_str + "_QE2D_EQE", QE_EQE, recoQE_EQE, weight, 20, 0., 700., 20, 0., 700.);
  }
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_Q2_" + pi_type_str, recoQE_Q2, weight, 20., 0., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_KEPi0_" + pi_type_str, recoQE_KEPi0, weight, 20., 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_KEPi1_" + pi_type_str, recoQE_KEPi1, weight, 20., 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_AngPi_" + pi_type_str, recoQE_AngPi, weight, 20., 0., 180.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_nu_" + pi_type_str, recoQE_nu, weight, 20., 0., 500.);
  JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_dEQE_" + pi_type_str, recoQE_EQE-recoQE_KEPi0-M_pion, weight, 20., -500., 500.);
  if (recoQE_Q2>0.05 && recoQE_Q2<0.15){
    JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_nu1_" + pi_type_str, recoQE_nu, weight, 20., 0., 500.);
  }
  if (recoQE_Q2>0.15 && recoQE_Q2<0.25){
    JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_nu2_" + pi_type_str, recoQE_nu, weight, 20., 0., 500.);
  }
  if (recoQE_Q2>0.25 && recoQE_Q2<0.35){
    JSFillHist(beam_selec_str, beam_selec_str + "_QEreco_nu3_" + pi_type_str, recoQE_nu, weight, 20., 0., 500.);
  }
}

void PionQE::Run_Daughter(TString daughter_sec_str, const vector<Daughter>& daughters){

  JSFillHist(daughter_sec_str, "N_daughters", daughters.size(), 1., 10., -0.5, 9.5);
  int N_true_proton_PID = 0;
  int N_true_piplus_PID	= 0;
  int N_true_piminus_PID = 0;
  for(unsigned int i = 0; i < daughters.size(); i++){
    Daughter this_daughter = daughters.at(i);
    double this_trklen = this_daughter.allTrack_alt_len();
    double this_trkscore = this_daughter.PFP_trackScore();
    double this_chi2_pion = Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), 211);
    double this_chi2_proton = Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), 2212);
    double this_michelscore = 0;
    if (this_daughter.allTrack_vertex_nHits()) this_michelscore = this_daughter.allTrack_vertex_michel_score()/this_daughter.allTrack_vertex_nHits();
    JSFillHist(daughter_sec_str, "daughters_chi2_pion", this_chi2_pion, 1., 1000., 0., 1000.);
    JSFillHist(daughter_sec_str, "daughters_chi2_proton", this_chi2_proton, 1., 1000., 0., 1000.);
    JSFillHist(daughter_sec_str, "daughters_trkscore", this_trkscore, 1., 1000., 0., 1.);
    JSFillHist(daughter_sec_str, "daughters_trklen", this_trklen, 1., 200., 0., 200.);
    JSFillHist(daughter_sec_str, "daughters_michelscore", this_michelscore, 1., 1000., 0., 1.);
    
    if(!IsData){
      int this_true_PDG = this_daughter.PFP_true_byHits_PDG();
      JSFillHist(daughter_sec_str, Form("daughters_chi2_pion_truePDG%d", this_true_PDG), this_chi2_pion, 1., 1000., 0., 1000.);
      JSFillHist(daughter_sec_str, Form("daughters_chi2_proton_truePDG%d", this_true_PDG), this_chi2_proton, 1., 1000., 0., 1000.);
      JSFillHist(daughter_sec_str, Form("daughters_trkscore_truePDG%d", this_true_PDG), this_trkscore, 1., 1000., 0., 1.);
      JSFillHist(daughter_sec_str, Form("daughters_trklen_truePDG%d", this_true_PDG), this_trklen, 1., 200., 0., 200.);
      JSFillHist(daughter_sec_str, Form("daughters_michelscore_truePDG%d", this_true_PDG), this_michelscore, 1., 1000., 0., 1.);
      
      if(this_daughter.PFP_true_byHits_PDG() == 2212) N_true_proton_PID++;
      else if(this_daughter.PFP_true_byHits_PDG() == 211) N_true_piplus_PID++;
      else if(this_daughter.PFP_true_byHits_PDG() == -211) N_true_piminus_PID++;
    }
  }

  JSFillHist(daughter_sec_str, "N_true_proton_PID", N_true_proton_PID, 1., 10., -0.5, 9.5);
  JSFillHist(daughter_sec_str, "N_true_piplus_PID", N_true_piplus_PID, 1., 10., -0.5, 9.5);
  JSFillHist(daughter_sec_str, "N_true_piminus_PID", N_true_piminus_PID, 1., 10., -0.5, 9.5);
  JSFillHist(daughter_sec_str, "N_true_piplus_PID_vs_N_true_proton_PID", N_true_piplus_PID, N_true_proton_PID, 1., 5., -0.5, 4.5, 5., -0.5, 4.5);
}

void PionQE::TrueDaughterStudy(const vector<TrueDaughter>& true_daughters){
  
}

std::vector<Daughter> PionQE::SelectLoosePions(const vector<Daughter>& in){

  vector<Daughter> out;
  double cut_trackScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_trk_len_upper = 180.;
  double cut_trk_len_lower = 10.;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_chi2 > cut_chi2_proton && this_in.allTrack_alt_len() < cut_trk_len_upper && this_in.allTrack_alt_len() > cut_trk_len_lower) out.push_back(this_in);
  }
  return out;
}

std::vector<Daughter> PionQE::SelectPions_trkscore(const vector<Daughter>& in, double cut_trackScore){

  vector<Daughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    if(this_in.PFP_trackScore() > cut_trackScore) out.push_back(this_in);
  }
  return out;
}

std::vector<Daughter> PionQE::SelectPions_chi2(const vector<Daughter>& in, double cut_chi2_proton){

  vector<Daughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_chi2 > cut_chi2_proton) out.push_back(this_in);
  }
  return out;
}

std::vector<Daughter> PionQE::SelectPions_trklen_upper(const vector<Daughter>& in, double cut_trk_len_upper){

  vector<Daughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.allTrack_alt_len() < cut_trk_len_upper) out.push_back(this_in);
  }
  return out;
}

std::vector<Daughter> PionQE::SelectPions_trklen_lower(const vector<Daughter>& in, double cut_trk_len_lower){

  vector<Daughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.allTrack_alt_len() > cut_trk_len_lower) out.push_back(this_in);
  }
  return out;
}

void PionQE::CalTrueQEVars(){

  if (pi_truetype != pitrue::kQE) return;

  double Eb = 4;
  TLorentzVector beam, pion, proton;

  beam.SetPxPyPzE(evt.true_beam_endPx*1000,
                  evt.true_beam_endPy*1000,
                  evt.true_beam_endPz*1000,
                  sqrt(pow(evt.true_beam_endP*1000,2)+M_pion*M_pion));

  for (size_t i = 0; i<evt.true_beam_daughter_PDG->size(); ++i){
    int pdg = evt.true_beam_daughter_PDG->at(i);
    double px = evt.true_beam_daughter_startPx->at(i)*1000;
    double py = evt.true_beam_daughter_startPy->at(i)*1000;
    double pz = evt.true_beam_daughter_startPz->at(i)*1000;
    double p = sqrt(px*px+py*py+pz*pz);

    if (abs(pdg) == 211){
      double E = sqrt(p*p+M_pion*M_pion);
      if (E>pion.E()){
        pion.SetPxPyPzE(px,py,pz,E);
      }
    }
    if (pdg == 2212){
      double E = sqrt(p*p+M_proton*M_proton);
      if (E>proton.E()){
        proton.SetPxPyPzE(px,py,pz,E);
      }
    }
  }

  if (beam.E() && pion.E()){
    QE_Q2 = -(beam - pion).Mag2()*1e-6; //GeV^2
    QE_KEPi0 = beam.E() - M_pion;
    QE_KEPi1 = pion.E() - M_pion;
    QE_AngPi = pion.Vect().Angle(beam.Vect())*180/TMath::Pi();
    QE_nu = beam.E() - pion.E();
    QE_EQE = (pow(M_proton,2)-pow(M_proton-Eb,2)-pow(M_pion,2)+2*(M_proton-Eb)*pion.E())/2/(M_proton-Eb-pion.E()+pion.Vect()*beam.Vect()/beam.Vect().Mag());
  }
}

void PionQE::CalRecoQEVars(Daughter pion){

  double KE_hypfit_like = KE_Hypfit_Likelihood(pion, 211);
  if (KE_hypfit_like<0) return;
  double E = KE_hypfit_like+M_pion;
  double p = sqrt(pow(E,2)-pow(M_pion,2));
  TVector3 p_vec_pion(pion.allTrack_endX() - pion.allTrack_startX(), pion.allTrack_endY() - pion.allTrack_startY(), pion.allTrack_endZ() - pion.allTrack_startZ()); 
  p_vec_pion = p_vec_pion.Unit();
  double px = p*p_vec_pion.X();
  double py = p*p_vec_pion.Y();
  double pz = p*p_vec_pion.Z();
  TLorentzVector Pion;
  Pion.SetPxPyPzE(px, py, pz, E);
  double beam_p = sqrt(pow(E_end_reco,2) - pow(M_pion,2));
  double beam_px = beam_p * evt.reco_beam_trackDirX;
  double beam_py = beam_p * evt.reco_beam_trackDirY;
  double beam_pz = beam_p * evt.reco_beam_trackDirZ;
  TLorentzVector Beam(beam_px, beam_py, beam_pz, E_end_reco);
  recoQE_Q2 = -(Beam - Pion).Mag2()*1e-6; //GeV^2
  recoQE_KEPi0 = Beam.E() - M_pion;
  recoQE_KEPi1 = Pion.E() - M_pion;
  recoQE_AngPi = Pion.Vect().Angle(Beam.Vect())*180/TMath::Pi();
  recoQE_nu = Beam.E() - Pion.E();
  double Eb = 4;
  recoQE_EQE = (pow(M_proton,2)-pow(M_proton-Eb,2)-pow(M_pion,2)+2*(M_proton-Eb)*Pion.E())/2/(M_proton-Eb-Pion.E()+Pion.Vect()*Beam.Vect()/Beam.Vect().Mag());
}

PionQE::PionQE(){

}

PionQE::~PionQE(){

}
