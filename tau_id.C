#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include "TStyle.h"
#include <vector>
//#include <cmath>

using namespace ROOT::VecOps;
using rvec_f = RVec<float>;
using rvec_b = RVec<bool>;
using rvec_i = RVec<int>;
using rvec_u = RVec<unsigned char>;

template <typename T>
void plot(T hist, TString name){
  TCanvas * c = new TCanvas("c",Form("c_%s", name.Data()));
  gStyle->SetOptStat(0);
  hist->Write();
  hist->DrawClone();
  c->Print(Form("%s.pdf",name.Data()));
}

template <typename T>
void plot_fit(T hist, TString name){
  TCanvas * c = new TCanvas("c",Form("c_%s", name.Data()));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatFontSize(0.03);
  hist->Write();
  hist->DrawClone();
  c->Print(Form("%s.pdf",name.Data()));
}

rvec_i good_idx(rvec_i g){
  vector<int> out;
  for(int i = 0; i < g.size(); i++){
    if( g[i] ) out.push_back( i );
  }
  return out; 
}

rvec_i identify_tau_jet(rvec_i j_idx, rvec_i tau_idx){
  vector<int> matched_idx;
  for(int i = 0; i < j_idx.size(); i++){
    for(int j = 0; j < tau_idx.size(); j++){
      if( j_idx[i] != tau_idx[j]) continue;
      matched_idx.push_back(j_idx[i]);
    }
  }
  return matched_idx;
}

int size_of_vec(rvec_i a){
  int n = a.size();
  return n;
}

int goodmuon_idx(rvec_i good){
  int idx = -1;
  for(int i = 0; i < good.size(); i++){
    if(good[i]) idx = i;
    break;
  }
  return idx;
}

int goodtaujet_idx(rvec_i good, int goodmu, rvec_f mu_eta, rvec_f mu_phi, rvec_f jet_eta, rvec_f jet_phi){
  int idx = -1;
  for(int i = 0; i < good.size(); i++){
    float dR = DeltaR(mu_eta[goodmu],jet_eta[i],mu_phi[goodmu],jet_phi[i]);
    if(good[i] && dR > 0.5){
      idx = i;
      break;
    }
  }
  return idx;
}

float calculate_mt(int goodmu, rvec_f mu_pt, float met_pt, rvec_f mu_phi, float met_phi){
  float muon_pt = mu_pt[goodmu];
  float muon_phi = mu_phi[goodmu];
  float dPhi = DeltaPhi(muon_phi, met_phi);
  float mt = sqrt(2*muon_pt*met_pt*(1-cos(dPhi)));
  return mt;
}

float calculate_dz(int goodmu, rvec_f mu_pt, rvec_f mu_eta, rvec_f mu_phi, rvec_f mu_mass, int good_j, rvec_f j_pt, rvec_f j_eta, rvec_f j_phi, rvec_f j_mass, float met_pt, float met_phi){
  int i = goodmu;
  int j = good_j;
  TLorentzVector a, b;
  a.SetPtEtaPhiM(mu_pt[i],mu_eta[i],mu_phi[i],mu_mass[i]);
  b.SetPtEtaPhiM(j_pt[j],j_eta[j],j_phi[j],j_mass[j]);
  
  float z = min(a.Phi(), b.Phi());
  float dphi = a.DeltaPhi(b);
  float zeta = z + 0.5*dphi;
  float diff_vis = zeta - (a+b).Phi();
  float diff_met = zeta - met_phi;
  float pt_vis_zeta = (a+b).Pt()*cos(diff_vis);
  float pt_met_zeta = met_pt*cos(diff_met);
  float D_zeta = pt_met_zeta - 0.85*pt_vis_zeta;

  return D_zeta;
}

int n_vetomuons(int goodmu_idx, rvec_i vetomuon){
  int n = 0;
  for(int i = 0; i < vetomuon.size(); i++){
    if(vetomuon[i]){
      if(i!=goodmu_idx) n++;
    }
  }
  return n;
}

float calculate_mass(int good_a, rvec_f a_pt, rvec_f a_eta, rvec_f a_phi, rvec_f a_mass, int good_b, rvec_f b_pt, rvec_f b_eta, rvec_f b_phi, rvec_f b_mass){
  float inv_mass;
  int i = good_a;
  int j = good_b;
  TLorentzVector a, b;
  a.SetPtEtaPhiM(a_pt[i],a_eta[i],a_phi[i],a_mass[i]); 
  b.SetPtEtaPhiM(b_pt[j],b_eta[j],b_phi[j],b_mass[j]);
  inv_mass = (a+b).M();
  
  return inv_mass;
}

int z_mass_idx(rvec_f z_mass, rvec_i good_idx){
  int idx = -1;
  float tmp = 9999;
  for(int i = 0; i < good_idx.size(); i++){
    if( abs(z_mass[i] - 91.2) < tmp){
      tmp = z_mass[i] - 91.2;
      idx = i;
    }
  }
  return idx;
}

int tag_n_probe(int taucand_idx, rvec_i goodtau, rvec_i tau_jetid, rvec_u tau_vs_jetid){
  int matched = 0;
  for(int i = 0; i < goodtau.size(); i++){
    if(goodtau[i] && (tau_vs_jetid[i] & 8) && (taucand_idx == tau_jetid[i])){
      matched = 1;
    }
  }
  return matched;
}

void tau_id(){
  //ROOT::RDataFrame df("Events", "/xrootd/store/data/Run2016B_ver2/SingleMuon/NANOAOD/Nano25Oct2019_ver2-v1/20000/0289BB82-DF74-A74F-A9F7-AE57B5041EE9.root");
  ROOT::RDataFrame df("Events", "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/034659D3-DA7F-9C43-B604-FADBF92587EC.root");
//  ROOT::RDataFrame df("Events", {"/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/034659D3-DA7F-9C43-B604-FADBF92587EC.root",
//          "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/10F251FC-3AFA-E148-84F2-8908B099A9AB.root",
//          "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/15DB1CA0-7691-6B48-AA5B-B49EEB7C2AD0.root",
//          "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/1840125F-6253-F142-94ED-DD293C3387D3.root",
//          "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/1AF74F4E-72DF-804E-B7D5-01EE53053AAE.root",
//          "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/21245745-A7A3-6144-99C5-D4255B38E428.root"}); 
//  ROOT::RDataFrame df("Events", {"/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/034659D3-DA7F-9C43-B604-FADBF92587EC.root",
//                                 "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/10F251FC-3AFA-E148-84F2-8908B099A9AB.root",
//                                 "/xrootd/store/data/Run2018D/SingleMuon/NANOAOD/Nano25Oct2019-v1/100000/15DB1CA0-7691-6B48-AA5B-B49EEB7C2AD0.root"});
  //ROOT::RDataFrame df("Events",input);
  // HLT Trigger
  auto df_S1_HLT = df.Filter("HLT_IsoMu24", "HLT Trigger IsoMu24");
  
  // Muon Selection and Veto
  auto df_S1_Muon = df_S1_HLT.Define("goodmuon","Muon_pt > 25 && abs(Muon_eta) < 2.4 && Muon_mediumId && Muon_pfRelIso04_all < 0.15 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045")
                             .Define("goodmuon_idx",goodmuon_idx,{"goodmuon"})
                             .Define("vetomuon","Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_mediumId && Muon_pfRelIso04_all < 0.15 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045")
                             .Filter("Sum(goodmuon) >= 1","Muon Selection");
  
  auto df_S1_Muon_veto = df_S1_Muon.Define("n_vetomuons",n_vetomuons,{"goodmuon_idx","vetomuon"})
                                   .Filter("n_vetomuons == 0","Muon veto");
  
  auto df_S1_Electron_veto = df_S1_Muon_veto.Define("vetoelectron","Electron_pt > 10 && abs(Electron_eta) < 2.5 && abs(Electron_dxy) < 0.045 && abs(Electron_dz) < 0.2 && Electron_pfRelIso03_all < 0.10 && Electron_mvaFall17V2noIso_WP90")
                                            .Filter("Sum(vetoelectron) == 0","Electron veto");
  
  // Tau candidate selection
  auto df_S1_Tau_cand = df_S1_Electron_veto.Define("goodtaujet","Jet_pt > 20 && abs(Jet_eta) < 2.3")
                                           .Define("goodtaujet_idx",goodtaujet_idx,{"goodtaujet","goodmuon_idx","Muon_eta","Muon_phi","Jet_eta","Jet_phi"})
                                           .Filter("goodtaujet_idx >= 0","Tau Candidate Selection");

 
  // Background Suppression
  auto df_S1_bkg = df_S1_Tau_cand.Define("mt",calculate_mt,{"goodmuon_idx","Muon_pt","MET_pt","Muon_phi","MET_phi"})
                                 .Define("Dz",calculate_dz,{"goodmuon_idx","Muon_pt","Muon_eta","Muon_phi","Muon_mass","goodtaujet_idx","Jet_pt","Jet_eta","Jet_phi","Jet_mass","MET_pt","MET_phi"})
                                 .Define("abs_dEta","abs(Muon_eta[goodmuon_idx] - Jet_eta[goodtaujet_idx])");

  // Z mass reconstruction
  auto df_S1_Z = df_S1_bkg.Define("mu_jet_mass",calculate_mass,{"goodmuon_idx", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "goodtaujet_idx", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"});

  auto df_S1_abs_dEta_sel = df_S1_Z.Filter("abs_dEta < 1.5","Delta eta cut"); 
  auto df_S1_mt_sel = df_S1_abs_dEta_sel.Filter("mt < 60","mt selection");
  auto df_S1_Dz_sel = df_S1_mt_sel.Filter("Dz > -25","Dz Selection"); 

  // Tau selection
  auto df_S1_Tau = df_S1_Dz_sel.Define("goodtaus","Tau_pt > 20 && abs(Tau_eta) < 2.3 && Tau_idDecayModeNewDMs && (Tau_idDeepTau2017v2p1VSe & 2) && (Tau_idDeepTau2017v2p1VSmu & 8)")
                               .Define("one_prong","Tau_decayMode == 0 || Tau_decayMode == 1")
                               .Define("three_prong","Tau_decayMode == 10 || Tau_decayMode == 11")
                               .Define("goodtaus_tot","goodtaus && (one_prong || three_prong)")
                               .Filter("Sum(goodtaus_tot) >= 1","Tau Selection")
                               .Define("probed",tag_n_probe,{"goodtaujet_idx","goodtaus_tot","Tau_jetIdx","Tau_idDeepTau2017v2p1VSe"});

  auto df_S1_probe_Tau = df_S1_Tau.Filter("probed==1","loose WP");
  
  // Histograms
  auto h_mu_jet_mass = df_S1_Z.Histo1D({"h_mu_jet_mass","h_mu_jet_mass",20,50,150},"mu_jet_mass");
  auto h_tau_cand = df_S1_Z.Define("tau_cand_pt","Jet_pt[goodtaujet_idx]").Histo1D({"h_tau_cand_pt","h_tau_cand_pt",20,0,100},"tau_cand_pt");
  auto h_abs_dEta = df_S1_bkg.Histo1D({"h_abs_dEta","h_abs_dEta",25,0,5},"abs_dEta"); 
  auto h_mt = df_S1_bkg.Histo1D({"h_mt","h_mt",15,0,150},"mt");
  auto h_Dz = df_S1_bkg.Histo1D({"h_Dz","h_Dz",40,-250,150},"Dz");
  auto h_z_mass = df_S1_probe_Tau.Histo1D({"h_z_mass_l","h_z_mass_l",70,40,110},"mu_jet_mass");

  // Asymmetric Gaussian Fitting
  TF1 *f_fit = new TF1("f", "[0] / abs((x < [1]) ? [2] : [3]) * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]))", 40., 120.);
  f_fit->SetParameters(100.,70.,-10.,10.);
  f_fit->SetParNames("Constant","Mean","Sigma1","Sigma2");
  h_z_mass->Fit(f_fit);

  TF1 *f_std = new TF1("f_std","exp(-x^2/2)/sqrt(2*pi)",40.,120.);



  int total_evt = df_S1_Tau.Count().GetValue();
  cout<<"Total Event = "<<total_evt<<endl;
  
  // Saving histograms
  TFile f("output.root", "recreate");
  plot( h_mu_jet_mass, "h_mu_jet_mass");
  plot( h_tau_cand, "h_tau_cand");
  plot( h_abs_dEta, "h_abs_dEta");
  plot( h_mt, "h_mt");
  plot( h_Dz, "h_Dz");
  plot_fit( h_z_mass, "h_z_mass_l");
  
  f.Close();
  
//  auto report = df_S1_Tau.Report();
//  report->Print();
}
