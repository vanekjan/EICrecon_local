/*
  TOF matching efficiency
  Macro for production of tof_eff_Dmp_run16.root input for the data-driven fast simulator
*/

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TGraph.h"
#include"TF1.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"
#include"TLine.h"
#include"TChain.h"
#include"TLorentzVector.h"

using namespace std;

void Plot_histos_pi_rejection(int e_energy = 18, int p_energy = 275)
{
  if( !(e_energy == 18 && p_energy == 275) && !(e_energy == 10 && p_energy == 100) && !(e_energy == 5 && p_energy == 41))
  {
    cout<<"Invalid beam energies."<<endl;

    return;
  }
/*
  const int nQ2bins = 2;
  float const Q2_bins[nQ2bins+1] = { 1., 2., 4. };

  const int nyInelParBins = 2;
  float const y_bins[nQ2bins+1] = { 0., 0.001, 0.002 };
*/

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //____________________________________________________
  //draw pi eCALion plot

  double array_mom_bins[8] = {0.1, 0.2, 0.5, 1., 2., 5., 10., 20.};
  double array_pi_false_rate_85[8] = {0.6332832672698294, 0.7495818158985306, 0.00930384575910461, 0.001692827694491846, 0.0001898238241173789, 0.00020018016214593134, 0.000536412900269677, 0.0006430092230696459};
  double array_pi_false_rate_95[8] = {0.8838860654090215, 0.9228366502089709, 0.02228665375203912, 0.0036237053948322794, 0.00048096291971113834, 0.0006523112180272588, 0.0022770026949322946, 0.0018829746368912276};


  TGraph *g_pi_false_rate_85 = new TGraph(8, array_mom_bins, array_pi_false_rate_85);

  TGraph *g_pi_false_rate_95 = new TGraph(8, array_mom_bins, array_pi_false_rate_95);

  TCanvas *pi_false_rate_can = new TCanvas("pi_false_rate_can", "pi_false_rate_can", 1200, 1000);

  pi_false_rate_can->cd();

  gPad->SetLogy();

  g_pi_false_rate_85->SetNameTitle("g_pi_false_rate_85", "g_pi_false_rate_85");
  g_pi_false_rate_85->GetXaxis()->SetTitle("p (GeV/#it{c})");
  g_pi_false_rate_85->GetXaxis()->CenterTitle();
  g_pi_false_rate_85->GetYaxis()->SetTitle("Rejection factor");
  g_pi_false_rate_85->GetYaxis()->CenterTitle();
  g_pi_false_rate_85->SetMarkerStyle(20);
  g_pi_false_rate_85->SetMarkerColor(kRed);
  g_pi_false_rate_85->SetLineWidth(2);
  g_pi_false_rate_85->SetLineColor(kRed);
  g_pi_false_rate_85->Draw("l a");

  g_pi_false_rate_95->SetNameTitle("g_pi_false_rate_95", "g_pi_false_rate_95");
  g_pi_false_rate_95->SetMarkerStyle(20);
  g_pi_false_rate_95->SetMarkerColor(kBlue);
  g_pi_false_rate_95->SetLineWidth(2);
  g_pi_false_rate_95->SetLineColor(kBlue);
  g_pi_false_rate_95->Draw("l same");

  TLegend *eCAL_leg = new TLegend(0.5, 0.65, 0.79, 0.75);
  eCAL_leg->AddEntry(g_pi_false_rate_85, "eCAL 85% eff.", "l");
  eCAL_leg->AddEntry(g_pi_false_rate_95, "eCAL 95% eff.", "l");
  eCAL_leg->SetBorderSize(0);
  eCAL_leg->SetFillColorAlpha(0, 0.01);
  eCAL_leg->Draw("same");

  pi_false_rate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/pi_false_rate.png");

  //___________________________________________________


  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 11;
  float const mom_bins[nMomBins+1] = { 0,0.5,1,1.5,2,3,4,5,6,7,10, 18 };

  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -4, -3, -2, -1, 0 };

  //load all files
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output.root", e_energy, p_energy, e_energy, p_energy), "READ");
  //TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i_p_scat_bins-output.root", e_energy, p_energy, e_energy, p_energy), "READ");

  //eta distributions in multiple Q^2 and inelasticity bins
  TH1D *h_eta_scat_ele[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_pi_minus[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_pi_minus_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_pi_minus_eCAL_85[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_pi_minus_eCAL_85_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_K_minus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_K_minus_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_anti_proton[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_anti_proton_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  //p distributions in multiple eta bins
  //use just wide Q^2 and y selection

  //eta bins defined higher + 1 eta bin for pfRICH acceptance window
  TH1D *h_p_scat_ele[nEtaBins+1];

  TH1D *h_p_pi_minus[nEtaBins+1];

  TH1D *h_p_pi_minus_eCAL_85[nEtaBins+1];
  TH1D *h_p_pi_minus_eCAL_95[nEtaBins+1];

  TH1D *h_p_pi_minus_eCAL_85_pfRICH[nEtaBins+1];
  TH1D *h_p_pi_minus_eCAL_95_pfRICH[nEtaBins+1];

  TH1D *h_p_pi_minus_pfRICH[nEtaBins+1];



  //reco tracks histograms
  //with eCAL
  TH1D* h_energy_RC = (TH1D*)inFile->Get("h_energy_RC");

  TH1D* h_momentum_RC = (TH1D*)inFile->Get("h_momentum_RC");

  TH1D *h_E_over_p_RC = (TH1D*)inFile->Get("h_E_over_p_RC");

  TH1D* h_Q2_RC = (TH1D*)inFile->Get("h_Q2_RC");

  TH1D *h_y_inelPar_RC = (TH1D*)inFile->Get("h_y_inelPar_RC");

  //reconstructed scattered electron
  TH1D *h_eta_scat_ele_RC_eCAL[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_scat_ele_RC_eCAL_E_over_p[nMomBins][nQ2bins][nyInelParBins];

  //momentum distributions in eta bins
  TH1D *h_p_scat_ele_RC_eCAL[nEtaBins+1];
  TH1D *h_p_scat_ele_RC_eCAL_E_over_p[nEtaBins+1];

  //background - negative charged particles
  //TH1D *h_eta_neg_ch_part_RC[nMomBins][nQ2bins][nyInelParBins];

  //negative charged particles after E/p cut
  //TH1D *h_eta_neg_ch_part_RC_eCAL_cut[nMomBins][nQ2bins][nyInelParBins];

  //negative charged par after pfRCIH veto
  //TH1D *h_eta_neg_ch_part_RC_pFRICH_cut[nMomBins][nQ2bins][nyInelParBins];

  //negative charged particles after eCAL+pfRCIH veto
  //TH1D *h_eta_neg_ch_part_RC_eCAL_pfRICH_cut[nMomBins][nQ2bins][nyInelParBins];

  //MC-RC matching
  TH1D *h_scat_ele_purity_eCAL[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_scat_ele_purity_eCAL_E_over_p[nMomBins][nQ2bins][nyInelParBins];

  //with pfRICH
  TH1D* h_momentum_RC_pfRICH = (TH1D*)inFile->Get("h_momentum_RC_pfRICH");

  TH1D* h_Q2_RC_pfRICH = (TH1D*)inFile->Get("h_Q2_RC_pfRICH");

  TH1D *h_y_inelPar_RC_pfRICH = (TH1D*)inFile->Get("h_y_inelPar_RC_pfRICH");

  TH1D *h_eta_scat_ele_RC_lead_p[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_scat_ele_RC_pfRICH[nMomBins][nQ2bins][nyInelParBins];


  TH1D *h_pi_pfRICH_PID_eff_RC[4]; //PID efficiency - no PID baseline, good PID and mis-PID
  TH1D *h_K_pfRICH_PID_eff_RC[4];
  TH1D *h_p_pfRICH_PID_eff_RC[4];

  TH1D *h_pi_pfRICH_PID_eff_MC[4]; //PID efficiency - no PID baseline, good PID and mis-PID
  TH1D *h_K_pfRICH_PID_eff_MC[4];
  TH1D *h_p_pfRICH_PID_eff_MC[4];

  TH1D *h_pi_pfRICH_PID_eff_MC_RC_match[4]; //PID efficiency - no PID baseline, good PID and mis-PID
  TH1D *h_K_pfRICH_PID_eff_MC_RC_match[4];
  TH1D *h_p_pfRICH_PID_eff_MC_RC_match[4];

  TH2D *h_pi_pfRICH_MC_p_RC_p = (TH2D*)inFile->Get("h_pi_pfRICH_MC_p_RC_p");
  TH2D *h_K_pfRICH_MC_p_RC_p = (TH2D*)inFile->Get("h_K_pfRICH_MC_p_RC_p");
  TH2D *h_p_pfRICH_MC_p_RC_p = (TH2D*)inFile->Get("h_p_pfRICH_MC_p_RC_p");



  //TH1D *h_K_purity_pfRICH_MC[nMomBins][nQ2bins][nyInelParBins];
  //TH1D *h_K_purity_pfRICH_RC[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_scat_ele_purity_lead_p[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_scat_ele_purity_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  //momentum distributions in eta bins
  TH1D *h_p_scat_ele_RC_lead_p[nEtaBins+1];
  TH1D *h_p_scat_ele_RC_pfRICH[nEtaBins+1];

  //_____________________________________________________________________________________________


  TCanvas *energy_can = new TCanvas("energy_can", "energy_can", 1200, 1000);

  energy_can->cd();

  h_energy_RC->GetXaxis()->SetTitle("E_{e} (GeV)");
  h_energy_RC->GetXaxis()->CenterTitle();
  h_energy_RC->GetYaxis()->SetTitle("Counts");
  h_energy_RC->GetYaxis()->CenterTitle();
  h_energy_RC->SetMarkerStyle(20);
  h_energy_RC->SetMarkerColor(kRed);
  h_energy_RC->SetLineColor(kRed);
  h_energy_RC->Draw("p e");

  TLegend *Leg_energy = new TLegend(0.3, 0.65, 0.59, 0.75 );
  Leg_energy->SetTextFont(42);
  Leg_energy->AddEntry(h_energy_RC, "Reco. scattered e^{-}");
  Leg_energy->SetBorderSize(0);
  Leg_energy->SetFillColorAlpha(0, 0.01);
  Leg_energy->Draw("same");

  TPaveText *Text_energy = new TPaveText(0.3, 0.45, 0.59, 0.65, "NDC");
  Text_energy->SetTextFont(42);
  Text_energy->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
  Text_energy->AddText("1 < Q^{2} < 20 GeV^{2}/c^{2}");
  Text_energy->AddText("0.01 < y < 0.95");
  //cent_text->AddText("#Lambda^{0} and #bar{#Lambda^{0}}");
  Text_energy->SetFillColorAlpha(0, 0.01);
  Text_energy->Draw("same");

  energy_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Electron_energy_RC.png", e_energy, p_energy));


  TCanvas *E_over_p_can = new TCanvas("E_over_p_can", "E_over_p_can", 1200, 1000);

  E_over_p_can->cd();

  h_E_over_p_RC->GetXaxis()->SetTitle("E/p");
  h_E_over_p_RC->GetXaxis()->CenterTitle();
  h_E_over_p_RC->GetYaxis()->SetTitle("Counts");
  h_E_over_p_RC->GetYaxis()->CenterTitle();
  h_E_over_p_RC->SetMarkerStyle(20);
  h_E_over_p_RC->SetMarkerColor(kRed);
  h_E_over_p_RC->SetLineColor(kRed);
  h_E_over_p_RC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  E_over_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/E_over_p_RC.png", e_energy, p_energy));


  TCanvas *momentum_can = new TCanvas("momentum_can", "momentum_can", 1200, 1000);

  momentum_can->cd();

  h_momentum_RC->GetXaxis()->SetTitle("p (GeV/c)");
  h_momentum_RC->GetXaxis()->CenterTitle();
  h_momentum_RC->GetYaxis()->SetTitle("Counts");
  h_momentum_RC->GetYaxis()->CenterTitle();
  h_momentum_RC->SetMarkerStyle(20);
  h_momentum_RC->SetMarkerColor(kRed);
  h_momentum_RC->SetLineColor(kRed);
  h_momentum_RC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  momentum_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Electron_momentum_RC.png", e_energy, p_energy));



  TCanvas *Q2_can = new TCanvas("Q2_can", "Q2_can", 1200, 1000);

  Q2_can->cd();

  h_Q2_RC->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_Q2_RC->GetXaxis()->CenterTitle();
  h_Q2_RC->GetYaxis()->SetTitle("Counts");
  h_Q2_RC->GetYaxis()->CenterTitle();
  h_Q2_RC->SetMarkerStyle(20);
  h_Q2_RC->SetMarkerColor(kRed);
  h_Q2_RC->SetLineColor(kRed);
  h_Q2_RC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  Q2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Electron_Q2_RC.png", e_energy, p_energy));


  TCanvas *y_inelPar_can = new TCanvas("y_inelPar_can", "y_inelPar_can", 1200, 1000);

  y_inelPar_can->cd();

  h_y_inelPar_RC->GetXaxis()->SetTitle("y");
  h_y_inelPar_RC->GetXaxis()->CenterTitle();
  h_y_inelPar_RC->GetYaxis()->SetTitle("Counts");
  h_y_inelPar_RC->GetYaxis()->CenterTitle();
  h_y_inelPar_RC->SetMarkerStyle(20);
  h_y_inelPar_RC->SetMarkerColor(kRed);
  h_y_inelPar_RC->SetLineColor(kRed);
  h_y_inelPar_RC->SetMinimum(0);
  h_y_inelPar_RC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  y_inelPar_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Electron_y_inelPar_RC.png", e_energy, p_energy));

  //_________________________________________________________________________________________________

  //pfRICH
  TCanvas *momentum_pfRICH_can = new TCanvas("momentum_pfRICH_can", "momentum_pfRICH_can", 1200, 1000);

  momentum_pfRICH_can->cd();

  h_momentum_RC_pfRICH->GetXaxis()->SetTitle("p (GeV/c)");
  h_momentum_RC_pfRICH->GetXaxis()->CenterTitle();
  h_momentum_RC_pfRICH->GetYaxis()->SetTitle("Counts");
  h_momentum_RC_pfRICH->GetYaxis()->CenterTitle();
  h_momentum_RC_pfRICH->SetMarkerStyle(20);
  h_momentum_RC_pfRICH->SetMarkerColor(kRed);
  h_momentum_RC_pfRICH->SetLineColor(kRed);
  h_momentum_RC_pfRICH->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  momentum_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/Electron_momentum_RC_pfRICH.png", e_energy, p_energy));



  TCanvas *Q2_pfRICH_can = new TCanvas("Q2_pfRICH_can", "Q2_pfRICH_can", 1200, 1000);

  Q2_pfRICH_can->cd();

  h_Q2_RC_pfRICH->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_Q2_RC_pfRICH->GetXaxis()->CenterTitle();
  h_Q2_RC_pfRICH->GetYaxis()->SetTitle("Counts");
  h_Q2_RC_pfRICH->GetYaxis()->CenterTitle();
  h_Q2_RC_pfRICH->SetMarkerStyle(20);
  h_Q2_RC_pfRICH->SetMarkerColor(kRed);
  h_Q2_RC_pfRICH->SetLineColor(kRed);
  h_Q2_RC_pfRICH->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  Q2_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/Electron_Q2_RC_pfRICH.png", e_energy, p_energy));


  TCanvas *y_inelPar_pfRICH_can = new TCanvas("y_inelPar_pfRICH_can", "y_inelPar_pfRICH_can", 1200, 1000);

  y_inelPar_pfRICH_can->cd();

  h_y_inelPar_RC_pfRICH->GetXaxis()->SetTitle("y");
  h_y_inelPar_RC_pfRICH->GetXaxis()->CenterTitle();
  h_y_inelPar_RC_pfRICH->GetYaxis()->SetTitle("Counts");
  h_y_inelPar_RC_pfRICH->GetYaxis()->CenterTitle();
  h_y_inelPar_RC_pfRICH->SetMarkerStyle(20);
  h_y_inelPar_RC_pfRICH->SetMarkerColor(kRed);
  h_y_inelPar_RC_pfRICH->SetLineColor(kRed);
  h_y_inelPar_RC_pfRICH->SetMinimum(0);
  h_y_inelPar_RC_pfRICH->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  y_inelPar_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/Electron_y_inelPar_RC_pfRICH.png", e_energy, p_energy));

  //___________________________________________________________________________________________________________


  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {
        //MC histograms
        TCanvas *eta_scat_ele_can = new TCanvas(Form("eta_scat_ele_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), Form("eta_scat_ele_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), 2000, 1500);

        eta_scat_ele_can->cd();


        h_eta_scat_ele[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));


        h_eta_pi_minus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_pfRICH_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_eCAL_85_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_eCAL_95_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_K_minus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_K_minus_pfRICH_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_anti_proton[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_anti_proton_pfRICH_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        if(h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Integral() == 0 ) continue;
        if( h_eta_pi_minus[mom_bin][Q2bin][y_bin]->Integral() < 15 ) continue;

        if(h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Integral() > 0 && h_eta_pi_minus[mom_bin][Q2bin][y_bin]->Integral() > 15)
        {
          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("p e");

          //cout<<h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->GetMaximum()<<endl;

          TCanvas *eta_scat_ele_pi_minus_can = new TCanvas(Form("eta_scat_ele_pi_minus_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), Form("eta_scat_ele_pi_minus_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), 2000, 1500);

          eta_scat_ele_pi_minus_can->cd();

          gPad->SetLogy();

          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
          //for lin scale
          //if( h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum() < h_eta_pi_minus[mom_bin][Q2bin][y_bin]->GetMaximum() ) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0., h_eta_pi_minus[mom_bin][Q2bin][y_bin]->GetMaximum()*1.7);
          //else h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0.5, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum()*1.7);
          //for log scale

          if(h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum() != 0 )
          {
            h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-8, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum()*1e8);
          }
          else if(h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->GetMaximum() != 0)
          {
            h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-8, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum()*1e8);
          }
          else
          {
            h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_pi_minus[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-8, h_eta_pi_minus[mom_bin][Q2bin][y_bin]->GetMaximum()*1e10);
          }
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerColor(kRed);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetLineColor(kRed);
          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMinimum(0);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("p e");


          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(25);
          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->SetMarkerStyle(24);
          h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->SetMarkerColor(8);
          h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->SetLineColor(8);
          h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin]->SetMarkerStyle(28);
          h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin]->SetMarkerColor(6);
          h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin]->SetLineColor(6);
          h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(47);
          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(kBlue);
          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(kBlue);
          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(33);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(kOrange+2);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(kOrange+2);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e same");




          TLegend *Leg_bins_anti = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_anti->SetTextFont(42);
          Leg_bins_anti->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "Scattered e^{-}");
          Leg_bins_anti->AddEntry(h_eta_pi_minus[mom_bin][Q2bin][y_bin], "#pi^{-}");
          Leg_bins_anti->AddEntry(h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin], "#pi^{-} eCAL 85% eff.");
          Leg_bins_anti->AddEntry(h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin], "#pi^{-} eCAL 95% eff.");
          Leg_bins_anti->AddEntry(h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin], "#pi^{-} eCAL(85%)+pfRICH");
          Leg_bins_anti->AddEntry(h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin], "#pi^{-} eCAL(95%)+pfRICH");
          Leg_bins_anti->SetBorderSize(0);
          Leg_bins_anti->SetFillColorAlpha(0, 0.01);
          Leg_bins_anti->Draw("same");

          TPaveText *Text_bins = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_bins->SetTextFont(42);
          Text_bins->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          if( mom_bin < nMomBins )
          {
            Text_bins->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_bins->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_bins->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }

          Text_bins->SetFillColorAlpha(0, 0.01);
          Text_bins->Draw("same");

          eta_scat_ele_pi_minus_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject/Eta_scat_e_vs_pi_minus_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));



          TCanvas *eta_scat_ele_pi_minus_pfRICH_can = new TCanvas(Form("eta_scat_ele_pi_minus_pfRICH_can_pfRICH_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), Form("eta_scat_ele_pi_minus_pfRICH_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), 2000, 1500);

          eta_scat_ele_pi_minus_pfRICH_can->cd();

          gPad->SetLogy();


          if( h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum() != 0) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-4, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum());
          else if( h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum() != 0) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-4, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum());
          else if( h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum() != 0) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-4, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum());
          else h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(h_eta_pi_minus[mom_bin][Q2bin][y_bin]->GetMaximum()*1e-4, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum() );


          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("p e");


          h_eta_pi_minus[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(21);
          h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(8);
          h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(8);
          h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(28);
          h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerColor(6);
          h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetLineColor(6);
          h_eta_K_minus[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(34);
          h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(kBlue);
          h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(kBlue);
          h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e same");

          h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerStyle(27);
          h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerColor(kOrange+2);
          h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetLineColor(kOrange+2);
          h_eta_anti_proton[mom_bin][Q2bin][y_bin]->Draw("p e same");


          h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(33);
          h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e same");




          TLegend *Leg_bins_anti_pfRICH = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_anti_pfRICH->SetTextFont(42);
          Leg_bins_anti_pfRICH->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "Scattered e^{-}");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_pi_minus[mom_bin][Q2bin][y_bin], "#pi^{-}");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin], "#pi^{-} pfRICH");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_K_minus[mom_bin][Q2bin][y_bin], "K^{-}");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin], "K^{-} pfRICH");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_anti_proton[mom_bin][Q2bin][y_bin], "#bar{p}");
          Leg_bins_anti_pfRICH->AddEntry(h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin], "#bar{p} pfRICH");
          Leg_bins_anti_pfRICH->SetBorderSize(0);
          Leg_bins_anti_pfRICH->SetFillColorAlpha(0, 0.01);
          Leg_bins_anti_pfRICH->Draw("same");

          Text_bins->Draw("same");

          eta_scat_ele_pi_minus_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject/Eta_scat_e_vs_pi_minus_pfRICH_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));
        }
        //_____________________________________________________________________________________________________________________




        //reco particles histograms
        //selection using eCAL

        //leading E cluster in eCAL
        h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_RC_eCAL_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));


        //background - negative charged particles
        //h_eta_neg_ch_part_RC[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_neg_ch_part_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        //mis-PID with eCAL
        h_scat_ele_purity_eCAL[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_scat_ele_purity_eCAL_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        if(h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->Integral() > 0)
        {

          TCanvas *eta_scat_ele_RC_eCAL_can = new TCanvas(Form("eta_scat_ele_RC_eCAL_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("eta_scat_ele_RC_eCAL_can_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          eta_scat_ele_RC_eCAL_can->cd();

          gPad->SetLogy();

          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetMaximum(h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->GetMaximum()*100);
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          //h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->SetMinimum(0);
          h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin]->Draw("p e");


          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(1);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("same p e ");


          TLegend *Leg_bins_RC = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_RC->SetTextFont(42);
          Leg_bins_RC->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "MC Scattered e^{-}");
          Leg_bins_RC->AddEntry(h_eta_scat_ele_RC_eCAL[mom_bin][Q2bin][y_bin], "Leading E cluster in eCAL");
          Leg_bins_RC->SetBorderSize(0);
          Leg_bins_RC->SetFillColorAlpha(0, 0.01);
          Leg_bins_RC->Draw("same");

          TPaveText *Text_bins_2 = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_bins_2->SetTextFont(42);
          Text_bins_2->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          if( mom_bin < nMomBins )
          {
            Text_bins_2->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_bins_2->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_bins_2->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }

          Text_bins_2->SetFillColorAlpha(0, 0.01);
          Text_bins_2->Draw("same");

          eta_scat_ele_RC_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Eta_scat_e_RC_eCAL_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


          //scattered electron purity histograms
          TCanvas *scat_ele_purity_eCAL_can = new TCanvas(Form("scat_ele_purity_eCAL_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("scat_ele_purity_eCAL_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          scat_ele_purity_eCAL_can->cd();

          h_scat_ele_purity_eCAL[mom_bin][Q2bin][y_bin]->Draw();

          double nElectrons = h_scat_ele_purity_eCAL[mom_bin][Q2bin][y_bin]->GetBinContent(2);
          double nBackground = h_scat_ele_purity_eCAL[mom_bin][Q2bin][y_bin]->GetBinContent(1);

          double e_purity = ( nElectrons/(nElectrons+nBackground) )*100;

          TPaveText *Text_match = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_match->SetTextFont(42);
          Text_match->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          Text_match->AddText(Form("Scat. e purity: %0.1f%%", e_purity));
          if( mom_bin < nMomBins )
          {
            Text_match->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_match->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_match->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }
          Text_match->SetFillColorAlpha(0, 0.01);
          Text_match->Draw("same");

          scat_ele_purity_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/scat_ele_purity_eCAL_can_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


        }//end if integral
        //_______________________________________________________________________________________________________________________


        //reconstructed scattered electron
        h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_RC_eCAL_E_over_p_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));


        //background - negative charged particles
        //h_eta_neg_ch_part_RC[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_neg_ch_part_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        //mis-PID with eCAL_E_over_p
        h_scat_ele_purity_eCAL_E_over_p[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_scat_ele_purity_eCAL_E_over_p_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        if(h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->Integral() > 0)
        {

          TCanvas *eta_scat_ele_RC_eCAL_E_over_p_can = new TCanvas(Form("eta_scat_ele_RC_eCAL_E_over_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("eta_scat_ele_RC_eCAL_E_over_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          eta_scat_ele_RC_eCAL_E_over_p_can->cd();

          gPad->SetLogy();

          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetMaximum(h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetMaximum()*100);
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          //h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->SetMinimum(0);
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->Draw("p e");


          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(1);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("same p e ");


          TLegend *Leg_bins_RC = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_RC->SetTextFont(42);
          Leg_bins_RC->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "MC Scattered e^{-}");
          Leg_bins_RC->AddEntry(h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin][Q2bin][y_bin], "Leading E cluster in eCAL_E_over_p");
          Leg_bins_RC->SetBorderSize(0);
          Leg_bins_RC->SetFillColorAlpha(0, 0.01);
          Leg_bins_RC->Draw("same");

          TPaveText *Text_bins_2 = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_bins_2->SetTextFont(42);
          Text_bins_2->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          if( mom_bin < nMomBins )
          {
            Text_bins_2->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_bins_2->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_bins_2->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }

          Text_bins_2->SetFillColorAlpha(0, 0.01);
          Text_bins_2->Draw("same");

          eta_scat_ele_RC_eCAL_E_over_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/Eta_scat_e_RC_eCAL_E_over_p_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


          //scattered electron purity histograms
          TCanvas *scat_ele_purity_eCAL_E_over_p_can = new TCanvas(Form("scat_ele_purity_eCAL_E_over_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("scat_ele_purity_eCAL_E_over_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          scat_ele_purity_eCAL_E_over_p_can->cd();

          h_scat_ele_purity_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->Draw();

          double nElectrons = h_scat_ele_purity_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetBinContent(2);
          double nBackground = h_scat_ele_purity_eCAL_E_over_p[mom_bin][Q2bin][y_bin]->GetBinContent(1);

          double e_purity = ( nElectrons/(nElectrons+nBackground) )*100;

          TPaveText *Text_match = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_match->SetTextFont(42);
          Text_match->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          Text_match->AddText(Form("Scat. e purity: %0.1f%%", e_purity));
          if( mom_bin < nMomBins )
          {
            Text_match->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_match->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_match->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }
          Text_match->SetFillColorAlpha(0, 0.01);
          Text_match->Draw("same");

          scat_ele_purity_eCAL_E_over_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco/scat_ele_purity_eCAL_E_over_p_can_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


        }//end if integral
        //_______________________________________________________________________________________________________________________



        //leading momentum tracks

        h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_RC_lead_p_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        //purity of leading momentum sample
        h_scat_ele_purity_lead_p[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_scat_ele_purity_lead_p_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        if(h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->Integral() > 0)
        {
          TCanvas *eta_scat_ele_RC_lead_p_can = new TCanvas(Form("eta_scat_ele_RC_lead_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("eta_scat_ele_RC_lead_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          eta_scat_ele_RC_lead_p_can->cd();

          gPad->SetLogy();

          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetMaximum(h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->GetMaximum()*500);
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          //h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->SetMinimum(0);
          h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin]->Draw("p e");


          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(1);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("same p e ");


          TLegend *Leg_bins_RC_lead_p = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_RC_lead_p->SetTextFont(42);
          Leg_bins_RC_lead_p->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "MC Scattered e^{-}");
          Leg_bins_RC_lead_p->AddEntry(h_eta_scat_ele_RC_lead_p[mom_bin][Q2bin][y_bin], "Lead. p RC particles");
          Leg_bins_RC_lead_p->SetBorderSize(0);
          Leg_bins_RC_lead_p->SetFillColorAlpha(0, 0.01);
          Leg_bins_RC_lead_p->Draw("same");

          TPaveText *Text_bins_3 = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_bins_3->SetTextFont(42);
          Text_bins_3->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          if( mom_bin < nMomBins )
          {
            Text_bins_3->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_bins_3->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_bins_3->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }

          Text_bins_3->SetFillColorAlpha(0, 0.01);
          Text_bins_3->Draw("same");

          eta_scat_ele_RC_lead_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/Eta_scat_e_RC_lead_p_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


          //scattered electron purity histograms
          TCanvas *scat_ele_purity_lead_p_can = new TCanvas(Form("scat_ele_purity_lead_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("scat_ele_purity_lead_p_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          scat_ele_purity_lead_p_can->cd();

          h_scat_ele_purity_lead_p[mom_bin][Q2bin][y_bin]->Draw();

          double nElectrons = h_scat_ele_purity_lead_p[mom_bin][Q2bin][y_bin]->GetBinContent(2);
          double nBackground = h_scat_ele_purity_lead_p[mom_bin][Q2bin][y_bin]->GetBinContent(1);

          double e_purity = ( nElectrons/(nElectrons+nBackground) )*100;

          TPaveText *Text_mathc_lead_p = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_mathc_lead_p->SetTextFont(42);
          Text_mathc_lead_p->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          Text_mathc_lead_p->AddText(Form("Scat. e purity: %0.1f%%", e_purity));
          if( mom_bin < nMomBins )
          {
            Text_mathc_lead_p->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_mathc_lead_p->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_mathc_lead_p->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }
          Text_mathc_lead_p->SetFillColorAlpha(0, 0.01);
          Text_mathc_lead_p->Draw("same");

          scat_ele_purity_lead_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/scat_ele_purity_lead_p_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


        }//end if integral
        //_______________________________________________________________________________________________________________________


        //selection using pfRICH

        h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        //mis-PID with eCAL
        h_scat_ele_purity_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_scat_ele_purity_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));

        if(h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->Integral() > 0)
        {
          TCanvas *eta_scat_ele_RC_pfRICH_can = new TCanvas(Form("eta_scat_ele_RC_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("eta_scat_ele_RC_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          eta_scat_ele_RC_pfRICH_can->cd();

          gPad->SetLogy();

          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetMaximum(h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->GetMaximum()*500);
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetLineColor(1);
          //h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->SetMinimum(0);
          h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin]->Draw("p e");


          //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(1);
          h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("same p e ");


          TLegend *Leg_bins_RC_pfRICH = new TLegend(0.2, 0.63, 0.54, 0.85 );
          Leg_bins_RC_pfRICH->SetTextFont(42);
          Leg_bins_RC_pfRICH->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "MC Scattered e^{-}");
          Leg_bins_RC_pfRICH->AddEntry(h_eta_scat_ele_RC_pfRICH[mom_bin][Q2bin][y_bin], "RC Scattered e^{-}");
          Leg_bins_RC_pfRICH->SetBorderSize(0);
          Leg_bins_RC_pfRICH->SetFillColorAlpha(0, 0.01);
          Leg_bins_RC_pfRICH->Draw("same");

          TPaveText *Text_bins_3 = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_bins_3->SetTextFont(42);
          Text_bins_3->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          if( mom_bin < nMomBins )
          {
            Text_bins_3->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_bins_3->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_bins_3->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }

          Text_bins_3->SetFillColorAlpha(0, 0.01);
          Text_bins_3->Draw("same");

          eta_scat_ele_RC_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/Eta_scat_e_RC_pfRICH_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


          //scattered electron purity histograms
          TCanvas *scat_ele_purity_pfRICH_can = new TCanvas(Form("scat_ele_purity_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("scat_ele_purity_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

          scat_ele_purity_pfRICH_can->cd();

          h_scat_ele_purity_pfRICH[mom_bin][Q2bin][y_bin]->Draw();

          double nElectrons = h_scat_ele_purity_pfRICH[mom_bin][Q2bin][y_bin]->GetBinContent(2);
          double nBackground = h_scat_ele_purity_pfRICH[mom_bin][Q2bin][y_bin]->GetBinContent(1);

          double e_purity = ( nElectrons/(nElectrons+nBackground) )*100;

          TPaveText *Text_mathc_pfRICH = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
          Text_mathc_pfRICH->SetTextFont(42);
          Text_mathc_pfRICH->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
          Text_mathc_pfRICH->AddText(Form("Scat. e purity: %0.1f%%", e_purity));
          if( mom_bin < nMomBins )
          {
            Text_mathc_pfRICH->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
            Text_mathc_pfRICH->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
            Text_mathc_pfRICH->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
          }
          Text_mathc_pfRICH->SetFillColorAlpha(0, 0.01);
          Text_mathc_pfRICH->Draw("same");

          scat_ele_purity_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/scat_ele_purity_pfRICH_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));


        }//end if integral
        //_______________________________________________________________________________________________________________________


      }//end for y bins
    }//end for Q2 bins
  }//end for mom bins


  for(unsigned int eta_bin = 0; eta_bin < nEtaBins+1; eta_bin++)
  {
    //pur MC
    //direct comparison of scattered e and pi- p spectra
    TCanvas *p_scat_ele_can = new TCanvas(Form("p_scat_ele_can_eta_%i", eta_bin), Form("p_scat_ele_can_eta_%i", eta_bin), 1200, 1000);
    p_scat_ele_can->cd();

    gPad->SetLogy();

    h_p_scat_ele[eta_bin] = (TH1D*)inFile->Get(Form("h_p_scat_ele_eta_%i", eta_bin));
    h_p_scat_ele[eta_bin]->Sumw2();
    h_p_scat_ele[eta_bin]->GetXaxis()->SetTitle("p (GeV)/#it{c}");
    h_p_scat_ele[eta_bin]->GetXaxis()->CenterTitle();
    h_p_scat_ele[eta_bin]->GetYaxis()->SetTitle("Counts");
    h_p_scat_ele[eta_bin]->GetYaxis()->CenterTitle();
    h_p_scat_ele[eta_bin]->SetMarkerStyle(20);
    h_p_scat_ele[eta_bin]->SetMarkerSize(2);
    h_p_scat_ele[eta_bin]->SetMarkerColor(kRed);
    h_p_scat_ele[eta_bin]->SetLineColor(kRed);
    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_pi_minus[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_eta_%i", eta_bin));
    h_p_pi_minus[eta_bin]->Sumw2();
    h_p_pi_minus[eta_bin]->SetMarkerStyle(25);
    h_p_pi_minus[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus[eta_bin]->SetLineColor(1);
    h_p_pi_minus[eta_bin]->Draw("p e same");

    p_scat_ele_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_scat_ele_eCAL_can = new TCanvas(Form("p_scat_ele_eCAL_can_eta_%i", eta_bin), Form("p_scat_ele_eCAL_can_eta_%i", eta_bin), 1200, 1000);
    p_scat_ele_eCAL_can->cd();

    gPad->SetLogy();

    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_pi_minus_eCAL_85[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_eCAL_85_eta_%i", eta_bin));
    h_p_pi_minus_eCAL_85[eta_bin]->Sumw2();
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerStyle(24);
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerColor(8);
    h_p_pi_minus_eCAL_85[eta_bin]->SetLineColor(8);
    h_p_pi_minus_eCAL_85[eta_bin]->Draw("p e same");

    h_p_pi_minus_eCAL_95[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_eCAL_95_eta_%i", eta_bin));
    h_p_pi_minus_eCAL_95[eta_bin]->Sumw2();
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerStyle(28);
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerColor(6);
    h_p_pi_minus_eCAL_95[eta_bin]->SetLineColor(6);
    h_p_pi_minus_eCAL_95[eta_bin]->Draw("p e same");

    p_scat_ele_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_eCAL_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_scat_ele_pfRICH_can = new TCanvas(Form("p_scat_ele_pfRICH_can_eta_%i", eta_bin), Form("p_scat_ele_pfRICH_can_eta_%i", eta_bin), 1200, 1000);
    p_scat_ele_pfRICH_can->cd();

    gPad->SetLogy();

    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_pi_minus_pfRICH[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_pfRICH_eta_%i", eta_bin));
    h_p_pi_minus_pfRICH[eta_bin]->Sumw2();
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerStyle(21);
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerColor(8);
    h_p_pi_minus_pfRICH[eta_bin]->SetLineColor(8);
    h_p_pi_minus_pfRICH[eta_bin]->Draw("p e same");

    p_scat_ele_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_pfRICH_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_scat_ele_pfRICH_eCAL_can = new TCanvas(Form("p_scat_ele_pfRICH_eCAL_can_eta_%i", eta_bin), Form("p_scat_ele_pfRICH_eCAL_can_eta_%i", eta_bin), 1200, 1000);
    p_scat_ele_pfRICH_eCAL_can->cd();

    gPad->SetLogy();

    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_pi_minus_eCAL_85_pfRICH[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_eCAL_85_pfRICH_eta_%i", eta_bin));
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->Sumw2();
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerStyle(24);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerColor(8);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetLineColor(8);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->Draw("p e same");


    h_p_pi_minus_eCAL_95_pfRICH[eta_bin] = (TH1D*)inFile->Get(Form("h_p_pi_minus_eCAL_95_pfRICH_eta_%i", eta_bin));
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->Sumw2();
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerStyle(28);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerColor(6);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetLineColor(6);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->Draw("p e same");

    p_scat_ele_pfRICH_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_pfRICH_eCAL_eta_%i.png", e_energy, p_energy, eta_bin));

    //____________________________________________

    // pi/e ratios
    TCanvas *p_pi_over_scat_ele_can = new TCanvas(Form("p_pi_over_scat_ele_can_eta_%i", eta_bin), Form("p_pi_over_scat_ele_can_eta_%i", eta_bin), 1200, 1000);
    p_pi_over_scat_ele_can->cd();

    gPad->SetLogy();

    h_p_pi_minus[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus[eta_bin]->GetXaxis()->SetTitle("p (GeV)/#it{c}");
    h_p_pi_minus[eta_bin]->GetXaxis()->CenterTitle();
    h_p_pi_minus[eta_bin]->GetYaxis()->SetTitle("#pi^{-}/scat e");
    h_p_pi_minus[eta_bin]->GetYaxis()->CenterTitle();
    h_p_pi_minus[eta_bin]->SetMarkerStyle(21);
    h_p_pi_minus[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus[eta_bin]->SetLineColor(1);
    h_p_pi_minus[eta_bin]->Draw("p e");

    p_pi_over_scat_ele_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_pi_over_scat_ele_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_pi_over_scat_ele_eCAL_can = new TCanvas(Form("p_pi_over_scat_ele_eCAL_can_eta_%i", eta_bin), Form("p_pi_over_scat_ele_eCAL_can_eta_%i", eta_bin), 1200, 1000);
    p_pi_over_scat_ele_eCAL_can->cd();

    gPad->SetLogy();

    h_p_pi_minus_eCAL_85[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus_eCAL_85[eta_bin]->GetXaxis()->SetTitle("p (GeV)/#it{c}");
    h_p_pi_minus_eCAL_85[eta_bin]->GetXaxis()->CenterTitle();
    h_p_pi_minus_eCAL_85[eta_bin]->GetYaxis()->SetTitle("#pi^{-}/scat e");
    h_p_pi_minus_eCAL_85[eta_bin]->GetYaxis()->CenterTitle();
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerStyle(21);
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_85[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus_eCAL_85[eta_bin]->SetLineColor(1);
    h_p_pi_minus_eCAL_85[eta_bin]->Draw("p e");


    h_p_pi_minus_eCAL_95[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerStyle(24);
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_95[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus_eCAL_95[eta_bin]->SetLineColor(1);
    h_p_pi_minus_eCAL_95[eta_bin]->Draw("p e same");


    p_pi_over_scat_ele_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_pi_over_scat_ele_eCAL_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_pi_over_scat_ele_pfRICH_can = new TCanvas(Form("p_pi_over_scat_ele_pfRICH_can_eta_%i", eta_bin), Form("p_pi_over_scat_ele_pfRICH_can_eta_%i", eta_bin), 1200, 1000);
    p_pi_over_scat_ele_pfRICH_can->cd();

    gPad->SetLogy();

    h_p_pi_minus_pfRICH[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus_pfRICH[eta_bin]->GetXaxis()->SetTitle("p (GeV)/#it{c}");
    h_p_pi_minus_pfRICH[eta_bin]->GetXaxis()->CenterTitle();
    h_p_pi_minus_pfRICH[eta_bin]->GetYaxis()->SetTitle("#pi^{-}/scat e");
    h_p_pi_minus_pfRICH[eta_bin]->GetYaxis()->CenterTitle();
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerStyle(21);
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_pfRICH[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus_pfRICH[eta_bin]->SetLineColor(1);
    h_p_pi_minus_pfRICH[eta_bin]->Draw("p e");

    p_pi_over_scat_ele_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_pi_over_scat_ele_pfRICH_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_pi_over_scat_ele_pfRICH_eCAL_can = new TCanvas(Form("p_pi_over_scat_ele_pfRICH_eCAL_can_eta_%i", eta_bin), Form("p_pi_over_scat_ele_pfRICH_eCAL_can_eta_%i", eta_bin), 1200, 1000);
    p_pi_over_scat_ele_pfRICH_eCAL_can->cd();

    gPad->SetLogy();

    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->GetXaxis()->SetTitle("p (GeV)/#it{c}");
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->GetXaxis()->CenterTitle();
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->GetYaxis()->SetTitle("#pi^{-}/scat e");
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->GetYaxis()->CenterTitle();
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerStyle(21);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->SetLineColor(1);
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->Draw("p e");


    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->Divide(h_p_scat_ele[eta_bin]);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerStyle(24);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetMarkerColor(1);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->SetLineColor(1);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->Draw("p e same");

    p_pi_over_scat_ele_pfRICH_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_pi_over_scat_ele_pfRICH_eCAL_eta_%i.png", e_energy, p_energy, eta_bin));

    //____________________________________________


    //eCAL scattered electron selection
    TCanvas *p_scat_ele_RC_eCAL_can = new TCanvas(Form("p_scat_ele_RC_eCAL_can_%i", eta_bin), Form("p_scat_ele_RC_eCAL_can_%i", eta_bin), 1200, 1000);
    p_scat_ele_RC_eCAL_can->cd();

    gPad->SetLogy();

    //MC scattered e
    h_p_scat_ele[eta_bin]->Draw("p e");

    //leading cluster tracks from eCAL
    h_p_scat_ele_RC_eCAL[eta_bin] = (TH1D*)inFile->Get(Form("h_p_scat_ele_RC_eCAL_eta_%i", eta_bin));
    h_p_scat_ele_RC_eCAL[eta_bin]->SetMarkerStyle(20);
    h_p_scat_ele_RC_eCAL[eta_bin]->SetMarkerSize(2);
    h_p_scat_ele_RC_eCAL[eta_bin]->SetMarkerColor(1);
    h_p_scat_ele_RC_eCAL[eta_bin]->SetLineColor(1);
    h_p_scat_ele_RC_eCAL[eta_bin]->Draw("p e same");

    p_scat_ele_RC_eCAL_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_RC_eCAL_eta_%i.png", e_energy, p_energy, eta_bin));


    TCanvas *p_scat_ele_RC_eCAL_E_over_p_can = new TCanvas(Form("p_scat_ele_RC_eCAL_E_over_p_can_%i", eta_bin), Form("p_scat_ele_RC_eCAL_E_over_p_can_%i", eta_bin), 1200, 1000);
    p_scat_ele_RC_eCAL_E_over_p_can->cd();

    gPad->SetLogy();

    //MC scattered e
    h_p_scat_ele[eta_bin]->Draw("p e");

    //leading cluster tracks from eCAL with E/p cut
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin] = (TH1D*)inFile->Get(Form("h_p_scat_ele_RC_eCAL_E_over_p_eta_%i", eta_bin));
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin]->SetMarkerStyle(20);
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin]->SetMarkerSize(2);
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin]->SetMarkerColor(1);
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin]->SetLineColor(1);
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin]->Draw("p e same");

    p_scat_ele_RC_eCAL_E_over_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_RC_eCAL_E_over_p_eta_%i.png", e_energy, p_energy, eta_bin));

    //______________________________________________________________________

    //pfRICH/tracking scattered electron selection
    TCanvas *p_scat_ele_RC_lead_p_can = new TCanvas(Form("p_scat_ele_RC_lead_p_can_%i", eta_bin), Form("p_scat_ele_RC_lead_p_can_%i", eta_bin), 1200, 1000);
    p_scat_ele_RC_lead_p_can->cd();

    gPad->SetLogy();

    //MC scattered e
    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_scat_ele_RC_lead_p[eta_bin] = (TH1D*)inFile->Get(Form("h_p_scat_ele_RC_lead_p_eta_%i", eta_bin));
    h_p_scat_ele_RC_lead_p[eta_bin]->SetMarkerStyle(20);
    h_p_scat_ele_RC_lead_p[eta_bin]->SetMarkerSize(2);
    h_p_scat_ele_RC_lead_p[eta_bin]->SetMarkerColor(1);
    h_p_scat_ele_RC_lead_p[eta_bin]->SetLineColor(1);
    h_p_scat_ele_RC_lead_p[eta_bin]->Draw("p e same");

    p_scat_ele_RC_lead_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_RC_lead_p_eta_%i.png", e_energy, p_energy, eta_bin));



    TCanvas *p_scat_ele_RC_pfRICH_can = new TCanvas(Form("p_scat_ele_RC_pfRICH_can_%i", eta_bin), Form("p_scat_ele_RC_pfRICH_can_%i", eta_bin), 1200, 1000);
    p_scat_ele_RC_pfRICH_can->cd();

    gPad->SetLogy();

    //MC scattered e
    h_p_scat_ele[eta_bin]->Draw("p e");

    h_p_scat_ele_RC_pfRICH[eta_bin] = (TH1D*)inFile->Get(Form("h_p_scat_ele_RC_pfRICH_eta_%i", eta_bin));
    h_p_scat_ele_RC_pfRICH[eta_bin]->SetMarkerStyle(20);
    h_p_scat_ele_RC_pfRICH[eta_bin]->SetMarkerSize(2);
    h_p_scat_ele_RC_pfRICH[eta_bin]->SetMarkerColor(1);
    h_p_scat_ele_RC_pfRICH[eta_bin]->SetLineColor(1);
    h_p_scat_ele_RC_pfRICH[eta_bin]->Draw("p e same");

    p_scat_ele_RC_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/p_scat_ele_RC_pfRICH_eta_%i.png", e_energy, p_energy, eta_bin));

    //______________________________________________________________________

  }




  //fitPiTOF->Write();
  //fitKTOF->Write();

  inFile->Close();

  return;

}
