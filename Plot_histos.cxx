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
/*
{"0.1": {"0.85": 0.6332832672698294, "0.95": 0.8838860654090215},
"0.2": {"0.85": 0.7495818158985306, "0.95": 0.9228366502089709},
"0.5": {"0.85": 0.00930384575910461, "0.95": 0.02228665375203912},
"1.0": {"0.85": 0.001692827694491846, "0.95": 0.0036237053948322794},
"2.0": {"0.85": 0.0001898238241173789, "0.95":0.00048096291971113834},
"5.0": {"0.85": 0.00020018016214593134,"0.95": 0.0006523112180272588},
"10.0": {"0.85": 0.000536412900269677, "0.95": 0.0022770026949322946},
"20": {"0.85": 0.0006430092230696459, "0.95": 0.0018829746368912276}}
*/
void Plot_histos(int e_energy = 18, int p_energy = 275)
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
  //draw pi rejection plot

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
  eCAL_leg->AddEntry(g_pi_false_rate_85, "eCAL 85% eff.");
  eCAL_leg->AddEntry(g_pi_false_rate_95, "eCAL 95% eff.");
  eCAL_leg->SetBorderSize(0);
  eCAL_leg->SetFillColorAlpha(0, 0.01);
  eCAL_leg->Draw("same");

  pi_false_rate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/pi_false_rate.png");

  //___________________________________________________


  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 8;
  float const mom_bins[nMomBins+1] = { 0,0.5,1,1.5,2,3,7,10, 18 };

  //load all files
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output.root", e_energy, p_energy, e_energy, p_energy), "READ");
  //TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i_p_scat_bins-output.root", e_energy, p_energy, e_energy, p_energy), "READ");

  TH1D* h_energy_MC = (TH1D*)inFile->Get("h_energy_MC");
  TH1D* h_energy_zoom_MC = (TH1D*)inFile->Get("h_energy_zoom_MC");

  TH1D* h_momentum_MC = (TH1D*)inFile->Get("h_momentum_MC");

  TH1D* h_Q2_MC = (TH1D*)inFile->Get("h_Q2_MC");
  TH1D* h_Q2_zoom_MC = (TH1D*)inFile->Get("h_Q2_zoom_MC");

  TH1D *h_y_inelPar_MC = (TH1D*)inFile->Get("h_y_inelPar_MC");
  TH1D *h_y_inelPar_zoom_MC = (TH1D*)inFile->Get("h_y_inelPar_zoom_MC");

  //eta distributions in multiple Q^2 and inelasticity bins
  TH1D *h_eta_scat_ele[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_ele[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_plus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_K_plus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_proton[nMomBins][nQ2bins][nyInelParBins];

  TH1D *h_eta_positron[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_K_minus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_anti_proton[nMomBins][nQ2bins][nyInelParBins];
  //_____________________________________________________________________________________________



  TCanvas *energy_can = new TCanvas("energy_can", "energy_can", 1200, 1000);

  energy_can->cd();

  h_energy_MC->GetXaxis()->SetTitle("E_{e} (GeV)");
  h_energy_MC->GetXaxis()->CenterTitle();
  h_energy_MC->GetYaxis()->SetTitle("Counts");
  h_energy_MC->GetYaxis()->CenterTitle();
  h_energy_MC->SetMarkerStyle(20);
  h_energy_MC->SetMarkerColor(kRed);
  h_energy_MC->SetLineColor(kRed);
  h_energy_MC->Draw("p e");

  TLegend *Leg_energy = new TLegend(0.5, 0.65, 0.79, 0.75 );
  Leg_energy->SetTextFont(42);
  Leg_energy->AddEntry(h_energy_MC, "Scattered e^{-}");
  Leg_energy->SetBorderSize(0);
  Leg_energy->SetFillColorAlpha(0, 0.01);
  Leg_energy->Draw("same");

  TPaveText *Text_energy = new TPaveText(0.5, 0.45, 0.79, 0.65, "NDC");
  Text_energy->SetTextFont(42);
  Text_energy->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
  Text_energy->AddText("1 < Q^{2} < 20 GeV^{2}/c^{2}");
  Text_energy->AddText("0.01 < y < 0.95");
  //cent_text->AddText("#Lambda^{0} and #bar{#Lambda^{0}}");
  Text_energy->SetFillColorAlpha(0, 0.01);
  Text_energy->Draw("same");

  energy_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_energy.png", e_energy, p_energy));


  TCanvas *energy_zoom_can = new TCanvas("energy_zoom_can", "energy_zoom_can", 1200, 1000);

  energy_zoom_can->cd();

  h_energy_zoom_MC->GetXaxis()->SetTitle("E_{e} (GeV)");
  h_energy_zoom_MC->GetXaxis()->CenterTitle();
  h_energy_zoom_MC->GetYaxis()->SetTitle("Counts");
  h_energy_zoom_MC->GetYaxis()->CenterTitle();
  h_energy_zoom_MC->SetMarkerStyle(20);
  h_energy_zoom_MC->SetMarkerColor(kRed);
  h_energy_zoom_MC->SetLineColor(kRed);
  h_energy_zoom_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  energy_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_zoom_energy.png", e_energy, p_energy));


  TCanvas *momentum_can = new TCanvas("momentum_can", "momentum_can", 1200, 1000);

  momentum_can->cd();

  h_momentum_MC->GetXaxis()->SetTitle("p (GeV/c)");
  h_momentum_MC->GetXaxis()->CenterTitle();
  h_momentum_MC->GetYaxis()->SetTitle("Counts");
  h_momentum_MC->GetYaxis()->CenterTitle();
  h_momentum_MC->SetMarkerStyle(20);
  h_momentum_MC->SetMarkerColor(kRed);
  h_momentum_MC->SetLineColor(kRed);
  h_momentum_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  momentum_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_momentum.png", e_energy, p_energy));



  TCanvas *Q2_can = new TCanvas("Q2_can", "Q2_can", 1200, 1000);

  Q2_can->cd();

  h_Q2_MC->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_Q2_MC->GetXaxis()->CenterTitle();
  h_Q2_MC->GetYaxis()->SetTitle("Counts");
  h_Q2_MC->GetYaxis()->CenterTitle();
  h_Q2_MC->SetMarkerStyle(20);
  h_Q2_MC->SetMarkerColor(kRed);
  h_Q2_MC->SetLineColor(kRed);
  h_Q2_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  Q2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_Q2.png", e_energy, p_energy));


  TCanvas *Q2_zoom_can = new TCanvas("Q2_zoom_can", "Q2_zoom_can", 1200, 1000);

  Q2_zoom_can->cd();

  h_Q2_zoom_MC->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_Q2_zoom_MC->GetXaxis()->CenterTitle();
  h_Q2_zoom_MC->GetYaxis()->SetTitle("Counts");
  h_Q2_zoom_MC->GetYaxis()->CenterTitle();
  h_Q2_zoom_MC->SetMarkerStyle(20);
  h_Q2_zoom_MC->SetMarkerColor(kRed);
  h_Q2_zoom_MC->SetLineColor(kRed);
  h_Q2_zoom_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  Q2_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_Q2_zoom.png", e_energy, p_energy));



  TCanvas *y_inelPar_can = new TCanvas("y_inelPar_can", "y_inelPar_can", 1200, 1000);

  y_inelPar_can->cd();

  h_y_inelPar_MC->GetXaxis()->SetTitle("y");
  h_y_inelPar_MC->GetXaxis()->CenterTitle();
  h_y_inelPar_MC->GetYaxis()->SetTitle("Counts");
  h_y_inelPar_MC->GetYaxis()->CenterTitle();
  h_y_inelPar_MC->SetMarkerStyle(20);
  h_y_inelPar_MC->SetMarkerColor(kRed);
  h_y_inelPar_MC->SetLineColor(kRed);
  h_y_inelPar_MC->SetMinimum(0);
  h_y_inelPar_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  y_inelPar_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_y_inelPar.png", e_energy, p_energy));


  TCanvas *y_inelPar_zoom_can = new TCanvas("y_inelPar_zoom_can", "y_inelPar_zoom_can", 1200, 1000);

  y_inelPar_zoom_can->cd();

  h_y_inelPar_zoom_MC->GetXaxis()->SetTitle("y");
  h_y_inelPar_zoom_MC->GetXaxis()->CenterTitle();
  h_y_inelPar_zoom_MC->GetYaxis()->SetTitle("Counts");
  h_y_inelPar_zoom_MC->GetYaxis()->CenterTitle();
  h_y_inelPar_zoom_MC->SetMarkerStyle(20);
  h_y_inelPar_zoom_MC->SetMarkerColor(kRed);
  h_y_inelPar_zoom_MC->SetLineColor(kRed);
  h_y_inelPar_zoom_MC->SetMinimum(0);;
  h_y_inelPar_zoom_MC->Draw("p e");

  Leg_energy->Draw("same");

  Text_energy->Draw("same");

  y_inelPar_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Electron_y_inelPar_zoom.png", e_energy, p_energy));

  //gPad->SetLogy();

  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {

        TCanvas *eta_scat_ele_can = new TCanvas(Form("eta_scat_ele_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), Form("eta_scat_ele_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), 2000, 1500);

        eta_scat_ele_can->cd();

        gPad->SetLogy();

        h_eta_scat_ele[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        if(h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Integral() == 0 ) continue;

        h_eta_ele[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_ele_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_pi_plus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_plus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_K_plus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_K_plus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_proton[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_proton_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        h_eta_positron[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_positron_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_pi_minus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_K_minus[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));
        h_eta_anti_proton[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));


        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->CenterTitle();
        //for lin scale
        //if( h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum() < h_eta_pi_plus[mom_bin][Q2bin][y_bin]->GetMaximum() ) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0., h_eta_pi_plus[mom_bin][Q2bin][y_bin]->GetMaximum()*1.7);
        //else h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0.5, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum()*1.7);
        //for log scale
        if( h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum() < h_eta_pi_plus[mom_bin][Q2bin][y_bin]->GetMaximum() ) h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0.5, h_eta_pi_plus[mom_bin][Q2bin][y_bin]->GetMaximum()*100);
        else h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetYaxis()->SetRangeUser(0.5, h_eta_scat_ele[mom_bin][Q2bin][y_bin]->GetMaximum()*200);
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerStyle(20);
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMarkerColor(kRed);
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetLineColor(kRed);
        //h_eta_scat_ele[mom_bin][Q2bin][y_bin]->SetMinimum(0);
        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("p e");


        h_eta_positron[mom_bin][Q2bin][y_bin]->SetMarkerStyle(24);
        h_eta_positron[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_positron[mom_bin][Q2bin][y_bin]->SetMarkerColor(kBlue);
        h_eta_positron[mom_bin][Q2bin][y_bin]->SetLineColor(kBlue);
        h_eta_positron[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_pi_plus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(25);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin]->SetLineColor(1);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_K_plus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(28);
        h_eta_K_plus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_K_plus[mom_bin][Q2bin][y_bin]->SetMarkerColor(6);
        h_eta_K_plus[mom_bin][Q2bin][y_bin]->SetLineColor(6);
        h_eta_K_plus[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_proton[mom_bin][Q2bin][y_bin]->SetMarkerStyle(27);
        h_eta_proton[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_proton[mom_bin][Q2bin][y_bin]->SetMarkerColor(kOrange+2);
        h_eta_proton[mom_bin][Q2bin][y_bin]->SetLineColor(kOrange+2);
        h_eta_proton[mom_bin][Q2bin][y_bin]->Draw("p e same");

        TLegend *Leg_bins = new TLegend(0.2, 0.63, 0.54, 0.85 );
        Leg_bins->SetTextFont(42);
        Leg_bins->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "Scattered e^{-}");
        Leg_bins->AddEntry(h_eta_positron[mom_bin][Q2bin][y_bin], "e^{+}");
        Leg_bins->AddEntry(h_eta_pi_plus[mom_bin][Q2bin][y_bin], "#pi^{+}");
        Leg_bins->AddEntry(h_eta_K_plus[mom_bin][Q2bin][y_bin], "K^{+}");
        Leg_bins->AddEntry(h_eta_proton[mom_bin][Q2bin][y_bin], "p");
        Leg_bins->SetBorderSize(0);
        Leg_bins->SetFillColorAlpha(0, 0.01);
        Leg_bins->Draw("same");

        TPaveText *Text_bins = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
        Text_bins->SetTextFont(42);
        Text_bins->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        if( mom_bin < nMomBins )
        {
          Text_bins->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
          Text_bins->AddText(Form("%0.1f < Q^{2} < %0.1f  GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
          Text_bins->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        }

        Text_bins->SetFillColorAlpha(0, 0.01);
        Text_bins->Draw("same");

        eta_scat_ele_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Eta_scat_e_vs_pos_part_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy,  mom_bin , Q2bin, y_bin));

        //h_eta_proton[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_eta_proton_mom_%i_Q2_%i_y_%i", mom_bin , Q2bin, y_bin));

        //__________________________________________________________________________________________________________________________________________


        TCanvas *eta_scat_ele_anti_can = new TCanvas(Form("eta_scat_ele_anti_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), Form("eta_scat_ele_anti_can_mom_%i_Q2_%i_yBin_%i", mom_bin , Q2bin, y_bin), 2000, 1500);

        eta_scat_ele_anti_can->cd();

        gPad->SetLogy();

        h_eta_scat_ele[mom_bin][Q2bin][y_bin]->Draw("p e");


        h_eta_ele[mom_bin][Q2bin][y_bin]->SetMarkerStyle(24);
        h_eta_ele[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_ele[mom_bin][Q2bin][y_bin]->SetMarkerColor(kBlue);
        h_eta_ele[mom_bin][Q2bin][y_bin]->SetLineColor(kBlue);
        h_eta_ele[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(25);
        h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetMarkerColor(1);
        h_eta_pi_minus[mom_bin][Q2bin][y_bin]->SetLineColor(1);
        h_eta_pi_minus[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerStyle(28);
        h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetMarkerColor(6);
        h_eta_K_minus[mom_bin][Q2bin][y_bin]->SetLineColor(6);
        h_eta_K_minus[mom_bin][Q2bin][y_bin]->Draw("p e same");


        h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerStyle(27);
        h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerSize(2);
        h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetMarkerColor(kOrange+2);
        h_eta_anti_proton[mom_bin][Q2bin][y_bin]->SetLineColor(kOrange+2);
        h_eta_anti_proton[mom_bin][Q2bin][y_bin]->Draw("p e same");

        TLegend *Leg_bins_anti = new TLegend(0.2, 0.63, 0.54, 0.85 );
        Leg_bins_anti->SetTextFont(42);
        Leg_bins_anti->AddEntry(h_eta_scat_ele[mom_bin][Q2bin][y_bin], "Scattered e^{-}");
        Leg_bins_anti->AddEntry(h_eta_ele[mom_bin][Q2bin][y_bin], "e^{-}");
        Leg_bins_anti->AddEntry(h_eta_pi_plus[mom_bin][Q2bin][y_bin], "#pi^{-}");
        Leg_bins_anti->AddEntry(h_eta_K_plus[mom_bin][Q2bin][y_bin], "K^{-}");
        Leg_bins_anti->AddEntry(h_eta_anti_proton[mom_bin][Q2bin][y_bin], "#bar{p}");
        Leg_bins_anti->SetBorderSize(0);
        Leg_bins_anti->SetFillColorAlpha(0, 0.01);
        Leg_bins_anti->Draw("same");

        Text_bins->Draw("same");

        eta_scat_ele_anti_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/Eta_scat_e_vs_neg_part_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));

      }

    }
  }



  //fitPiTOF->Write();
  //fitKTOF->Write();

  inFile->Close();

  return;

}
