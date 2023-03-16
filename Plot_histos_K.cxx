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

void Plot_histos_K(int e_energy = 18, int p_energy = 275)
{
  if( !(e_energy == 18 && p_energy == 275) && !(e_energy == 10 && p_energy == 100) && !(e_energy == 5 && p_energy == 41))
  {
    cout<<"Invalid beam energies."<<endl;

    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineScalePS(1);//for pdf lines

  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 24;
  float const mom_bins[nMomBins+1] = {0.,1.0,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6,6.5,7,7.5,8,10,12};

  //pfRICH eta acceptance mc_mom.Eta() > -3.8 && mc_mom.Eta() < -1.5
  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -3.8, -3, -2.5, -2, -1.5};

  //load all files
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output_K.root", e_energy, p_energy, e_energy, p_energy), "READ");

  TH1F* h_hadron_zh = (TH1F*)inFile->Get("h_hadron_zh");

  //MC histograms
  //p bins + p integrated
  //TH1F *h_eta_K_plus_lead_MC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_lead_MC[nQ2bins][nyInelParBins];

  //TH1F *h_p_K_plus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_K_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_pi_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_p_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];

  TH1F *h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[nQ2bins][nyInelParBins];
  TH1F *h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins]; //wider bins for efficiency

  TH1F *h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];

  TH1F *h_p_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins]; //wider bins for efficiency


  //MC -> RC histograms
  TH1F *h_eta_K_minus_lead_RC[nQ2bins][nyInelParBins];

  //TH1F *h_p_K_plus_lead_RC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_K_minus_lead_RC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_pi_minus_lead_RC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_p_minus_lead_RC_Q2_y[nQ2bins][nyInelParBins];

  TH1F *h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[nQ2bins][nyInelParBins];
  TH1F *h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins]; //wider bins for efficiency

  TH1F *h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];

  TH1F *h_p_K_minus_lead_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_K_minus_lead_RC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins]; //wider bins for efficiency


  //____________________________________________________________

  //K-
  TH1F *h_K_minus_lead_purity_pfRICH_MC_tot[nQ2bins][nyInelParBins];

  //_______________________________________________________________

  TCanvas *z_h_can = new TCanvas("z_h_can", "z_h_can", 1200, 1000);
  z_h_can->cd();

  h_hadron_zh->Draw();

  z_h_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/h_hadron_zh.png",e_energy, p_energy));


  TCanvas *K_minus_lead_eta_pfRICH_one_can = new TCanvas("K_minus_lead_eta_pfRICH_one_can", "K_minus_lead_eta_pfRICH_one_can", 8000, 6000);
  K_minus_lead_eta_pfRICH_one_can->Divide(4,4);

  TCanvas *Lead_p_pfRICH_one_can = new TCanvas("Lead_p_pfRICH_one_can", "Lead_p_pfRICH_one_can", 8000, 6000);
  Lead_p_pfRICH_one_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_one_can = new TCanvas("K_minus_lead_p_pfRICH_one_can", "K_minus_lead_p_pfRICH_one_can", 8000, 6000);
  K_minus_lead_p_pfRICH_one_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_compare_one_can = new TCanvas("K_minus_lead_p_pfRICH_compare_one_can", "K_minus_lead_p_pfRICH_compare_one_can", 8000, 6000);
  K_minus_lead_p_pfRICH_compare_one_can->Divide(4,4);


  TCanvas *K_minus_lead_p_pfRICH_purity_one_can = new TCanvas("K_minus_lead_p_pfRICH_purity_one_can", "K_minus_lead_p_pfRICH_purity_one_can", 8000, 6000);
  K_minus_lead_p_pfRICH_purity_one_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_purity_wide_one_can = new TCanvas("K_minus_lead_p_pfRICH_purity_wide_one_can", "K_minus_lead_p_pfRICH_purity_wide_one_can", 8000, 6000);
  K_minus_lead_p_pfRICH_purity_wide_one_can->Divide(4,4);

  //MC->RC canvases
  TCanvas *K_minus_lead_eta_pfRICH_one_RC_can = new TCanvas("K_minus_lead_eta_pfRICH_one_RC_can", "K_minus_lead_eta_pfRICH_one_RC_can", 8000, 6000);
  K_minus_lead_eta_pfRICH_one_RC_can->Divide(4,4);

  TCanvas *Lead_p_pfRICH_one_RC_can = new TCanvas("Lead_p_pfRICH_one_RC_can", "Lead_p_pfRICH_one_RC_can", 8000, 6000);
  Lead_p_pfRICH_one_RC_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_one_RC_can = new TCanvas("K_minus_lead_p_pfRICH_one_RC_can", "K_minus_lead_p_pfRICH_one_RC_can", 8000, 6000);
  K_minus_lead_p_pfRICH_one_RC_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_compare_one_RC_can = new TCanvas("K_minus_lead_p_pfRICH_compare_one_RC_can", "K_minus_lead_p_pfRICH_compare_one_RC_can", 8000, 6000);
  K_minus_lead_p_pfRICH_compare_one_RC_can->Divide(4,4);


  TCanvas *K_minus_lead_p_pfRICH_purity_one_RC_can = new TCanvas("K_minus_lead_p_pfRICH_purity_one_RC_can", "K_minus_lead_p_pfRICH_purity_one_RC_can", 8000, 6000);
  K_minus_lead_p_pfRICH_purity_one_RC_can->Divide(4,4);

  TCanvas *K_minus_lead_p_pfRICH_purity_wide_one_RC_can = new TCanvas("K_minus_lead_p_pfRICH_purity_wide_one_RC_can", "K_minus_lead_p_pfRICH_purity_wide_one_RC_can", 8000, 6000);
  K_minus_lead_p_pfRICH_purity_wide_one_RC_can->Divide(4,4);


  for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
  {
    for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
    {

      TPaveText *Text_K_MC = new TPaveText(0.15, 0.65, 0.4, 0.85, "NDC");
      Text_K_MC->SetTextFont(42);
      Text_K_MC->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
      //Text_K_MC->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
      Text_K_MC->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
      Text_K_MC->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
      Text_K_MC->AddText("z_{h} > 0.2");
      Text_K_MC->SetFillColorAlpha(0, 0.01);

/*
      //K purity histograms
      //K-
      h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_K_minus_lead_purity_pfRICH_MC_tot_Q2_%i_y_%i" , Q2bin, y_bin));

      //pT and eta (witihn pfRICH acceptance) integrated purity
      double K_minus_lead_purity = 0;

      if( h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin]->Integral() > 0 )
      {
        TCanvas *K_minus_lead_purity_pfRICH_RC_can = new TCanvas(Form("K_minus_lead_purity_pfRICH_RC_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_minus_lead_purity_pfRICH_RC_can_Q2_%i_y_%i" , Q2bin, y_bin), 2000, 1500);

        K_minus_lead_purity_pfRICH_RC_can->cd();


        h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin]->Draw();

        double nK = h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin]->GetBinContent(2);
        double nBackground = h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin]->GetBinContent(1);

        K_minus_lead_purity = ( nK/(nK+nBackground) )*100;


        TPaveText *Text_K_minus_lead_purity = new TPaveText(0.55, 0.6, 0.84, 0.85, "NDC");
        Text_K_minus_lead_purity->SetTextFont(42);
        Text_K_minus_lead_purity->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        Text_K_minus_lead_purity->AddText(Form("K^{-} pfRICH purity: %0.1f%%", K_minus_lead_purity));
        //Text_K_minus_lead_purity->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_minus_lead_purity->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_minus_lead_purity->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_minus_lead_purity->AddText("z_{h} > 0.2");
        Text_K_minus_lead_purity->SetFillColorAlpha(0, 0.01);
        Text_K_minus_lead_purity->Draw("same");

        //if( mom_bin == nMomBins )
        //{
          //K_minus_lead_purity_pfRICH_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_purity_pfRICH_RC_p_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin, Q2bin, y_bin));
          //K_minus_lead_purity_pfRICH_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_purity_pfRICH_RC_p_%i_Q2_%i_yBin_%i.pdf", e_energy, p_energy, mom_bin, Q2bin, y_bin));
        //}

      }
      //________________________________________________________________________________________________________________________________________

*/

      //K MC histos
      //K-
      h_eta_K_minus_lead_MC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_pi_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_p_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_MC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_as_K_minus_lead_MC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_pi_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_p_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));


      if( h_eta_K_minus_lead_MC[Q2bin][y_bin]->Integral() > 0 )
      {

        K_minus_lead_eta_pfRICH_one_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_eta_K_minus_lead_MC[Q2bin][y_bin]->GetMaximum()*100);
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->SetMarkerStyle(34);
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->SetMarkerColor(1);
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->SetLineColor(1);
        h_eta_K_minus_lead_MC[Q2bin][y_bin]->Draw("p e");


        Text_K_MC->Draw("same");

        TLegend *Leg_bins_K = new TLegend(0.41, 0.7, 0.89, 0.85);
        Leg_bins_K->SetTextFont(42);
        Leg_bins_K->AddEntry(h_eta_K_minus_lead_MC[Q2bin][y_bin], "MC K^{-}");
        Leg_bins_K->SetBorderSize(0);
        Leg_bins_K->SetFillColorAlpha(0, 0.01);
        Leg_bins_K->Draw("same");

        //K_minus_lead_eta_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_eta_p_%i_Q2_%i_y_%i.png" ,e_energy, p_energy, mom_bin, Q2bin, y_bin));
        //K_minus_lead_eta_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_eta_p_%i_Q2_%i_y_%i.pdf" ,e_energy, p_energy, mom_bin, Q2bin, y_bin));

        //_____________________________________________________________________________

        Lead_p_pfRICH_one_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetLineColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->Draw("p e");


        h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(24);
        h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerColor(2);
        h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetLineColor(2);
        h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin]->Draw("p e same");

        h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(25);
        h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerColor(4);
        h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetLineColor(4);
        h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin]->Draw("p e same");

        TPaveText *Text_K_MC_p = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
        Text_K_MC_p->SetTextFont(42);
        Text_K_MC_p->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        //Text_K_MC_p->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_MC_p->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_MC_p->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_MC_p->AddText("-3.8 < #eta < -1.5");
        Text_K_MC_p->AddText("z_{h} > 0.2");
        Text_K_MC_p->SetFillColorAlpha(0, 0.01);
        Text_K_MC_p->Draw("same");

        TLegend *Leg_bins_lead_p = new TLegend(0.41, 0.7, 0.89, 0.85);
        Leg_bins_lead_p->SetTextFont(42);
        Leg_bins_lead_p->AddEntry(h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin], "MC K^{-}");
        Leg_bins_lead_p->AddEntry(h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin], "MC #pi^{-}");
        Leg_bins_lead_p->AddEntry(h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin], "MC #bar{p}");
        Leg_bins_lead_p->SetBorderSize(0);
        Leg_bins_lead_p->SetFillColorAlpha(0, 0.01);
        Leg_bins_lead_p->Draw("same");


        //TCanvas *K_minus_lead_p_pfRICH_can = new TCanvas(Form("K_minus_lead_p_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_minus_lead_p_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), 2000, 1500);
        //K_minus_lead_p_pfRICH_can->cd();

        K_minus_lead_p_pfRICH_one_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();


        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetLineColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->Draw("p e");


        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");




        Text_K_MC_p->Draw("same");

        TLegend *Leg_bins_K_p = new TLegend(0.41, 0.7, 0.89, 0.85);
        Leg_bins_K_p->SetTextFont(42);
        Leg_bins_K_p->AddEntry(h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin], "MC K^{-}");
        Leg_bins_K_p->AddEntry(h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_p->SetBorderSize(0);
        Leg_bins_K_p->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_p->Draw("same");




        //gPad->SetLogy();

        K_minus_lead_p_pfRICH_compare_one_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e");


        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(28);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(6);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(6);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");

        h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(24);
        h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(2);
        h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(2);
        h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");

        h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(25);
        h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(4);
        h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(4);
        h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");


        Text_K_MC_p->Draw("same");

        TLegend *Leg_bins_K_p_compare = new TLegend(0.41, 0.6, 0.89, 0.85);
        Leg_bins_K_p_compare->SetTextFont(42);
        Leg_bins_K_p_compare->AddEntry(h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_p_compare->AddEntry(h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} as K^{-}");
        Leg_bins_K_p_compare->AddEntry(h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH #pi^{-} as K^{-}");
        Leg_bins_K_p_compare->AddEntry(h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH #bar{p} as K^{-}");
        Leg_bins_K_p_compare->SetBorderSize(0);
        Leg_bins_K_p_compare->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_p_compare->Draw("same");




        //____________________________________________________________________________________________________________________________

        //add K+ here if needed

        //purity as a function of p
        K_minus_lead_p_pfRICH_purity_one_can->cd((4*Q2bin+y_bin)+1);


        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin] = (TH1F*)h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Clone(Form("h_p_K_as_K_minus_lead_MC_pfRICH_clone_Q2_%i_y_%i" , Q2bin, y_bin));
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Divide(h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin], h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin], 1,1,"b");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->SetTitle("Purity");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->SetRangeUser(0,1.9);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Smooth(5);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Draw("p e");



        TPaveText *Text_K_MC_p_purity = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
        Text_K_MC_p_purity->SetTextFont(42);
        Text_K_MC_p_purity->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        //Text_K_MC_p_purity->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_MC_p_purity->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_MC_p_purity->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_MC_p_purity->AddText("-3.8 < #eta < -1.5");
        Text_K_MC_p_purity->AddText("z_{h} > 0.2");
        Text_K_MC_p_purity->SetFillColorAlpha(0, 0.01);
        Text_K_MC_p_purity->Draw("same");

        TLegend *Leg_bins_K_purity = new TLegend(0.41, 0.6, 0.89, 0.85);
        Leg_bins_K_purity->SetTextFont(42);
        Leg_bins_K_purity->AddEntry(h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_purity->SetBorderSize(0);
        Leg_bins_K_purity->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_purity->Draw("same");

        //______________________________________________________________________________


        K_minus_lead_p_pfRICH_purity_wide_one_can->cd((4*Q2bin+y_bin)+1);


        //h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_clone[Q2bin][y_bin] = (TH1F*)h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin]->Clone(Form("h_p_K_as_K_minus_lead_MC_pfRICH_clone_Q2_%i_y_%i" , Q2bin, y_bin));
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Divide(h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin], h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin], 1,1,"b");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->SetTitle("Purity");
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->SetRangeUser(0,1.9);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetLineColor(8);
        //h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Smooth(2);
        h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Draw("p e");

        Text_K_MC_p_purity->Draw("same");

        Leg_bins_K_purity->Draw("same");

        //______________________________________________________________________________

      }//end if itegral MC

      //_________________________________________________________________________________________________________________________

      //MC->RC
      //MC->RC histos
      h_eta_K_minus_lead_RC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_minus_lead_RC_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_RC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_pi_minus_lead_RC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_p_minus_lead_RC_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_lead_RC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin));

      h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_as_K_minus_lead_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_K_as_K_minus_lead_RC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_pi_as_K_minus_lead_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));
      h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_p_p_as_K_minus_lead_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));


      //use MC or RC in if?
      //if( h_eta_K_minus_lead_RC[Q2bin][y_bin]->Integral() > 0 )
      if( h_eta_K_minus_lead_MC[Q2bin][y_bin]->Integral() > 0 )
      {
/*
        Lead_p_pfRICH_one_RC_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerColor(1);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetLineColor(1);
        h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->Draw("p e");


        h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(24);
        h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerColor(2);
        h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetLineColor(2);
        h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin]->Draw("p e same");

        h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(25);
        h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerColor(4);
        h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin]->SetLineColor(4);
        h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin]->Draw("p e same");

        TPaveText *Text_K_RC_p = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
        Text_K_RC_p->SetTextFont(42);
        Text_K_RC_p->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        //Text_K_RC_p->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_RC_p->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_RC_p->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_RC_p->AddText("-3.8 < #eta < -1.5");
        Text_K_RC_p->AddText("z_{h} > 0.2");
        Text_K_RC_p->SetFillColorAlpha(0, 0.01);
        Text_K_RC_p->Draw("same");

        TLegend *Leg_bins_lead_p = new TLegend(0.41, 0.7, 0.89, 0.85);
        Leg_bins_lead_p->SetTextFont(42);
        Leg_bins_lead_p->AddEntry(h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin], "MC K^{-}");
        Leg_bins_lead_p->AddEntry(h_p_pi_minus_lead_RC_Q2_y[Q2bin][y_bin], "MC #pi^{-}");
        Leg_bins_lead_p->AddEntry(h_p_p_minus_lead_RC_Q2_y[Q2bin][y_bin], "MC #bar{p}");
        Leg_bins_lead_p->SetBorderSize(0);
        Leg_bins_lead_p->SetFillColorAlpha(0, 0.01);
        Leg_bins_lead_p->Draw("same");


        //TCanvas *K_minus_lead_p_pfRICH_RC_can = new TCanvas(Form("K_minus_lead_p_pfRICH_RC_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_minus_lead_p_pfRICH_RC_can_Q2_%i_y_%i" , Q2bin, y_bin), 2000, 1500);
        //K_minus_lead_p_pfRICH_RC_can->cd();
*/
        K_minus_lead_p_pfRICH_one_RC_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->SetLineColor(1);
        h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin]->Draw("p e");


        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");


        TPaveText *Text_K_RC_p = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
        Text_K_RC_p->SetTextFont(42);
        Text_K_RC_p->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        //Text_K_RC_p->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_RC_p->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_RC_p->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_RC_p->AddText("-3.8 < #eta < -1.5");
        Text_K_RC_p->AddText("z_{h} > 0.2");
        Text_K_RC_p->SetFillColorAlpha(0, 0.01);
        Text_K_RC_p->Draw("same");
        Text_K_RC_p->Draw("same");

        TLegend *Leg_bins_K_p = new TLegend(0.41, 0.7, 0.89, 0.85);
        Leg_bins_K_p->SetTextFont(42);
        Leg_bins_K_p->AddEntry(h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin], "MC K^{-}");
        Leg_bins_K_p->AddEntry(h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_p->SetBorderSize(0);
        Leg_bins_K_p->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_p->Draw("same");




        //gPad->SetLogy();

        K_minus_lead_p_pfRICH_compare_one_RC_can->cd( (4*Q2bin+y_bin)+1 );

        gPad->SetLogy();

        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->GetYaxis()->SetRangeUser(1e-1, h_p_K_minus_lead_RC_Q2_y[Q2bin][y_bin]->GetMaximum()*100);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e");


        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(28);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(6);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(6);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");

        h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(24);
        h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(2);
        h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(2);
        h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");

        h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerStyle(25);
        h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerColor(4);
        h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->SetLineColor(4);
        h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Draw("p e same");


        Text_K_RC_p->Draw("same");

        TLegend *Leg_bins_K_p_compare = new TLegend(0.41, 0.6, 0.89, 0.85);
        Leg_bins_K_p_compare->SetTextFont(42);
        Leg_bins_K_p_compare->AddEntry(h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_p_compare->AddEntry(h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH K^{-} as K^{-}");
        Leg_bins_K_p_compare->AddEntry(h_p_pi_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH #pi^{-} as K^{-}");
        Leg_bins_K_p_compare->AddEntry(h_p_p_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], "pfRICH #bar{p} as K^{-}");
        Leg_bins_K_p_compare->SetBorderSize(0);
        Leg_bins_K_p_compare->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_p_compare->Draw("same");




        //____________________________________________________________________________________________________________________________

        //add K+ here if needed

        //purity as a function of p
        K_minus_lead_p_pfRICH_purity_one_RC_can->cd((4*Q2bin+y_bin)+1);


        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin] = (TH1F*)h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Clone(Form("h_p_K_as_K_minus_lead_RC_pfRICH_clone_Q2_%i_y_%i" , Q2bin, y_bin));
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Divide(h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin], h_p_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin], 1,1,"b");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->SetTitle("Purity");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->GetYaxis()->SetRangeUser(0,1.9);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->SetLineColor(8);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Smooth(5);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin]->Draw("p e");



        TPaveText *Text_K_RC_p_purity = new TPaveText(0.15, 0.6, 0.4, 0.85, "NDC");
        Text_K_RC_p_purity->SetTextFont(42);
        Text_K_RC_p_purity->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        //Text_K_RC_p_purity->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_RC_p_purity->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_RC_p_purity->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_RC_p_purity->AddText("-3.8 < #eta < -1.5");
        Text_K_RC_p_purity->AddText("z_{h} > 0.2");
        Text_K_RC_p_purity->SetFillColorAlpha(0, 0.01);
        Text_K_RC_p_purity->Draw("same");

        TLegend *Leg_bins_K_purity = new TLegend(0.41, 0.6, 0.89, 0.85);
        Leg_bins_K_purity->SetTextFont(42);
        Leg_bins_K_purity->AddEntry(h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin], "pfRICH K^{-} candidates");
        Leg_bins_K_purity->SetBorderSize(0);
        Leg_bins_K_purity->SetFillColorAlpha(0, 0.01);
        Leg_bins_K_purity->Draw("same");

        //______________________________________________________________________________


        K_minus_lead_p_pfRICH_purity_wide_one_RC_can->cd((4*Q2bin+y_bin)+1);


        //h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_clone[Q2bin][y_bin] = (TH1F*)h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH[Q2bin][y_bin]->Clone(Form("h_p_K_as_K_minus_lead_RC_pfRICH_clone_Q2_%i_y_%i" , Q2bin, y_bin));
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Divide(h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin], h_p_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin], 1,1,"b");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetXaxis()->SetRangeUser(0,12);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->SetTitle("Purity");
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->GetYaxis()->SetRangeUser(0,1.9);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerStyle(34);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerColor(8);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetMarkerSize(3.5);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->SetLineColor(8);
        //h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Smooth(2);
        h_p_K_as_K_minus_lead_RC_Q2_y_pfRICH_eff[Q2bin][y_bin]->Draw("p e");

        Text_K_RC_p_purity->Draw("same");

        Leg_bins_K_purity->Draw("same");

        //______________________________________________________________________________

      }


    }//end for y bins
  }//end for Q2 bins


  //MC plots
  K_minus_lead_eta_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_eta_one.png" ,e_energy, p_energy));
  K_minus_lead_eta_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_eta_one.pdf" ,e_energy, p_energy));


  Lead_p_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/Lead_p_one.png" ,e_energy, p_energy));
  Lead_p_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/Lead_p_one.pdf" ,e_energy, p_energy));


  K_minus_lead_p_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one.pdf" ,e_energy, p_energy));

  K_minus_lead_p_pfRICH_compare_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_compare.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_compare_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_compare.pdf" ,e_energy, p_energy));

  K_minus_lead_p_pfRICH_purity_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_one.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_purity_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_one.pdf" ,e_energy, p_energy));


  K_minus_lead_p_pfRICH_purity_wide_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_wide_bins_one.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_purity_wide_one_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_wide_bins_one.pdf" ,e_energy, p_energy));


  //MC->RC plots
  K_minus_lead_p_pfRICH_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_RC.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_RC.pdf" ,e_energy, p_energy));

  K_minus_lead_p_pfRICH_compare_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_RC_compare.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_compare_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_one_RC_compare.pdf" ,e_energy, p_energy));

  K_minus_lead_p_pfRICH_purity_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_one_RC.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_purity_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_one_RC.pdf" ,e_energy, p_energy));

  K_minus_lead_p_pfRICH_purity_wide_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_wide_bins_one_RC.png" ,e_energy, p_energy));
  K_minus_lead_p_pfRICH_purity_wide_one_RC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_lead_p_purity_wide_bins_one_RC.pdf" ,e_energy, p_energy));



  //______________________________________________________________________________________________________________________________________________________________


  inFile->Close();

  return;

}
