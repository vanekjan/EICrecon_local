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
void Plot_histos_K(int e_energy = 18, int p_energy = 275)
{
  if( !(e_energy == 18 && p_energy == 275) && !(e_energy == 10 && p_energy == 100) && !(e_energy == 5 && p_energy == 41))
  {
    cout<<"Invalid beam energies."<<endl;

    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 11;
  float const mom_bins[nMomBins+1] = { 0,0.5,1,1.5,2,3,4,5,6,7,10, 18 };


  //pfRICH eta acceptance mc_mom.Eta() > -3.8 && mc_mom.Eta() < -1.5
  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -3.8, -3, -2.5, -2, -1.5};

  //load all files
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output_K.root", e_energy, p_energy, e_energy, p_energy), "READ");

  //K spectra
  TH1F *h_eta_K_plus_MC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC[nQ2bins][nyInelParBins];

  TH1F *h_p_K_plus_MC[nEtaBins+1];
  TH1F *h_p_K_minus_MC[nEtaBins+1];


  TH1F *h_eta_K_plus_MC_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_pfRICH[nQ2bins][nyInelParBins];

  TH1F *h_p_K_plus_MC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_MC_pfRICH[nEtaBins+1];


  //K purity

  TH1F *h_K_minus_purity_pfRICH_MC[nQ2bins][nyInelParBins];

  TH1F *h_K_minus_purity_p_eta_pfRICH_MC[nEtaBins+1];


  TH1F *h_K_plus_purity_pfRICH_MC[nQ2bins][nyInelParBins];

  TH1F *h_K_plus_purity_p_eta_pfRICH_MC[nEtaBins+1];



  for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
  {
    for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
    {

      TPaveText *Text_K_MC = new TPaveText(0.25, 0.65, 0.5, 0.85, "NDC");
      Text_K_MC->SetTextFont(42);
      Text_K_MC->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
      //Text_K_MC->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
      Text_K_MC->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
      Text_K_MC->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
      Text_K_MC->SetFillColorAlpha(0, 0.01);

      //K spectra

      //K-
      h_eta_K_minus_MC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_minus_MC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_eta_K_minus_MC_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_minus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));

      if( h_eta_K_minus_MC[Q2bin][y_bin]->Integral() > 0 )
      {
        TCanvas *K_minus_eta_pfRICH_can = new TCanvas(Form("K_minus_eta_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_minus_eta_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), 1200, 1000);
        K_minus_eta_pfRICH_can->cd();

        gPad->SetLogy();


        h_eta_K_minus_MC[Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_K_minus_MC[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_K_minus_MC[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_K_minus_MC[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_eta_K_minus_MC[Q2bin][y_bin]->SetMarkerStyle(28);
        h_eta_K_minus_MC[Q2bin][y_bin]->SetMarkerColor(4);
        h_eta_K_minus_MC[Q2bin][y_bin]->SetLineColor(4);
        h_eta_K_minus_MC[Q2bin][y_bin]->Draw("p e");


        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_eta_K_minus_MC_pfRICH[Q2bin][y_bin]->Draw("p e same");

        Text_K_MC->Draw("same");

        K_minus_eta_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_eta_Q2_%i_y_%i.png" ,e_energy, p_energy, Q2bin, y_bin));

      }

      //K+
      h_eta_K_plus_MC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_plus_MC_Q2_%i_y_%i" , Q2bin, y_bin));
      h_eta_K_plus_MC_pfRICH[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_eta_K_plus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin));

      if( h_eta_K_plus_MC[Q2bin][y_bin]->Integral() > 0 )
      {

        TCanvas *K_plus_eta_pfRICH_can = new TCanvas(Form("K_plus_eta_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_plus_eta_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), 1200, 1000);
        K_plus_eta_pfRICH_can->cd();

        gPad->SetLogy();

        h_eta_K_plus_MC[Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_K_plus_MC[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_K_plus_MC[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_K_plus_MC[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_eta_K_plus_MC[Q2bin][y_bin]->SetMarkerStyle(28);
        h_eta_K_plus_MC[Q2bin][y_bin]->SetMarkerColor(4);
        h_eta_K_plus_MC[Q2bin][y_bin]->SetLineColor(4);
        h_eta_K_plus_MC[Q2bin][y_bin]->Draw("p e");


        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->GetXaxis()->SetTitle("#eta");
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->GetXaxis()->CenterTitle();
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->GetYaxis()->SetTitle("Counts");
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->GetYaxis()->CenterTitle();
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->SetMarkerStyle(34);
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->SetMarkerColor(8);
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->SetLineColor(8);
        h_eta_K_plus_MC_pfRICH[Q2bin][y_bin]->Draw("p e same");

        Text_K_MC->Draw("same");

        K_plus_eta_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_plus_eta__Q2_%i_y_%i.png" ,e_energy, p_energy, Q2bin, y_bin));

      }


      //K purity histograms

      //K-
      h_K_minus_purity_pfRICH_MC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_K_minus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin));

      if( h_K_minus_purity_pfRICH_MC[Q2bin][y_bin]->Integral() > 0 )
      {
        TCanvas *K_minus_purity_pfRICH_can = new TCanvas(Form("K_minus_purity_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_minus_purity_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), 1200, 1000);

        K_minus_purity_pfRICH_can->cd();


        h_K_minus_purity_pfRICH_MC[Q2bin][y_bin]->Draw();

        double nK = h_K_minus_purity_pfRICH_MC[Q2bin][y_bin]->GetBinContent(2);
        double nBackground = h_K_minus_purity_pfRICH_MC[Q2bin][y_bin]->GetBinContent(1);

        double K_purity = ( nK/(nK+nBackground) )*100;

        TPaveText *Text_K_minus_purity = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
        Text_K_minus_purity->SetTextFont(42);
        Text_K_minus_purity->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        Text_K_minus_purity->AddText(Form("K^{-} pfRICH purity: %0.1f%%", K_purity));
        //Text_K_minus_purity->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_minus_purity->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_minus_purity->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_minus_purity->SetFillColorAlpha(0, 0.01);
        Text_K_minus_purity->Draw("same");

        K_minus_purity_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_purity_pfRICH_MC_Q2_%i_yBin_%i.png", e_energy, p_energy, Q2bin, y_bin));

      }


      //K-
      h_K_plus_purity_pfRICH_MC[Q2bin][y_bin] = (TH1F*)inFile->Get(Form("h_K_plus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin));

      if( h_K_plus_purity_pfRICH_MC[Q2bin][y_bin]->Integral() > 0 )
      {
        TCanvas *K_plus_purity_pfRICH_can = new TCanvas(Form("K_plus_purity_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), Form("K_plus_purity_pfRICH_can_Q2_%i_y_%i" , Q2bin, y_bin), 1200, 1000);

        K_plus_purity_pfRICH_can->cd();


        h_K_plus_purity_pfRICH_MC[Q2bin][y_bin]->Draw();

        double nK = h_K_plus_purity_pfRICH_MC[Q2bin][y_bin]->GetBinContent(2);
        double nBackground = h_K_plus_purity_pfRICH_MC[Q2bin][y_bin]->GetBinContent(1);

        double K_purity = ( nK/(nK+nBackground) )*100;

        TPaveText *Text_K_plus_purity = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
        Text_K_plus_purity->SetTextFont(42);
        Text_K_plus_purity->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        Text_K_plus_purity->AddText(Form("K^{+} pfRICH purity: %0.1f%%", K_purity));
        //Text_K_plus_purity->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
        Text_K_plus_purity->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
        Text_K_plus_purity->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        Text_K_plus_purity->SetFillColorAlpha(0, 0.01);
        Text_K_plus_purity->Draw("same");

        K_plus_purity_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_plus_purity_pfRICH_MC_Q2_%i_yBin_%i.png", e_energy, p_energy , Q2bin, y_bin));

      }
    }
  }


  for(unsigned int eta_bin = 0; eta_bin < nEtaBins+1; eta_bin++)
  {

    //K purity histograms

    //K-
    h_K_minus_purity_p_eta_pfRICH_MC[eta_bin] = (TH1F*)inFile->Get(Form("h_K_minus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin));

    if(h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->Integral() != 0)
    {
      TCanvas *K_minus_purity_pfRICH_MC_can = new TCanvas(Form("K_minus_purity_pfRICH_MC_can_eta_%i" , eta_bin), Form("K_minus_purity_pfRICH_MC_can_eta_%i", eta_bin), 1200, 1000);

      K_minus_purity_pfRICH_MC_can->cd();

      h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->Draw();

      double nK = h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->GetBinContent(2);
      double nBackground = h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->GetBinContent(1);

      double K_minus_purity = ( nK/(nK+nBackground) )*100;

      TPaveText *Text_K_minus_minus_purity_MC = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
      Text_K_minus_minus_purity_MC->SetTextFont(42);
      Text_K_minus_minus_purity_MC->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
      Text_K_minus_minus_purity_MC->AddText(Form("K^{-} pfRICH purity: %0.1f%%", K_minus_purity));
      //Text_K_minus_minus_purity_MC->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
      if(eta_bin < nEtaBins) Text_K_minus_minus_purity_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[eta_bin], eta_bins[eta_bin+1]));
      else Text_K_minus_minus_purity_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[0], eta_bins[nEtaBins]));
      Text_K_minus_minus_purity_MC->SetFillColorAlpha(0, 0.01);
      Text_K_minus_minus_purity_MC->Draw("same");

      K_minus_purity_pfRICH_MC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_purity_pfRICH_MC_eta_%i.png", e_energy, p_energy, eta_bin));

    }


    //K+
    h_K_plus_purity_p_eta_pfRICH_MC[eta_bin] = (TH1F*)inFile->Get(Form("h_K_plus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin));

    if(h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->Integral() != 0 )
    {
      TCanvas *K_plus_purity_pfRICH_MC_can = new TCanvas(Form("K_plus_purity_pfRICH_MC_can_eta_%i" , eta_bin), Form("K_plus_purity_pfRICH_MC_can_eta_%i" , eta_bin), 1200, 1000);

      K_plus_purity_pfRICH_MC_can->cd();

      h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->Draw();

      double nK = h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->GetBinContent(2);
      double nBackground = h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->GetBinContent(1);

      double K_plus_purity = ( nK/(nK+nBackground) )*100;

      TPaveText *Text_K_plus_minus_purity_MC = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
      Text_K_plus_minus_purity_MC->SetTextFont(42);
      Text_K_plus_minus_purity_MC->AddText(Form("RC ep DIS %ix%i GeV", e_energy, p_energy));
      Text_K_plus_minus_purity_MC->AddText(Form("K pfRICH purity: %0.1f%%", K_plus_purity));
      //Text_K_plus_minus_purity_MC->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
      if(eta_bin < nEtaBins) Text_K_plus_minus_purity_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[eta_bin], eta_bins[eta_bin+1]));
      else Text_K_plus_minus_purity_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[0], eta_bins[nEtaBins]));
      Text_K_plus_minus_purity_MC->SetFillColorAlpha(0, 0.01);
      Text_K_plus_minus_purity_MC->Draw("same");

      K_plus_purity_pfRICH_MC_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_plus_purity_pfRICH_MC_eta_%i.png", e_energy, p_energy, eta_bin));

    }



    TPaveText *Text_K_p_MC = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
    Text_K_p_MC->SetTextFont(42);
    Text_K_p_MC->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
    //Text_K_p_MC->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins, mom_bins[mom_bin+1]));
    if(eta_bin < nEtaBins) Text_K_p_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[eta_bin], eta_bins[eta_bin+1]));
    else Text_K_p_MC->AddText(Form("%0.1f < #eta < %0.1f GeV/c", eta_bins[0], eta_bins[nEtaBins]));
    Text_K_p_MC->SetFillColorAlpha(0, 0.01);

    //K p spectra

    //K-
    h_p_K_minus_MC[eta_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_MC_eta_%i" , eta_bin));
    h_p_K_minus_MC_pfRICH[eta_bin] = (TH1F*)inFile->Get(Form("h_p_K_minus_MC_pfRICH_eta_%i" , eta_bin));

    if( h_p_K_minus_MC[eta_bin]->Integral() > 0 )
    {
      TCanvas *K_minus_p_pfRICH_can = new TCanvas(Form("K_minus_p_pfRICH_can_eta_%i" , eta_bin), Form("K_minus_p_pfRICH_can_eta_%i" , eta_bin), 1200, 1000);
      K_minus_p_pfRICH_can->cd();

      gPad->SetLogy();


      h_p_K_minus_MC[eta_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
      h_p_K_minus_MC[eta_bin]->GetXaxis()->CenterTitle();
      h_p_K_minus_MC[eta_bin]->GetYaxis()->SetTitle("Counts");
      h_p_K_minus_MC[eta_bin]->GetYaxis()->CenterTitle();
      h_p_K_minus_MC[eta_bin]->SetMarkerStyle(28);
      h_p_K_minus_MC[eta_bin]->SetMarkerColor(4);
      h_p_K_minus_MC[eta_bin]->SetLineColor(4);
      h_p_K_minus_MC[eta_bin]->Draw("p e");


      h_p_K_minus_MC_pfRICH[eta_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
      h_p_K_minus_MC_pfRICH[eta_bin]->GetXaxis()->CenterTitle();
      h_p_K_minus_MC_pfRICH[eta_bin]->GetYaxis()->SetTitle("Counts");
      h_p_K_minus_MC_pfRICH[eta_bin]->GetYaxis()->CenterTitle();
      h_p_K_minus_MC_pfRICH[eta_bin]->SetMarkerStyle(34);
      h_p_K_minus_MC_pfRICH[eta_bin]->SetMarkerColor(8);
      h_p_K_minus_MC_pfRICH[eta_bin]->SetLineColor(8);
      h_p_K_minus_MC_pfRICH[eta_bin]->Draw("p e same");

      Text_K_p_MC->Draw("same");

      K_minus_p_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_minus_p_eta_%i.png" , e_energy, p_energy, eta_bin));

    }

    //K+
    h_p_K_plus_MC[eta_bin] = (TH1F*)inFile->Get(Form("h_p_K_plus_MC_eta_%i" , eta_bin));
    h_p_K_plus_MC_pfRICH[eta_bin] = (TH1F*)inFile->Get(Form("h_p_K_plus_MC_pfRICH_eta_%i" , eta_bin));

    if( h_p_K_plus_MC[eta_bin]->Integral() > 0 )
    {

      TCanvas *K_plus_p_pfRICH_can = new TCanvas(Form("K_plus_p_pfRICH_can_eta_%i" , eta_bin), Form("K_plus_p_pfRICH_can_eta_%i" , eta_bin), 1200, 1000);
      K_plus_p_pfRICH_can->cd();

      gPad->SetLogy();

      h_p_K_plus_MC[eta_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
      h_p_K_plus_MC[eta_bin]->GetXaxis()->CenterTitle();
      h_p_K_plus_MC[eta_bin]->GetYaxis()->SetTitle("Counts");
      h_p_K_plus_MC[eta_bin]->GetYaxis()->CenterTitle();
      h_p_K_plus_MC[eta_bin]->SetMarkerStyle(28);
      h_p_K_plus_MC[eta_bin]->SetMarkerColor(4);
      h_p_K_plus_MC[eta_bin]->SetLineColor(4);
      h_p_K_plus_MC[eta_bin]->Draw("p e");


      h_p_K_plus_MC_pfRICH[eta_bin]->GetXaxis()->SetTitle("p (GeV/#it{c})");
      h_p_K_plus_MC_pfRICH[eta_bin]->GetXaxis()->CenterTitle();
      h_p_K_plus_MC_pfRICH[eta_bin]->GetYaxis()->SetTitle("Counts");
      h_p_K_plus_MC_pfRICH[eta_bin]->GetYaxis()->CenterTitle();
      h_p_K_plus_MC_pfRICH[eta_bin]->SetMarkerStyle(34);
      h_p_K_plus_MC_pfRICH[eta_bin]->SetMarkerColor(8);
      h_p_K_plus_MC_pfRICH[eta_bin]->SetLineColor(8);
      h_p_K_plus_MC_pfRICH[eta_bin]->Draw("p e same");

      Text_K_p_MC->Draw("same");

      K_plus_p_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_plus_p_eta_%i.png", e_energy, p_energy, eta_bin));

    }


  }




  //fitPiTOF->Write();
  //fitKTOF->Write();

  inFile->Close();

  return;

}
