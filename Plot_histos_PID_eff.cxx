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

void Plot_histos_PID_eff(int e_energy = 18, int p_energy = 275)
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
  gStyle->SetLineScalePS(1);//for pdf lines

  //____________________________________________________


  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 11;
  float const mom_bins[nMomBins+1] = { 0,0.5,1,1.5,2,3,4,5,6,7,10, 18 };

  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -4, -3, -2, -1, 0 };

  //load all files
  //TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output_QA.root", e_energy, p_energy, e_energy, p_energy), "READ");
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output_QA_e_pi_smearing.root", e_energy, p_energy, e_energy, p_energy), "READ");
  //TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i_p_scat_bins-output.root", e_energy, p_energy, e_energy, p_energy), "READ");


  TFile *outFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/output/%ix%i/DIS_%ix%i-PID_table_QA.root", e_energy, p_energy, e_energy, p_energy), "recreate");



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

  //e/pi PID efficiency histograms
  TH1D *h_e_pfRICH_PID_eff_MC[3];
  TH1D *h_e_pi_pfRICH_PID_eff_MC[3];

  TH1D *h_e_pfRICH_PID_eff_RC[3];
  TH1D *h_e_pi_pfRICH_PID_eff_RC[3];

  TH1D *h_e_pfRICH_PID_eff_MC_RC[3];
  TH1D *h_e_pi_pfRICH_PID_eff_MC_RC[3];


  //e/pi PID TOF efficiency histograms
  TH1D *h_e_pfRICH_TOF_PID_eff_MC[3];
  TH1D *h_e_pi_pfRICH_TOF_PID_eff_MC[3];

  //TH1D *h_e_pfRICH_TOF_PID_eff_RC[3];
  //TH1D *h_e_pi_pfRICH_TOF_PID_eff_RC[3];

  //TH1D *h_e_pfRICH_TOF_PID_eff_MC_RC[3];
  //TH1D *h_e_pi_pfRICH_TOF_PID_eff_MC_RC[3];


  //baseline spectra for PID efficiency calculation
  //for QA
  TCanvas *pfRICH_ref_spectra = new TCanvas("pfRICH_ref_spectra", "pfRICH_ref_spectra", 3600, 3000);
  pfRICH_ref_spectra->Divide(3,3);

  TPaveText *pfRICH_ref_spectra_pi_text = new TPaveText(0.55, 0.75, 0.75, 0.8, "NDC");
  pfRICH_ref_spectra_pi_text->AddText("#pi");
  pfRICH_ref_spectra_pi_text->SetTextFont(42);
  pfRICH_ref_spectra_pi_text->SetFillColorAlpha(0, 0.01);

  TPaveText *pfRICH_ref_spectra_K_text = new TPaveText(0.55, 0.75, 0.75, 0.8, "NDC");
  pfRICH_ref_spectra_K_text->AddText("K");
  pfRICH_ref_spectra_K_text->SetTextFont(42);
  pfRICH_ref_spectra_K_text->SetFillColorAlpha(0, 0.01);

  TPaveText *pfRICH_ref_spectra_p_text = new TPaveText(0.55, 0.75, 0.75, 0.8, "NDC");
  pfRICH_ref_spectra_p_text->AddText("p");
  pfRICH_ref_spectra_p_text->SetTextFont(42);
  pfRICH_ref_spectra_p_text->SetFillColorAlpha(0, 0.01);

  TPaveText *pfRICH_PID_eff_pi_text[3];
  TPaveText *pfRICH_PID_eff_K_text[3];
  TPaveText *pfRICH_PID_eff_p_text[3];


  TPaveText *pfRICH_e_pi_PID_eff_e_text[2];
  TPaveText *pfRICH_e_pi_PID_eff_pi_text[2];


  for(unsigned int PID_bin = 0; PID_bin < 4; PID_bin++)
  {
    h_pi_pfRICH_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_pi_pfRICH_PID_eff_RC_%i", PID_bin)); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_K_pfRICH_PID_eff_RC_%i", PID_bin));
    h_p_pfRICH_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_p_pfRICH_PID_eff_RC_%i", PID_bin));

    h_pi_pfRICH_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_pi_pfRICH_PID_eff_MC_%i", PID_bin)); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_K_pfRICH_PID_eff_MC_%i", PID_bin));
    h_p_pfRICH_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_p_pfRICH_PID_eff_MC_%i", PID_bin));

    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin] = (TH1D*)inFile->Get(Form("h_pi_pfRICH_PID_eff_MC_RC_match_%i", PID_bin)); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin] = (TH1D*)inFile->Get(Form("h_K_pfRICH_PID_eff_MC_RC_match_%i", PID_bin));
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin] = (TH1D*)inFile->Get(Form("h_p_pfRICH_PID_eff_MC_RC_match_%i", PID_bin));

    if(PID_bin > 0)
    {
      //PID_bin go from 1 to 3, because PID_bin 0 is baseline without any PID
      //need to shift index for TPaveText by -1
      if(PID_bin == 2) pfRICH_PID_eff_pi_text[PID_bin-1] = new TPaveText(0.25, 0.45, 0.54, 0.6, "NDC");
      else pfRICH_PID_eff_pi_text[PID_bin-1] = new TPaveText(0.55, 0.15, 0.84, 0.3, "NDC");
      pfRICH_PID_eff_pi_text[PID_bin-1]->SetTextFont(42);
      pfRICH_PID_eff_pi_text[PID_bin-1]->SetFillColorAlpha(0, 0.01);

      if(PID_bin == 2) pfRICH_PID_eff_K_text[PID_bin-1] = new TPaveText(0.25, 0.45, 0.54, 0.6, "NDC");
      else pfRICH_PID_eff_K_text[PID_bin-1] = new TPaveText(0.55, 0.15, 0.84, 0.3, "NDC");
      pfRICH_PID_eff_K_text[PID_bin-1]->SetTextFont(42);
      pfRICH_PID_eff_K_text[PID_bin-1]->SetFillColorAlpha(0, 0.01);

      pfRICH_PID_eff_p_text[PID_bin-1] = new TPaveText(0.55, 0.45, 0.84, 0.6, "NDC");
      pfRICH_PID_eff_p_text[PID_bin-1]->SetTextFont(42);
      pfRICH_PID_eff_p_text[PID_bin-1]->SetFillColorAlpha(0, 0.01);

      if(PID_bin == 1)
      {
        pfRICH_PID_eff_pi_text[PID_bin-1]->AddText("MC #pi identified as #pi");
        pfRICH_PID_eff_K_text[PID_bin-1]->AddText("MC K identified as K");
        pfRICH_PID_eff_p_text[PID_bin-1]->AddText("MC p identified as p");
      }

      if(PID_bin == 2)
      {
        pfRICH_PID_eff_pi_text[PID_bin-1]->AddText("MC #pi mis-identified as K");
        pfRICH_PID_eff_K_text[PID_bin-1]->AddText("MC K mis-identified as #pi");
        pfRICH_PID_eff_p_text[PID_bin-1]->AddText("MC p mis-identified as #pi");

      }

      if(PID_bin == 3)
      {
        pfRICH_PID_eff_pi_text[PID_bin-1]->AddText("MC #pi mis-identified as p");
        pfRICH_PID_eff_K_text[PID_bin-1]->AddText("MC K mis-identified as p");
        pfRICH_PID_eff_p_text[PID_bin-1]->AddText("MC p mis-identified as K");
      }

    }



    //e/pi PID table
    if(PID_bin < 3)
    {
      //RICH PID
      h_e_pfRICH_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_PID_eff_MC_%i", PID_bin));
      h_e_pi_pfRICH_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_PID_eff_MC_%i", PID_bin));

      h_e_pfRICH_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_PID_eff_RC_%i", PID_bin));
      h_e_pi_pfRICH_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_PID_eff_RC_%i", PID_bin));

      h_e_pfRICH_PID_eff_MC_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_PID_eff_MC_RC_%i", PID_bin));
      h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_PID_eff_MC_RC_%i", PID_bin));

      //TOF PID
      h_e_pfRICH_TOF_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_TOF_PID_eff_MC_%i", PID_bin));
      h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_TOF_PID_eff_MC_%i", PID_bin));

      //h_e_pfRICH_TOF_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_TOF_PID_eff_RC_%i", PID_bin));
      //h_e_pi_pfRICH_TOF_PID_eff_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_TOF_PID_eff_RC_%i", PID_bin));

      //h_e_pfRICH_TOF_PID_eff_MC_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin));
      //h_e_pi_pfRICH_TOF_PID_eff_MC_RC[PID_bin] = (TH1D*)inFile->Get(Form("h_e_pi_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin));

      if(PID_bin > 0)
      {
        pfRICH_e_pi_PID_eff_e_text[PID_bin-1] = new TPaveText(0.55, 0.15, 0.84, 0.2, "NDC");
        pfRICH_e_pi_PID_eff_e_text[PID_bin-1]->SetTextFont(42);
        pfRICH_e_pi_PID_eff_e_text[PID_bin-1]->SetFillColorAlpha(0, 0.01);

        pfRICH_e_pi_PID_eff_pi_text[PID_bin-1] = new TPaveText(0.55, 0.15, 0.84, 0.2, "NDC");
        pfRICH_e_pi_PID_eff_pi_text[PID_bin-1]->SetTextFont(42);
        pfRICH_e_pi_PID_eff_pi_text[PID_bin-1]->SetFillColorAlpha(0, 0.01);

        if(PID_bin == 1)
        {
          pfRICH_e_pi_PID_eff_e_text[PID_bin-1]->AddText("e #rightarrow e");
          pfRICH_e_pi_PID_eff_pi_text[PID_bin-1]->AddText("pi #rightarrow pi");
        }

        if(PID_bin == 2)
        {
          pfRICH_e_pi_PID_eff_e_text[PID_bin-1]->AddText("e #rightarrow pi");
          pfRICH_e_pi_PID_eff_pi_text[PID_bin-1]->AddText("pi #rightarrow e");
        }
      }

    }

  }


  //_____________________________________________________________________________________________


  //pi/K/p PID efficiencies
  TCanvas *pfRICH_PID_eff_MC = new TCanvas("pfRICH_PID_eff_MC", "pfRICH_PID_eff_MC", 3600, 3000);
  pfRICH_PID_eff_MC->Divide(3,3);

  TCanvas *pfRICH_PID_eff_RC = new TCanvas("pfRICH_PID_eff_RC", "pfRICH_PID_eff_RC", 3600, 3000);
  pfRICH_PID_eff_RC->Divide(3,3);

  TCanvas *pfRICH_PID_eff_MC_RC_match = new TCanvas("pfRICH_PID_eff_MC_RC_match", "pfRICH_PID_eff_MC_RC_match", 3600, 3000);
  pfRICH_PID_eff_MC_RC_match->Divide(3,3);






  //bin 0 is the baseline without PID
  for(unsigned int PID_bin = 1; PID_bin < 4; PID_bin++)
  {



    pfRICH_PID_eff_MC->cd(PID_bin);

    if(PID_bin == 1) h_pi_pfRICH_PID_eff_MC[0]->Sumw2();
    h_pi_pfRICH_PID_eff_MC[PID_bin]->Sumw2();
    h_pi_pfRICH_PID_eff_MC[PID_bin]->Divide(h_pi_pfRICH_PID_eff_MC[PID_bin], h_pi_pfRICH_PID_eff_MC[0], 1, 1, "b");
    //h_pi_pfRICH_PID_eff_MC[PID_bin]->Divide(h_pi_pfRICH_PID_eff_MC[0]);
    h_pi_pfRICH_PID_eff_MC[PID_bin]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_pi_pfRICH_PID_eff_MC[PID_bin]->GetXaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_MC[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_pi_pfRICH_PID_eff_MC[PID_bin]->GetYaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_MC[PID_bin]->SetMarkerStyle(20);
    h_pi_pfRICH_PID_eff_MC[PID_bin]->SetLineColor(1);
    h_pi_pfRICH_PID_eff_MC[PID_bin]->Draw("p e");

    pfRICH_PID_eff_pi_text[PID_bin-1]->Draw("same");



    pfRICH_PID_eff_MC->cd(PID_bin+3);

    int PID_bin_K;

    if( PID_bin == 1) PID_bin_K = 2;
    if( PID_bin == 2) PID_bin_K = 1;
    if( PID_bin == 3) PID_bin_K = 3;

    if(PID_bin == 1) h_K_pfRICH_PID_eff_MC[0]->Sumw2();
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->Sumw2();
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->Divide(h_K_pfRICH_PID_eff_MC[PID_bin_K], h_K_pfRICH_PID_eff_MC[0], 1, 1, "b");
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->GetXaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->GetYaxis()->SetTitle("Efficiency");
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->GetYaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->SetMarkerStyle(20);
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->SetLineColor(1);
    h_K_pfRICH_PID_eff_MC[PID_bin_K]->Draw("p e");

    pfRICH_PID_eff_K_text[PID_bin_K-1]->Draw("same");


    pfRICH_PID_eff_MC->cd(PID_bin+6);


    if(PID_bin == 1) h_p_pfRICH_PID_eff_MC[0]->Sumw2();



    int PID_bin_p;

    if( PID_bin == 1) PID_bin_p = 2;
    if( PID_bin == 2) PID_bin_p = 3;
    if( PID_bin == 3) PID_bin_p = 1;

    h_p_pfRICH_PID_eff_MC[PID_bin_p]->Sumw2();
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->Divide(h_p_pfRICH_PID_eff_MC[PID_bin_p], h_p_pfRICH_PID_eff_MC[0], 1, 1, "b");
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->GetXaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->GetYaxis()->SetTitle("Efficiency");
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->GetYaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->SetMarkerStyle(20);
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->SetLineColor(1);
    h_p_pfRICH_PID_eff_MC[PID_bin_p]->Draw("p e");

    pfRICH_PID_eff_p_text[PID_bin_p-1]->Draw("same");


    //__________________________________________________________________________________

    pfRICH_PID_eff_RC->cd(PID_bin);



    if(PID_bin == 1) h_pi_pfRICH_PID_eff_RC[0]->Sumw2();
    h_pi_pfRICH_PID_eff_RC[PID_bin]->Sumw2();
    h_pi_pfRICH_PID_eff_RC[PID_bin]->Divide(h_pi_pfRICH_PID_eff_RC[PID_bin], h_pi_pfRICH_PID_eff_RC[0], 1, 1, "b");
    h_pi_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
    h_pi_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_pi_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_RC[PID_bin]->Draw("p e");

    pfRICH_PID_eff_pi_text[PID_bin-1]->Draw("same");


    pfRICH_PID_eff_RC->cd(PID_bin+3);

    if(PID_bin == 1) h_K_pfRICH_PID_eff_RC[0]->Sumw2();
    h_K_pfRICH_PID_eff_RC[PID_bin]->Sumw2();
    h_K_pfRICH_PID_eff_RC[PID_bin]->Divide(h_K_pfRICH_PID_eff_RC[PID_bin], h_K_pfRICH_PID_eff_RC[0], 1, 1, "b");
    h_K_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
    h_K_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_K_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_RC[PID_bin]->Draw("p e");

    pfRICH_PID_eff_K_text[PID_bin-1]->Draw("same");


    pfRICH_PID_eff_RC->cd(PID_bin+6);

    if(PID_bin == 1) h_p_pfRICH_PID_eff_RC[0]->Sumw2();
    h_p_pfRICH_PID_eff_RC[PID_bin]->Sumw2();
    h_p_pfRICH_PID_eff_RC[PID_bin]->Divide(h_p_pfRICH_PID_eff_RC[PID_bin], h_p_pfRICH_PID_eff_RC[0], 1, 1, "b");
    h_p_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
    h_p_pfRICH_PID_eff_RC[PID_bin]->GetXaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_p_pfRICH_PID_eff_RC[PID_bin]->GetYaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_RC[PID_bin]->Draw("p e");

    pfRICH_PID_eff_p_text[PID_bin-1]->Draw("same");

    //_____________________________________________________________________________________________________________

    pfRICH_PID_eff_MC_RC_match->cd(PID_bin);

    if(PID_bin == 1) h_pi_pfRICH_PID_eff_MC_RC_match[0]->Sumw2();
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->Sumw2();
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->Divide(h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin], h_pi_pfRICH_PID_eff_MC_RC_match[0], 1, 1, "b");
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->SetTitle("p_{MC_RC} (GeV/#it{c})");
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->CenterTitle();
    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin]->Draw("p e");

    pfRICH_PID_eff_pi_text[PID_bin-1]->Draw("same");



    pfRICH_PID_eff_MC_RC_match->cd(PID_bin+3);

    if(PID_bin == 1) h_K_pfRICH_PID_eff_MC_RC_match[0]->Sumw2();
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->Sumw2();
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->Divide(h_K_pfRICH_PID_eff_MC_RC_match[PID_bin], h_K_pfRICH_PID_eff_MC_RC_match[0], 1, 1, "b");
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->SetTitle("p_{MC_RC} (GeV/#it{c})");
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->CenterTitle();
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin]->Draw("p e");

    pfRICH_PID_eff_K_text[PID_bin-1]->Draw("same");


    pfRICH_PID_eff_MC_RC_match->cd(PID_bin+6);

    if(PID_bin == 1) h_p_pfRICH_PID_eff_MC_RC_match[0]->Sumw2();
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->Sumw2();
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->Divide(h_p_pfRICH_PID_eff_MC_RC_match[PID_bin], h_p_pfRICH_PID_eff_MC_RC_match[0], 1, 1, "b");
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->SetTitle("p_{MC_RC} (GeV/#it{c})");
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetXaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->SetTitle("Efficiency");
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->GetYaxis()->CenterTitle();
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin]->Draw("p e");

    pfRICH_PID_eff_p_text[PID_bin-1]->Draw("same");


    //__________________________________________________________________________________

    //plot baseline spectra for PID efficiency calculation


    pfRICH_ref_spectra->cd(PID_bin);

    gPad->SetLogy();

    if(PID_bin == 1)
    {
      h_pi_pfRICH_PID_eff_MC[0]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
      h_pi_pfRICH_PID_eff_MC[0]->GetXaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_MC[0]->GetYaxis()->SetTitle("Counts");
      h_pi_pfRICH_PID_eff_MC[0]->GetYaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_MC[0]->Draw();

      pfRICH_ref_spectra_pi_text->Draw("same");

    }


    if(PID_bin == 2)
    {
      h_pi_pfRICH_PID_eff_RC[0]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
      h_pi_pfRICH_PID_eff_RC[0]->GetXaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_RC[0]->GetYaxis()->SetTitle("Counts");
      h_pi_pfRICH_PID_eff_RC[0]->GetYaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_RC[0]->Draw();

      pfRICH_ref_spectra_pi_text->Draw("same");

    }

    if(PID_bin == 3)
    {
      h_pi_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->SetTitle("p_{MC-RC} (GeV/#it{c})");
      h_pi_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->SetTitle("Counts");
      h_pi_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->CenterTitle();
      h_pi_pfRICH_PID_eff_MC_RC_match[0]->Draw();

      pfRICH_ref_spectra_pi_text->Draw("same");
    }


    pfRICH_ref_spectra->cd(PID_bin+3);

    gPad->SetLogy();

    if(PID_bin == 1)
    {
      h_K_pfRICH_PID_eff_MC[0]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
      h_K_pfRICH_PID_eff_MC[0]->GetXaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_MC[0]->GetYaxis()->SetTitle("Counts");
      h_K_pfRICH_PID_eff_MC[0]->GetYaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_MC[0]->Draw();

      pfRICH_ref_spectra_K_text->Draw("same");

    }


    if(PID_bin == 2)
    {
      h_K_pfRICH_PID_eff_RC[0]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
      h_K_pfRICH_PID_eff_RC[0]->GetXaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_RC[0]->GetYaxis()->SetTitle("Counts");
      h_K_pfRICH_PID_eff_RC[0]->GetYaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_RC[0]->Draw();

      pfRICH_ref_spectra_K_text->Draw("same");

    }

    if(PID_bin == 3)
    {
      h_K_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->SetTitle("p_{MC-RC} (GeV/#it{c})");
      h_K_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->SetTitle("Counts");
      h_K_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->CenterTitle();
      h_K_pfRICH_PID_eff_MC_RC_match[0]->Draw();

      pfRICH_ref_spectra_K_text->Draw("same");
    }


    pfRICH_ref_spectra->cd(PID_bin+6);

    gPad->SetLogy();

    if(PID_bin == 1)
    {
      h_p_pfRICH_PID_eff_MC[0]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
      h_p_pfRICH_PID_eff_MC[0]->GetXaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_MC[0]->GetYaxis()->SetTitle("Counts");
      h_p_pfRICH_PID_eff_MC[0]->GetYaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_MC[0]->Draw();

      pfRICH_ref_spectra_p_text->Draw("same");

    }


    if(PID_bin == 2)
    {
      h_p_pfRICH_PID_eff_RC[0]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
      h_p_pfRICH_PID_eff_RC[0]->GetXaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_RC[0]->GetYaxis()->SetTitle("Counts");
      h_p_pfRICH_PID_eff_RC[0]->GetYaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_RC[0]->Draw();

      pfRICH_ref_spectra_p_text->Draw("same");

    }

    if(PID_bin == 3)
    {
      h_p_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->SetTitle("p_{MC-RC} (GeV/#it{c})");
      h_p_pfRICH_PID_eff_MC_RC_match[0]->GetXaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->SetTitle("Counts");
      h_p_pfRICH_PID_eff_MC_RC_match[0]->GetYaxis()->CenterTitle();
      h_p_pfRICH_PID_eff_MC_RC_match[0]->Draw();

      pfRICH_ref_spectra_p_text->Draw("same");
    }

  }

  pfRICH_ref_spectra->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_ref_spectra.png", e_energy, p_energy));

  pfRICH_PID_eff_MC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_PID_eff_MC.png", e_energy, p_energy));
  pfRICH_PID_eff_MC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_PID_eff_MC.pdf", e_energy, p_energy));
  pfRICH_PID_eff_RC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_PID_eff_RC.png", e_energy, p_energy));
  pfRICH_PID_eff_MC_RC_match->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_PID_eff_MC_RC_match.png", e_energy, p_energy));


  //plot 2D MC vs. RC p QA histograms
  TCanvas *pi_pfRICH_MC_p_RC_p_can = new TCanvas("pi_pfRICH_MC_p_RC_p_can", "pi_pfRICH_MC_p_RC_p_can", 1200, 1000);
  pi_pfRICH_MC_p_RC_p_can->cd();

  gPad->SetLogz();

  h_pi_pfRICH_MC_p_RC_p->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
  h_pi_pfRICH_MC_p_RC_p->GetXaxis()->CenterTitle();
  h_pi_pfRICH_MC_p_RC_p->GetYaxis()->SetTitle("p_{RC} (GeV/#it{c})");
  h_pi_pfRICH_MC_p_RC_p->GetYaxis()->CenterTitle();
  h_pi_pfRICH_MC_p_RC_p->Draw("colz");


  pfRICH_ref_spectra_pi_text->Draw("same");

  pi_pfRICH_MC_p_RC_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_pi_MC_p_vs_RC_p.png", e_energy, p_energy));


  TCanvas *K_pfRICH_MC_p_RC_p_can = new TCanvas("K_pfRICH_MC_p_RC_p_can", "K_pfRICH_MC_p_RC_p_can", 1200, 1000);
  K_pfRICH_MC_p_RC_p_can->cd();

  gPad->SetLogz();

  h_K_pfRICH_MC_p_RC_p->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
  h_K_pfRICH_MC_p_RC_p->GetXaxis()->CenterTitle();
  h_K_pfRICH_MC_p_RC_p->GetYaxis()->SetTitle("p_{RC} (GeV/#it{c})");
  h_K_pfRICH_MC_p_RC_p->GetYaxis()->CenterTitle();
  h_K_pfRICH_MC_p_RC_p->Draw("colz");

  pfRICH_ref_spectra_K_text->Draw("same");

  K_pfRICH_MC_p_RC_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_K_MC_p_vs_RC_p.png", e_energy, p_energy));


  TCanvas *p_pfRICH_MC_p_RC_p_can = new TCanvas("p_pfRICH_MC_p_RC_p_can", "p_pfRICH_MC_p_RC_p_can", 1200, 1000);
  p_pfRICH_MC_p_RC_p_can->cd();

  gPad->SetLogz();

  h_p_pfRICH_MC_p_RC_p->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
  h_p_pfRICH_MC_p_RC_p->GetXaxis()->CenterTitle();
  h_p_pfRICH_MC_p_RC_p->GetYaxis()->SetTitle("p_{RC} (GeV/#it{c})");
  h_p_pfRICH_MC_p_RC_p->GetYaxis()->CenterTitle();
  h_p_pfRICH_MC_p_RC_p->Draw("colz");

  pfRICH_ref_spectra_p_text->Draw("same");

  p_pfRICH_MC_p_RC_p_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_p_MC_p_vs_RC_p.png", e_energy, p_energy));


  //____________________________________________________________________________________________________________________________________________________________



  //e/pi RICH PID efficiencies
  TCanvas *pfRICH_e_pi_PID_eff_MC = new TCanvas("pfRICH_e_pi_PID_eff_MC", "pfRICH_e_pi_PID_eff_MC", 3600, 3000);
  pfRICH_e_pi_PID_eff_MC->Divide(2,2);

  TCanvas *pfRICH_e_pi_PID_eff_RC = new TCanvas("pfRICH_e_pi_PID_eff_RC", "pfRICH_e_pi_PID_eff_RC", 3600, 3000);
  pfRICH_e_pi_PID_eff_RC->Divide(2,2);

  TCanvas *pfRICH_e_pi_PID_eff_MC_RC_match = new TCanvas("pfRICH_e_pi_PID_eff_MC_RC_match", "pfRICH_e_pi_PID_eff_MC_RC_match", 3600, 3000);
  pfRICH_e_pi_PID_eff_MC_RC_match->Divide(2,2);


  //e/pi TOF PID efficiencies
  TCanvas *pfRICH_TOF_e_pi_PID_eff_MC = new TCanvas("pfRICH_TOF_e_pi_PID_eff_MC", "pfRICH_TOF_e_pi_PID_eff_MC", 3600, 3000);
  pfRICH_TOF_e_pi_PID_eff_MC->Divide(2,2);

  //TCanvas *pfRICH_TOF_e_pi_PID_eff_RC = new TCanvas("pfRICH_TOF_e_pi_PID_eff_RC", "pfRICH_TOF_e_pi_PID_eff_RC", 3600, 3000);
  //pfRICH_TOF_e_pi_PID_eff_RC->Divide(2,2);

  //TCanvas *pfRICH_TOF_e_pi_PID_eff_MC_RC_match = new TCanvas("pfRICH_TOF_e_pi_PID_eff_MC_RC_match", "pfRICH_TOF_e_pi_PID_eff_MC_RC_match", 3600, 3000);
  //pfRICH_TOF_e_pi_PID_eff_MC_RC_match->Divide(2,2);



  int PID_bin_e_pi = -1;

  for(unsigned int PID_bin_e = 1; PID_bin_e < 3; PID_bin_e++)
  {
    //e/pi RICH PID efficiencies



    pfRICH_e_pi_PID_eff_MC->cd(PID_bin_e);

    if(PID_bin_e == 1) h_e_pfRICH_PID_eff_MC[0]->Sumw2();
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->Sumw2();
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->Divide(h_e_pfRICH_PID_eff_MC[PID_bin_e], h_e_pfRICH_PID_eff_MC[0], 1, 1, "b");
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->GetXaxis()->SetRangeUser(0, 8);
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->SetMarkerStyle(20);
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->SetMarkerSize(3);
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->SetMarkerColor(1);
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->SetLineColor(1);
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->Draw("p e");
    h_e_pfRICH_PID_eff_MC[PID_bin_e]->Write();

    pfRICH_e_pi_PID_eff_e_text[PID_bin_e-1]->Draw("same");



    pfRICH_e_pi_PID_eff_MC->cd(PID_bin_e+2);

    if( PID_bin_e == 1) PID_bin_e_pi = 2;
    if( PID_bin_e == 2) PID_bin_e_pi = 1;

    if(PID_bin_e == 1) h_e_pi_pfRICH_PID_eff_MC[0]->Sumw2();
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->Sumw2();
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->Divide(h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi], h_e_pi_pfRICH_PID_eff_MC[0], 1, 1, "b");
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->GetXaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->GetXaxis()->SetRangeUser(0, 8);
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->GetYaxis()->SetTitle("Efficiency");
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->GetYaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->SetMarkerStyle(20);
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->SetMarkerSize(3);
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->SetMarkerColor(1);
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->SetLineColor(1);
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->Draw("p e");
    h_e_pi_pfRICH_PID_eff_MC[PID_bin_e_pi]->Write();

    pfRICH_e_pi_PID_eff_pi_text[PID_bin_e_pi-1]->Draw("same");


    //__________________________________________________________________________________


    pfRICH_e_pi_PID_eff_RC->cd(PID_bin_e);

    if(PID_bin_e == 1) h_e_pfRICH_PID_eff_RC[0]->Sumw2();
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->Sumw2();
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->Divide(h_e_pfRICH_PID_eff_RC[PID_bin_e], h_e_pfRICH_PID_eff_RC[0], 1, 1, "b");
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_RC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_e_text[PID_bin_e-1]->Draw("same");


    pfRICH_e_pi_PID_eff_RC->cd(PID_bin_e+2);

    if(PID_bin_e == 1) h_e_pi_pfRICH_PID_eff_RC[0]->Sumw2();
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->Sumw2();
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->Divide(h_e_pi_pfRICH_PID_eff_RC[PID_bin_e], h_e_pi_pfRICH_PID_eff_RC[0], 1, 1, "b");
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->GetXaxis()->SetTitle("p_{RC} (GeV/#it{c})");
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_RC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_pi_text[PID_bin_e-1]->Draw("same");

    //__________________________________________________________________________________

    pfRICH_e_pi_PID_eff_MC_RC_match->cd(PID_bin_e);

    if(PID_bin_e == 1) h_e_pfRICH_PID_eff_MC_RC[0]->Sumw2();
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->Sumw2();
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->Divide(h_e_pfRICH_PID_eff_MC_RC[PID_bin_e], h_e_pfRICH_PID_eff_MC_RC[0], 1, 1, "b");
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetXaxis()->SetTitle("p_{MC_RC} (GeV/#it{c})");
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pfRICH_PID_eff_MC_RC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_e_text[PID_bin_e-1]->Draw("same");


    pfRICH_e_pi_PID_eff_MC_RC_match->cd(PID_bin_e+2);

    if(PID_bin_e == 1) h_e_pi_pfRICH_PID_eff_MC_RC[0]->Sumw2();
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->Sumw2();
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->Divide(h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e], h_e_pi_pfRICH_PID_eff_MC_RC[0], 1, 1, "b");
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetXaxis()->SetTitle("p_{MC_RC} (GeV/#it{c})");
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_pi_text[PID_bin_e-1]->Draw("same");

    //_________________________________________________________________________________________________________________________________________________

    //e/pi TOF PID efficiencies




    pfRICH_TOF_e_pi_PID_eff_MC->cd(PID_bin_e);

    if(PID_bin_e == 1) h_e_pfRICH_TOF_PID_eff_MC[0]->Sumw2();
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Sumw2();
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Divide(h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e], h_e_pfRICH_TOF_PID_eff_MC[0], 1, 1, "b");
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_e_text[PID_bin_e-1]->Draw("same");


    pfRICH_TOF_e_pi_PID_eff_MC->cd(PID_bin_e+2);

    if(PID_bin_e == 1) h_e_pi_pfRICH_TOF_PID_eff_MC[0]->Sumw2();
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Sumw2();
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Divide(h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e], h_e_pi_pfRICH_TOF_PID_eff_MC[0], 1, 1, "b");
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetXaxis()->SetTitle("p_{MC} (GeV/#it{c})");
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetXaxis()->CenterTitle();
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetYaxis()->SetTitle("Efficiency");
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->GetYaxis()->CenterTitle();
    h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin_e]->Draw("p e");

    pfRICH_e_pi_PID_eff_pi_text[PID_bin_e-1]->Draw("same");

    //__________________________________________________________________________________

  }

  pfRICH_e_pi_PID_eff_MC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_e_pi_PID_eff_MC.png", e_energy, p_energy));
  pfRICH_e_pi_PID_eff_MC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_e_pi_PID_eff_MC.pdf", e_energy, p_energy));

  pfRICH_e_pi_PID_eff_RC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_e_pi_PID_eff_RC.png", e_energy, p_energy));
  pfRICH_e_pi_PID_eff_MC_RC_match->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_e_pi_PID_eff_MC_RC.png", e_energy, p_energy));


  pfRICH_TOF_e_pi_PID_eff_MC->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/pi_reject_reco_pfRICH/pfRICH_TOF_e_pi_PID_eff_MC.png", e_energy, p_energy));



  outFile->Close();
  inFile->Close();

  return;

}
