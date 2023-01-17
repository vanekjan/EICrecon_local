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

  //load all files
  TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i-output.root", e_energy, p_energy, e_energy, p_energy), "READ");
  //TFile *inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/input/%ix%i/DIS_%ix%i_p_scat_bins-output.root", e_energy, p_energy, e_energy, p_energy), "READ");

  TH1D *h_K_MC_RC_match_pfRICH[nMomBins][nQ2bins][nyInelParBins];

  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {
        //K purity histograms
        TCanvas *K_MC_RC_match_pfRICH_can = new TCanvas(Form("K_MC_RC_match_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("K_MC_RC_match_pfRICH_can_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 1200, 1000);

        K_MC_RC_match_pfRICH_can->cd();

        h_K_MC_RC_match_pfRICH[mom_bin][Q2bin][y_bin] = (TH1D*)inFile->Get(Form("h_K_purity_pfRICH_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin));
        if( h_K_MC_RC_match_pfRICH[mom_bin][Q2bin][y_bin]->Integral() == 0 ) continue;
        h_K_MC_RC_match_pfRICH[mom_bin][Q2bin][y_bin]->Draw();

        double nK = h_K_MC_RC_match_pfRICH[mom_bin][Q2bin][y_bin]->GetBinContent(2);
        double nBackground = h_K_MC_RC_match_pfRICH[mom_bin][Q2bin][y_bin]->GetBinContent(1);

        double K_purity = ( nK/(nK+nBackground) )*100;

        TPaveText *Text_mathc_pfRICH = new TPaveText(0.55, 0.65, 0.84, 0.85, "NDC");
        Text_mathc_pfRICH->SetTextFont(42);
        Text_mathc_pfRICH->AddText(Form("MC ep DIS %ix%i GeV", e_energy, p_energy));
        Text_mathc_pfRICH->AddText(Form("K pfRICH purity: %0.1f%%", K_purity));
        if( mom_bin < nMomBins )
        {
          Text_mathc_pfRICH->AddText(Form("%0.1f < p < %0.1f GeV/c", mom_bins[mom_bin], mom_bins[mom_bin+1]));
          Text_mathc_pfRICH->AddText(Form("%0.1f < Q^{2} < %0.1f GeV^{2}/c^{2}", Q2_bins[Q2bin], Q2_bins[Q2bin+1]));
          Text_mathc_pfRICH->AddText(Form("%0.3f < y < %0.3f", y_bins[y_bin], y_bins[y_bin+1]));
        }
        Text_mathc_pfRICH->SetFillColorAlpha(0, 0.01);
        Text_mathc_pfRICH->Draw("same");

        K_MC_RC_match_pfRICH_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/EIC/EICrecon/figs/%ix%i/K_purity/K_purity_pfRICH_mom_%i_Q2_%i_yBin_%i.png", e_energy, p_energy, mom_bin , Q2bin, y_bin));
      }
    }
  }



  //fitPiTOF->Write();
  //fitKTOF->Write();

  inFile->Close();

  return;

}
