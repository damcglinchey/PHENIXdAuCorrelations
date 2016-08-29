////////////////////////////////////////////////////////////////////////////////
//
// Plot energy dependent C_n and v_n
// Input comes from output of plot_corrfuncs.C
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 29 Aug 2016
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include <iostream>
#include <vector>
#include <utility>

using namespace std;


void plot_energydep()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  //==========================================================================//
  // SET RUNNING CONDITIONS
  //==========================================================================//

  const char* inFile = "correlations.root";

  bool printPlots = true;

  const int NE = 3;
  int energy[] = {200, 62, 39, 20};
  int energyColor[] =  {kBlue, kRed, kGreen+2, kMagenta+2};

  //==========================================================================//
  // DECLARE VARIABLES
  //==========================================================================//

  TGraphErrors *gv2_CNTFVTXSFVTXN_mult[NE];
  TGraphErrors *gv2_CNTBBCSFVTXS_mult[NE];
  TGraphErrors *gv2_CNTBBCSFVTXN_mult[NE];


  char hname[500];

  //==========================================================================//
  // GET RESULTS FROM FILE
  //==========================================================================//
  cout << endl;
  cout << "--> Getting results from " << inFile << endl;

  TFile *fin = TFile::Open(inFile);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << inFile << endl;
    return;
  }

  for (int ie = 0; ie < NE; ie++)
  {
    sprintf(hname, "gv2_dAu%i_CNTFVTXSFVTXN_mult", energy[ie]);
    gv2_CNTFVTXSFVTXN_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTFVTXSFVTXN_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }

    sprintf(hname, "gv2_dAu%i_CNTBBCSFVTXS_mult", energy[ie]);
    gv2_CNTBBCSFVTXS_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTBBCSFVTXS_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }

    sprintf(hname, "gv2_dAu%i_CNTBBCSFVTXN_mult", energy[ie]);
    gv2_CNTBBCSFVTXN_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTBBCSFVTXN_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }

  }

  //==========================================================================//
  // PLOT OBJECTS
  //==========================================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  TH1F* haxis_mult = new TH1F("haxis_mult", ";<N_{PC1}>;v_{2}", 22, 0, 22);
  haxis_mult->SetMinimum(0);
  haxis_mult->SetMaximum(0.30);

  TLegend *leg3sub = new TLegend(0.55, 0.75, 0.88, 0.88);
  leg3sub->SetFillStyle(0);
  leg3sub->SetBorderSize(0);
  leg3sub->AddEntry(gv2_CNTFVTXSFVTXN_mult[0], "3 sub-event CNT-FVTXN-FVTXS", "P");
  leg3sub->AddEntry(gv2_CNTBBCSFVTXN_mult[0], "3 sub-event CNT-FVTXN-BBCS", "P");
  leg3sub->AddEntry(gv2_CNTBBCSFVTXS_mult[0], "3 sub-event CNT-BBCS-FVTXS", "P");


  TLegend *legenergy = new TLegend(0.55, 0.75, 0.88, 0.88);
  legenergy->SetFillStyle(0);
  legenergy->SetBorderSize(0);
  for (int ie = 0; ie < NE; ie++)
    legenergy->AddEntry(gv2_CNTBBCSFVTXS_mult[ie], Form("dAu %i GeV", energy[ie]), "P");



  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  //==========================================================================//
  // PLOT
  //==========================================================================//

  TCanvas *cv2mult = new TCanvas("cv2mult", "v2 mult", 800, 600);

  cv2mult->cd(1);
  haxis_mult->Draw();

  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTFVTXSFVTXN_mult[ie]->DrawClone("P");
    gv2_CNTBBCSFVTXS_mult[ie]->DrawClone("P");
    gv2_CNTBBCSFVTXN_mult[ie]->DrawClone("P");
  }

  leg3sub->Draw("same");



  TCanvas* cv2energy = new TCanvas("cv2energy", "v2 energy", 600, 900);
  cv2energy->Divide(1, 3);

  cv2energy->cd(1);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTFVTXSFVTXN_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTFVTXSFVTXN_mult[ie]->SetLineColor(energyColor[ie]);
    gv2_CNTFVTXSFVTXN_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--FVTXS--FVTXN");
  legenergy->Draw("same");

  cv2energy->cd(2);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTBBCSFVTXN_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTBBCSFVTXN_mult[ie]->SetLineColor(energyColor[ie]);
    gv2_CNTBBCSFVTXN_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--BBCS--FVTXN");

  cv2energy->cd(3);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTBBCSFVTXS_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTBBCSFVTXS_mult[ie]->SetLineColor(energyColor[ie]);
    gv2_CNTBBCSFVTXS_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--BBCS--FVTXS");

  //==========================================================================//
  // Print Plots
  //==========================================================================//
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing Plots" << endl;

    cv2mult->Print("pdfs/v2_multiplicity.pdf");
    cv2energy->Print("pdfs/v2_multiplicity_method.pdf");
  }

}