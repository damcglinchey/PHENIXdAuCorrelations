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

  bool printPlots = false;

  const int NE = 3;
  int energy[] = {200, 62, 39, 20};
  int energyColor[] =  {kBlue, kRed, kGreen+2, kMagenta+2};

  //==========================================================================//
  // DECLARE VARIABLES
  //==========================================================================//

  TGraphErrors *gv2_CNTFVTXSFVTXN_mult[NE];
  TGraphErrors *gv2_CNTBBCSFVTXS_mult[NE];
  TGraphErrors *gv2_CNTBBCSFVTXN_mult[NE];

  TGraphErrors *gv2_CNTFVTXSFVTXN_pT[NE];
  TGraphErrors *gv2_CNTBBCSFVTXS_pT[NE];

  TGraphErrors *gc2c1_deta[NE];


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
    //-- v2 vs multiplicity
    sprintf(hname, "gv2_dAu%i_CNTFVTXSFVTXN_mult", energy[ie]);
    gv2_CNTFVTXSFVTXN_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTFVTXSFVTXN_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gv2_CNTFVTXSFVTXN_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTFVTXSFVTXN_mult[ie]->SetLineColor(energyColor[ie]);

    sprintf(hname, "gv2_dAu%i_CNTBBCSFVTXS_mult", energy[ie]);
    gv2_CNTBBCSFVTXS_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTBBCSFVTXS_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gv2_CNTBBCSFVTXS_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTBBCSFVTXS_mult[ie]->SetLineColor(energyColor[ie]);

    sprintf(hname, "gv2_dAu%i_CNTBBCSFVTXN_mult", energy[ie]);
    gv2_CNTBBCSFVTXN_mult[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTBBCSFVTXN_mult[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gv2_CNTBBCSFVTXN_mult[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTBBCSFVTXN_mult[ie]->SetLineColor(energyColor[ie]);


    //-- v2 vs pt 0-5%
    sprintf(hname, "gv2_dAu%i_CNTFVTXSFVTXN_pT_c0", energy[ie]);
    gv2_CNTFVTXSFVTXN_pT[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTFVTXSFVTXN_pT[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gv2_CNTFVTXSFVTXN_pT[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTFVTXSFVTXN_pT[ie]->SetLineColor(energyColor[ie]);

    sprintf(hname, "gv2_dAu%i_CNTBBCSFVTXS_pT_c0", energy[ie]);
    gv2_CNTBBCSFVTXS_pT[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gv2_CNTBBCSFVTXS_pT[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gv2_CNTBBCSFVTXS_pT[ie]->SetMarkerColor(energyColor[ie]);
    gv2_CNTBBCSFVTXS_pT[ie]->SetLineColor(energyColor[ie]);


    //-- v2 vs delta eta
    sprintf(hname, "gc2c1_deta_dAu%i_c0", energy[ie]);
    gc2c1_deta[ie] = (TGraphErrors*) fin->Get(hname);
    if (!gc2c1_deta[ie])
    {
      cout << "ERROR!! Unable to find " << hname << " in " << inFile << endl;
      return;
    }
    gc2c1_deta[ie]->SetMarkerColor(energyColor[ie]);
    gc2c1_deta[ie]->SetLineColor(energyColor[ie]);
  }

  //==========================================================================//
  // PLOT OBJECTS
  //==========================================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  TH1F* haxis_mult = new TH1F("haxis_mult", ";<N_{PC1}>;v_{2}", 22, 0, 22);
  haxis_mult->SetMinimum(0);
  haxis_mult->SetMaximum(0.20);

  TLegend *leg3sub = new TLegend(0.55, 0.75, 0.88, 0.88);
  leg3sub->SetFillStyle(0);
  leg3sub->SetBorderSize(0);
  leg3sub->AddEntry(gv2_CNTFVTXSFVTXN_mult[0], "3 sub-event CNT-FVTXN-FVTXS", "P");
  leg3sub->AddEntry(gv2_CNTBBCSFVTXN_mult[0], "3 sub-event CNT-FVTXN-BBCS", "P");
  leg3sub->AddEntry(gv2_CNTBBCSFVTXS_mult[0], "3 sub-event CNT-BBCS-FVTXS", "P");


  TLegend *legenergy = new TLegend(0.15, 0.75, 0.4, 0.88);
  legenergy->SetFillStyle(0);
  legenergy->SetBorderSize(0);
  for (int ie = 0; ie < NE; ie++)
    legenergy->AddEntry(gv2_CNTBBCSFVTXS_mult[ie], Form("dAu %i GeV", energy[ie]), "P");


  TLegend *legept = new TLegend(0.15, 0.75, 0.35, 0.88);
  legept->SetFillStyle(0);
  legept->SetBorderSize(0);
  legept->SetTextAlign(22);
  legept->SetHeader("0--5%");
  for (int ie = 0; ie < NE; ie++)
    legept->AddEntry(gv2_CNTBBCSFVTXS_mult[ie], Form("dAu %i GeV", energy[ie]), "P");


  TLegend *legdeta = new TLegend(0.55, 0.75, 0.88, 0.88);
  legdeta->SetFillStyle(0);
  legdeta->SetBorderSize(0);
  for (int ie = 0; ie < NE; ie++)
    legdeta->AddEntry(gc2c1_deta[ie], Form("dAu %i GeV", energy[ie]), "P");

  TH1F* haxis_vnpt = new TH1F("haxis_vnpt", ";p_{T} [GeV/c];v_{2}", 100, 0, 5);
  haxis_vnpt->SetMinimum(0.01);
  haxis_vnpt->SetMaximum(0.29);
  haxis_vnpt->GetYaxis()->CenterTitle();
  haxis_vnpt->GetYaxis()->SetTitleFont(63);
  haxis_vnpt->GetYaxis()->SetTitleSize(20);
  haxis_vnpt->GetYaxis()->SetTitleOffset(1.7);
  haxis_vnpt->GetYaxis()->SetLabelFont(63);
  haxis_vnpt->GetYaxis()->SetLabelSize(14);
  haxis_vnpt->GetXaxis()->SetTitleFont(63);
  haxis_vnpt->GetXaxis()->SetTitleSize(20);
  haxis_vnpt->GetXaxis()->SetTitleOffset(1.7);
  haxis_vnpt->GetXaxis()->SetLabelFont(63);
  haxis_vnpt->GetXaxis()->SetLabelSize(14);




  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  TLatex lcent;
  lcent.SetNDC();

  //==========================================================================//
  // PLOT
  //==========================================================================//

  TCanvas *cv2mult = new TCanvas("cv2mult", "v2 mult", 800, 600);

  cv2mult->cd(1);
  haxis_mult->Draw();

  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTFVTXSFVTXN_mult[ie]->DrawClone("P");
    // gv2_CNTBBCSFVTXS_mult[ie]->DrawClone("P");
    // gv2_CNTBBCSFVTXN_mult[ie]->DrawClone("P");
  }

  // leg3sub->Draw("same");
  ltitle.DrawLatex(0.5, 0.95, "CNT--FVTXS--FVTXN");
  legenergy->Draw("same");


  TCanvas* cv2energy = new TCanvas("cv2energy", "v2 energy", 600, 900);
  cv2energy->Divide(1, 3);

  cv2energy->cd(1);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTFVTXSFVTXN_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--FVTXS--FVTXN");
  legenergy->Draw("same");

  cv2energy->cd(2);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTBBCSFVTXN_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--BBCS--FVTXN");

  cv2energy->cd(3);
  haxis_mult->Draw();
  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTBBCSFVTXS_mult[ie]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, "CNT--BBCS--FVTXS");



  TCanvas* cv2pt = new TCanvas("cv2pt", "v2 pT", 800, 800);
  cv2pt->cd(1);
  haxis_vnpt->GetXaxis()->SetRangeUser(0, 3);
  haxis_vnpt->Draw();

  for (int ie = 0; ie < NE; ie++)
  {
    gv2_CNTFVTXSFVTXN_pT[ie]->Draw("P");
  }

  legept->Draw("same");
  ltitle.DrawLatex(0.5, 0.95, "CNT--FVTXS--FVTXN");




  TCanvas* cdeta = new TCanvas("cdeta", "deta", 800, 800);

  cdeta->cd(1);
  gc2c1_deta[0]->GetYaxis()->SetRangeUser(0, 0.7);
  gc2c1_deta[0]->Draw("AP");
  for (int ie = 0; ie < NE; ie++)
    gc2c1_deta[ie]->Draw("P");

  legdeta->Draw("same");

  //==========================================================================//
  // Print Plots
  //==========================================================================//
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing Plots" << endl;

    cv2mult->Print("pdfs/v2_multiplicity.pdf");
    cv2energy->Print("pdfs/v2_multiplicity_method.pdf");
    cv2pt->Print("pdfs/v2_pT_energy.pdf");
    cdeta->Print("pdfs/c2c1_deta_energy.pdf");
  }

}