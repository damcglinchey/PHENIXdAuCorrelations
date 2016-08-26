////////////////////////////////////////////////////////////////////////////////
//
// Plot the diagnostic histograms output by CorrFuncdAuBESUltraLightReco
// for each Run 16 d+Au energy
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 25 Aug 2016
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>

#include <iostream>
#include <iomanip>

using namespace std;

void plot_energy_histos()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //============================================================================
  // SET RUNNING CONDITIONS
  //============================================================================

  const int NE = 4;
  int energy[] = {200, 62, 39, 20};

  bool printPlots = true;

  const int NC  =  6; // number of centrality bins
  int cl[] = {0,  5, 10, 20, 40,  60};
  int ch[] = {5, 10, 20, 40, 60, 100};
  int centColor[] = {kBlack, kBlue, kRed, kGreen + 2, kMagenta + 2, kOrange + 2};

  int energyMarker[] = {kFullCircle, kFullSquare, kFullDiamond, kFullStar};
  int energyColor[] =  {kBlue, kRed, kGreen+2, kMagenta+2};



  //============================================================================
  // DECLARE VARIABLES
  //============================================================================

  //-- Histograms
  TH1D* hcentrality[NE];
  TH1D* hQbbcs[NE];
  TH1D* hzvrtx[NE];

  TH1D* hncnt[NE];
  TH1D* hnfvtxclus[NE][2];
  TH1D* hnfvtxclus_tot[NE];
  TH1D* hnfvtxtrk[NE][2];
  TH1D* hnfvtxtrk_tot[NE];

  TH2D* hnpc1centrality[NE];
  TH1D* hnpc1_cent[NE][NC];

  //-- multiplicity
  double mean_npc1[NE][NC];
  double mean_npc1_e[NE][NC];
  TGraphErrors *gnpc1_cent[NE];

  //-- other
  char fname[500];
  TFile* fin;

  double centidx[NC];
  double centidx_e[NC] = {0};
  for (int ic = 0; ic < NC; ic++)
    centidx[ic] = ic + 0.5;

  //============================================================================
  // GET HISTOGRAMS
  //============================================================================
  cout << endl;
  cout << "--> Getting Histograms" << endl;

  for (int ie = 0; ie < NE; ie++)
  {
    sprintf(fname, "rootfiles/Histos_CorrFuncdAuBES_dAu%i.root", energy[ie]);
    fin = TFile::Open(fname);
    if (!fin)
    {
      cout << "ERROR!! Unable to open " << fname << endl;
      return;
    }
    cout << "--> For dAu" << energy[ie] << " Getting histograms from:"
         << fname << endl;


    hcentrality[ie] = (TH1D*) fin->Get("hcentrality");
    if (!hcentrality[ie])
    {
      cout << "ERROR!! Unable to find hcentrality in " << fname << endl;
      return;
    }
    hcentrality[ie]->SetName(Form("hcentrality_%i", ie));
    hcentrality[ie]->SetDirectory(0);
    hcentrality[ie]->SetLineColor(energyColor[ie]);

    hQbbcs[ie] = (TH1D*) fin->Get("hQbbcs");
    hQbbcs[ie]->SetName(Form("hQbbcs_%i", ie));
    hQbbcs[ie]->SetDirectory(0);
    hQbbcs[ie]->SetLineColor(energyColor[ie]);

    hzvrtx[ie] = (TH1D*) fin->Get("hzvrtx");
    hzvrtx[ie]->SetName(Form("hzvrtx_%i", ie));
    hzvrtx[ie]->SetDirectory(0);
    hzvrtx[ie]->SetLineColor(energyColor[ie]);

    hncnt[ie] = (TH1D*) fin->Get("hncnt");
    hncnt[ie]->SetName(Form("hncnt_%i", ie));
    hncnt[ie]->SetDirectory(0);
    hncnt[ie]->SetLineColor(energyColor[ie]);

    hnfvtxclus_tot[ie] = (TH1D*) fin->Get("hnfvtxclus_tot");
    hnfvtxclus_tot[ie]->SetName(Form("hnfvtxclus_tot_%i", ie));
    hnfvtxclus_tot[ie]->SetDirectory(0);
    hnfvtxclus_tot[ie]->SetLineColor(energyColor[ie]);

    hnfvtxtrk_tot[ie] = (TH1D*) fin->Get("hnfvtxtrk_tot");
    hnfvtxtrk_tot[ie]->SetName(Form("hnfvtxtrk_tot_%i", ie));
    hnfvtxtrk_tot[ie]->SetDirectory(0);
    hnfvtxtrk_tot[ie]->SetLineColor(energyColor[ie]);

    // Get N&S FVTX multiplicity
    for (int j = 0; j < 2; j++)
    {
      hnfvtxclus[ie][j] = (TH1D*) fin->Get(Form("hnfvtxclus_%i", j));
      hnfvtxclus[ie][j]->SetName(Form("hnfvtxclus_%i_%i", ie, j));
      hnfvtxclus[ie][j]->SetDirectory(0);
      hnfvtxclus[ie][j]->SetLineColor(energyColor[ie]);

      hnfvtxtrk[ie][j] = (TH1D*) fin->Get(Form("hnfvtxtrk_%i", j));
      hnfvtxtrk[ie][j]->SetName(Form("hnfvtxtrk_%i_%i", ie, j));
      hnfvtxtrk[ie][j]->SetDirectory(0);
      hnfvtxtrk[ie][j]->SetLineColor(energyColor[ie]);
    }


    //-- Get PC1 multiplicity in each centrality bin
    hnpc1centrality[ie] = (TH2D*) fin->Get("hnpc1_centrality");
    hnpc1centrality[ie]->SetName(Form("hnpc1centrality_%i", ie));
    hnpc1centrality[ie]->SetDirectory(0);
    hnpc1centrality[ie]->SetLineColor(energyColor[ie]);

    for (int ic = 0; ic < NC; ic++)
    {
      int bl = hnpc1centrality[ie]->GetXaxis()->FindBin(cl[ic]);
      int bh = hnpc1centrality[ie]->GetXaxis()->FindBin(ch[ic] - 1);

      hnpc1_cent[ie][ic] = (TH1D*) hnpc1centrality[ie]->ProjectionY(
                             Form("hnpc1_cent_%i_%i", ie, ic),
                             bl, bh);
      hnpc1_cent[ie][ic]->SetDirectory(0);
      hnpc1_cent[ie][ic]->SetLineColor(centColor[ic]);
      hnpc1_cent[ie][ic]->Scale(1. / ((double)(ch[ic] - cl[ic])));

      mean_npc1[ie][ic] = hnpc1_cent[ie][ic]->GetMean();
      mean_npc1_e[ie][ic] = hnpc1_cent[ie][ic]->GetRMS();

    }

    // gnpc1_cent[ie] = new TGraphErrors(NC, centidx, mean_npc1[ie],
    //                                   centidx_e, mean_npc1_e[ie]);
    gnpc1_cent[ie] = new TGraphErrors();
    for (int ic = 0; ic < NC; ic++)
    {
      gnpc1_cent[ie]->SetPoint(ic, centidx[ic] + (ie - 1) * 0.05,
                               mean_npc1[ie][ic]);
      gnpc1_cent[ie]->SetPointError(ic, 0, mean_npc1_e[ie][ic]);
    }
    gnpc1_cent[ie]->SetLineColor(energyColor[ie]);
    gnpc1_cent[ie]->SetMarkerStyle(energyMarker[ie]);
    gnpc1_cent[ie]->SetMarkerColor(energyColor[ie]);


  }

  //============================================================================
  // PLOT OBJECTS
  //============================================================================
  cout << endl;
  cout << "--> Plotting" << endl;

  TLegend *legE = new TLegend(0.5, 0.5, 0.85, 0.85);
  legE->SetFillStyle(0);
  legE->SetBorderSize(0);
  for (int ie = 0; ie < NE; ie++)
  {
    legE->AddEntry(hcentrality[ie],
                   Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy[ie]),
                   "L");
  }

  TLegend *legEP = new TLegend(0.5, 0.5, 0.85, 0.85);
  legEP->SetFillStyle(0);
  legEP->SetBorderSize(0);
  for (int ie = 0; ie < NE; ie++)
  {
    legEP->AddEntry(gnpc1_cent[ie],
                    Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy[ie]),
                    "P");
  }


  TH1F* haxis_cncent = new TH1F("haxis_cncent", "", NC, 0, NC);
  haxis_cncent->GetYaxis()->CenterTitle();
  haxis_cncent->GetYaxis()->SetTitleSize(0.08);
  haxis_cncent->GetYaxis()->SetTitleOffset(0.75);
  haxis_cncent->GetXaxis()->SetLabelSize(0.08);

  char ctitle[500];
  for (int ibin = 1; ibin <= haxis_cncent->GetNbinsX(); ibin++)
  {
    sprintf(ctitle, "%i - %i%%", cl[ibin - 1], ch[ibin - 1]);
    haxis_cncent->GetXaxis()->SetBinLabel(ibin, ctitle);
  }



  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  //============================================================================
  // PLOT
  //============================================================================

  TCanvas *ccheck = new TCanvas("ccheck", "checking", 1200, 1000);
  ccheck->Divide(3, 2);

  ccheck->cd(1);
  gPad->SetLogy();
  hcentrality[0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hcentrality[ie]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "Centrality");
  legE->Draw("same");

  ccheck->cd(2);
  gPad->SetLogy();
  hQbbcs[0]->GetXaxis()->SetRangeUser(0, 300);
  hQbbcs[0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hQbbcs[ie]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "BBC South Charge");

  ccheck->cd(3);
  hzvrtx[0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hzvrtx[ie]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "Z Vertex");

  ccheck->cd(4);
  gPad->SetLogy();
  hncnt[0]->GetXaxis()->SetRangeUser(0, 20);
  hncnt[0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hncnt[ie]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "# PHCentralTracks");

  ccheck->cd(5);
  gPad->SetLogy();
  hnfvtxclus_tot[0]->GetXaxis()->SetRangeUser(0, 800);
  hnfvtxclus_tot[0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hnfvtxclus_tot[ie]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "Total # FVTX Clusters");

  ccheck->cd(6);
  gPad->SetLogy();
  hnfvtxtrk[0][0]->GetXaxis()->SetRangeUser(0, 50);
  hnfvtxtrk[0][0]->Draw();

  for (int ie = 1; ie < NE; ie++)
    hnfvtxtrk[ie][0]->Draw("same");

  ltitle.DrawLatex(0.5, 0.95, "# South FVTX Tracks");






  TCanvas* cnpc1cent = new TCanvas("cnpc1cent", "npc1cent", 1200, 1000);
  cnpc1cent->Divide(2, 2);

  for (int ie = 0; ie < NE; ie++)
  {
    cnpc1cent->cd(ie + 1);
    gPad->SetLogy();
    hnpc1_cent[ie][0]->GetXaxis()->SetRangeUser(0, 50);
    hnpc1_cent[ie][0]->Draw();
    for (int ic = 1; ic < NC; ic++)
      hnpc1_cent[ie][ic]->Draw("same");

    ltitle.DrawLatex(0.5, 0.95, Form("dAu #sqrt{s_{_{NN}}}=%i GeV", energy[ie]));
  }


  TCanvas *cmult = new TCanvas("cmult", "mult", 800, 800);

  cmult->cd(1);
  haxis_cncent->SetMinimum(0);
  haxis_cncent->SetMaximum(25);
  haxis_cncent->SetTitle(";;<N_{PC1}>");
  haxis_cncent->GetYaxis()->SetTitleSize(0.04);
  haxis_cncent->GetXaxis()->SetLabelSize(0.04);
  haxis_cncent->Draw();

  for (int ie = 0; ie < NE; ie++)
    gnpc1_cent[ie]->Draw("P");

  legEP->Draw("same");

  //============================================================================
  // PRINT PLOTS
  //============================================================================
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing plots" << endl;

    ccheck->SaveAs("pdfs/histos_check.pdf");
    cnpc1cent->SaveAs("pdfs/NPC1_centrality.pdf");
    cmult->SaveAs("pdfs/multiplicity_energy.pdf");


    //-- print latex table
    cout << setprecision(4);
    cout << endl;
    cout << "centrality";
    for (int ie = 0; ie < NE; ie++)
    cout << " & " << setw(4) << energy[ie] << " GeV";
    cout << " \\\\" << endl;
    for (int ic = 0; ic < NC; ic++)
    {
      cout << cl[ic] << "--" << ch[ic] << "\\%";
      for (int ie = 0; ie < NE; ie++)
        cout << " & " << setw(8) << mean_npc1[ie][ic];
      cout << " \\\\" << endl;
    }


  }
}