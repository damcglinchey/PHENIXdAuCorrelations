////////////////////////////////////////////////////////////////////////////////
//
// Plot the d+Au correlation functions
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 1 Jul 2016
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

using namespace std;

void plot_corrfuncs()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  //==========================================================================//
  // SET RUNNING CONDITIONS
  //==========================================================================//

  //-- dAu 200
  const int NC  =  6; // number of centrality bins
  const int NZ  = 10; // number of z bins
  const int NPT =  7; // number of pT bins
  // const char *inFile = "rootfiles/CorrFuncdAuBES_dAu200.root";
  // int energy = 200;
  const char *inFile = "rootfiles/CorrFuncdAuBES_dAu62.root";
  int energy = 62;
  // const char *inFile = "rootfiles/CorrFuncdAuBES_dAu39.root";
  // int energy = 39;

  float ptl[] = {0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00};
  float pth[] = {0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 5.00};
  int cl[] = {0,  5, 10, 20, 40,  60};
  int ch[] = {5, 10, 20, 40, 60, 100};
  // float ptl[] = {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00};
  // float pth[] = {0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 5.00};

  bool printPlots = true;

  // correlations between detectors
  const int NCORRPT = 4; // number of correlations with pT binning
  const int NCORR = 8; // number of correlations
  const char *cCorr[] =
  {
    "CNTBBCS",
    "CNTBBCN",
    "CNTFVTXS",
    "CNTFVTXN",
    "BBCNBBCS",
    "BBCNFVTXS",
    "FVTXNBBCS",
    "FVTXNFVTXS",
  };


  // Rebin value
  const int REBIN = 6;

  const int NPAR = 5; // number of orders in the fourier fit
  int fitColor[NPAR] = {kBlue, kRed, kGreen + 2, kMagenta + 2, kYellow + 2};

  //==========================================================================//
  // DECLARE VARIABLES
  //==========================================================================//

  //-- histograms
  // pT dependent
  TH1D* dphi_FG12_pt[NCORRPT][NZ][NC][NPT];
  TH1D* dphi_BG12_pt[NCORRPT][NZ][NC][NPT];

  TH1D* dphi_FGsum_pt[NCORRPT][NC][NPT];
  TH1D* dphi_BGsum_pt[NCORRPT][NC][NPT];
  TH1D* dphi_corr_pt[NCORRPT][NC][NPT];

  // pT integrated
  TH1D* dphi_FG12[NCORR][NZ][NC];
  TH1D* dphi_BG12[NCORR][NZ][NC];

  TH1D* dphi_FGsum[NCORR][NC];
  TH1D* dphi_BGsum[NCORR][NC];
  TH1D* dphi_corr[NCORR][NC];

  //-- fitting

  TF1* fcorr = new TF1("fcorr",
                       "1 + 2*[0]*TMath::Cos(x) + 2*[1]*TMath::Cos(2*x) + 2*[2]*TMath::Cos(3*x) + 2*[3]*TMath::Cos(4*x) + 2*[4]*TMath::Cos(5*x)",
                       -1 * TMath::Pi() / 2., 3 * TMath::Pi() / 2.);
  fcorr->SetLineColor(kBlack);
  fcorr->SetLineStyle(2);

  // pT dependent cn's
  double cn_pt[NCORRPT][NC][NPT][NPAR];
  double cn_pt_e[NCORRPT][NC][NPT][NPAR];
  // pT integrated cn's
  double cn[NCORR][NC][NPAR];
  double cn_e[NCORR][NC][NPAR];


  //-- cn graphs
  TGraphErrors *gcn_pt[NCORRPT][NC][NPAR];
  TGraphErrors *gcn[NCORR][NPAR];

  //-- temp stuff
  char hname[500];


  //==========================================================================//
  // GET HISTOGRAMS
  //==========================================================================//
  cout << endl;
  cout << "--> Reading correlation histograms from " << inFile << endl;

  TFile *fin = TFile::Open(inFile);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << inFile << endl;
    return;
  }

  //-- Loop over all z vrtx directories
  for (int iz = 0; iz < NZ; iz++)
  {
    //-- change to the correct directory
    bool chdir = fin->cd(Form("c00_z%02i_r00", iz));
    if (!chdir)
    {
      cout << "ERROR!! Unable to change into directory for "
           << " iz: " << iz << endl;
      return;
    }

    //-- Loop over all the centrality bins
    for (int ic = 0; ic < NC; ic++)
    {

      //-- Get the pT dependent CNT-FVTX and CNT-BBC correlations
      for (int icorr = 0; icorr < NCORRPT; icorr++)
      {

        for (int ipt = 0; ipt < NPT; ipt++)
        {

          // FG12
          sprintf(hname,
                  "dphi_%s_%03i_%03i_%03.0f_%03.0f_FG12",
                  cCorr[icorr],
                  cl[ic], ch[ic],
                  ptl[ipt] * 100., pth[ipt] * 100);

          dphi_FG12_pt[icorr][iz][ic][ipt] = (TH1D*) gROOT->FindObject(hname);
          if (!dphi_FG12_pt[icorr][iz][ic][ipt])
          {
            cout << "ERROR!! Unable to find " << hname << " for "
                 << " ic:" << ic << " iz:" << iz << " ipt:" << ipt
                 << endl;
            return;
          }
          dphi_FG12_pt[icorr][iz][ic][ipt]->SetDirectory(0);
          sprintf(hname,
                  "dphi_FG12_pt_%s_%i_%i_%i", cCorr[icorr], iz, ic, ipt);
          dphi_FG12_pt[icorr][iz][ic][ipt]->SetName(hname);
          dphi_FG12_pt[icorr][iz][ic][ipt]->SetLineColor(kBlue);
          dphi_FG12_pt[icorr][iz][ic][ipt]->SetLineWidth(2);
          dphi_FG12_pt[icorr][iz][ic][ipt]->Rebin(REBIN);


          // BG12
          sprintf(hname,
                  "dphi_%s_%03i_%03i_%03.0f_%03.0f_BG12",
                  cCorr[icorr],
                  cl[ic], ch[ic],
                  ptl[ipt] * 100., pth[ipt] * 100);

          dphi_BG12_pt[icorr][iz][ic][ipt] = (TH1D*) gROOT->FindObject(hname);
          if (!dphi_BG12_pt[icorr][iz][ic][ipt])
          {
            cout << "ERROR!! Unable to find " << hname << " for "
                 << " ic:" << ic << " iz:" << iz << " ipt:" << ipt
                 << endl;
            return;
          }
          dphi_BG12_pt[icorr][iz][ic][ipt]->SetDirectory(0);
          sprintf(hname,
                  "dphi_BG12_pt_%s_%i_%i_%i", cCorr[icorr], iz, ic, ipt);
          dphi_BG12_pt[icorr][iz][ic][ipt]->SetName(hname);
          dphi_BG12_pt[icorr][iz][ic][ipt]->SetLineColor(kRed);
          dphi_BG12_pt[icorr][iz][ic][ipt]->Rebin(REBIN);


          //-- sum the FG & BG (pT dependent)
          if (iz == 0)
          {
            // pT dependent sums
            dphi_FGsum_pt[icorr][ic][ipt] =
              (TH1D*) dphi_FG12_pt[icorr][iz][ic][ipt]->Clone(
                Form("dphi_FGsum_pt_%s_%i_%i", cCorr[icorr], ic, ipt));
            dphi_FGsum_pt[icorr][ic][ipt]->SetDirectory(0);

            dphi_BGsum_pt[icorr][ic][ipt] =
              (TH1D*) dphi_BG12_pt[icorr][iz][ic][ipt]->Clone(
                Form("dphi_BGsum_pt_%s_%i_%i", cCorr[icorr], ic, ipt));
            dphi_BGsum_pt[icorr][ic][ipt]->SetDirectory(0);

          }
          else
          {
            dphi_FGsum_pt[icorr][ic][ipt]->Add(dphi_FG12_pt[icorr][iz][ic][ipt]);
            dphi_BGsum_pt[icorr][ic][ipt]->Add(dphi_BG12_pt[icorr][iz][ic][ipt]);
          }


          //-- sum the FG & BG (pT integrated)
          if (iz == 0 && ipt == 0)
          {
            // pT dependent sums
            dphi_FGsum[icorr][ic] =
              (TH1D*) dphi_FG12_pt[icorr][iz][ic][ipt]->Clone(
                Form("dphi_FGsum_%s_%i", cCorr[icorr], ic));
            dphi_FGsum[icorr][ic]->SetDirectory(0);

            dphi_BGsum[icorr][ic] =
              (TH1D*) dphi_BG12_pt[icorr][iz][ic][ipt]->Clone(
                Form("dphi_BGsum_%s_%i", cCorr[icorr], ic));
            dphi_BGsum[icorr][ic]->SetDirectory(0);

          }
          else
          {
            dphi_FGsum[icorr][ic]->Add(dphi_FG12_pt[icorr][iz][ic][ipt]);
            dphi_BGsum[icorr][ic]->Add(dphi_BG12_pt[icorr][iz][ic][ipt]);
          }


        } // ipt
      } // icorr



      //-- Get the pT integrated correlations
      for (int icorr = NCORRPT; icorr < NCORR; icorr++)
      {

        sprintf(hname, "dphi_%s_%03i_%03i_FG12", cCorr[icorr], cl[ic], ch[ic]);
        dphi_FG12[icorr][iz][ic] = (TH1D*) gROOT->FindObject(hname);
        if (!dphi_FG12[icorr][iz][ic])
        {
          cout << "ERROR!! Unable to find " << hname << " for "
               << " icorr:" << icorr << " ic:" << ic << " iz:" << iz
               << endl;
          return;
        }
        dphi_FG12[icorr][iz][ic]->SetDirectory(0);
        sprintf(hname, "dphi_FG12_%s_%i_%i", cCorr[icorr], ic, iz);
        dphi_FG12[icorr][iz][ic]->SetName(hname);
        dphi_FG12[icorr][iz][ic]->SetLineColor(kBlue);
        dphi_FG12[icorr][iz][ic]->SetLineWidth(2);
        dphi_FG12[icorr][iz][ic]->Rebin(REBIN);


        sprintf(hname, "dphi_%s_%03i_%03i_BG12", cCorr[icorr], cl[ic], ch[ic]);
        dphi_BG12[icorr][iz][ic] = (TH1D*) gROOT->FindObject(hname);
        if (!dphi_BG12[icorr][iz][ic])
        {
          cout << "ERROR!! Unable to find " << hname << " for "
               << " icorr:" << icorr << " ic:" << ic << " iz:" << iz
               << endl;
          return;
        }
        dphi_BG12[icorr][iz][ic]->SetDirectory(0);
        sprintf(hname, "dphi_BG12_%s_%i_%i", cCorr[icorr], ic, iz);
        dphi_BG12[icorr][iz][ic]->SetName(hname);
        dphi_BG12[icorr][iz][ic]->SetLineColor(kRed);
        dphi_BG12[icorr][iz][ic]->Rebin(REBIN);


        // Sum the FG & BG
        if (iz == 0)
        {
          dphi_FGsum[icorr][ic] =
            (TH1D*) dphi_FG12[icorr][iz][ic]->Clone(
              Form("dphi_FGsum_%s_%i", cCorr[icorr], ic));
          dphi_FGsum[icorr][ic]->SetDirectory(0);

          dphi_BGsum[icorr][ic] =
            (TH1D*) dphi_BG12[icorr][iz][ic]->Clone(
              Form("dphi_BGsum_%s_%i", cCorr[icorr], ic));
          dphi_BGsum[icorr][ic]->SetDirectory(0);
        }
        else
        {
          dphi_FGsum[icorr][ic]->Add(dphi_FG12[icorr][iz][ic]);
          dphi_BGsum[icorr][ic]->Add(dphi_BG12[icorr][iz][ic]);
        }

      } // icorr

    } // ic

    //-- change back to the main directory
    fin->cd();

  } // iz

  fin->Close();
  delete fin;

  //==========================================================================//
  // CALCULATE THE CORRELATION FUNCTIONS
  //==========================================================================//
  cout << endl;
  cout << "--> Calculating the correlation functions" << endl;

  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    for (int ic = 0; ic < NC; ic++)
    {
      // pT dependent (icorr < NCORRPT)
      if ( icorr < NCORRPT )
      {
        for (int ipt = 0; ipt < NPT; ipt++)
        {
          cout << " icorr:" << icorr << " ic:" << ic << " ipt:" << ipt << endl;

          dphi_corr_pt[icorr][ic][ipt] =
            (TH1D*) dphi_FGsum_pt[icorr][ic][ipt]->Clone(
              Form("dphi_corr_pt_%i_%i_%i", icorr, ic, ipt));
          dphi_corr_pt[icorr][ic][ipt]->SetDirectory(0);
          dphi_corr_pt[icorr][ic][ipt]->Divide(dphi_BGsum_pt[icorr][ic][ipt]);
          dphi_corr_pt[icorr][ic][ipt]->Scale(
            dphi_BGsum_pt[icorr][ic][ipt]->Integral() / dphi_FGsum_pt[icorr][ic][ipt]->Integral());
          dphi_corr_pt[icorr][ic][ipt]->SetMarkerStyle(kOpenCircle);
          dphi_corr_pt[icorr][ic][ipt]->SetMarkerColor(kBlack);
        } // ipt
      } // if icorr < NCORRPT


      // pT integrated
      cout << " icorr:" << icorr << " ic:" << ic << endl;
      dphi_corr[icorr][ic] =
        (TH1D*) dphi_FGsum[icorr][ic]->Clone(
          Form("dphi_corr_%i_%i", icorr, ic));
      dphi_corr[icorr][ic]->Sumw2();
      dphi_corr[icorr][ic]->SetDirectory(0);
      dphi_corr[icorr][ic]->Divide(dphi_BGsum[icorr][ic]);
      dphi_corr[icorr][ic]->Scale(
        dphi_BGsum[icorr][ic]->Integral() / dphi_FGsum[icorr][ic]->Integral());
      dphi_corr[icorr][ic]->SetMarkerStyle(kOpenCircle);
      dphi_corr[icorr][ic]->SetMarkerColor(kBlack);
    } // ic

  } // icorr


  //==========================================================================//
  // MAKE GRAPHS OF CORRELATION COEFFICIENTS
  //==========================================================================//
  cout << endl;
  cout << "--> Making Graphs of correlation coefficients" << endl;



  //==========================================================================//
  // FIT CORRELATIONS
  //==========================================================================//
  cout << endl;
  cout << "--> Fitting correlation functions" << endl;


  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    for (int ic = 0; ic < NC; ic++)
    {
      //-- pT dependent
      if (icorr < NCORRPT)
      {
        for (int ipt = 0; ipt < NPT; ipt++)
        {
          dphi_corr_pt[icorr][ic][ipt]->Fit(fcorr, "RQ0N");
          for (int i = 0; i < NPAR; i++)
          {
            cn_pt[icorr][ic][ipt][i] = fcorr->GetParameter(i);
            cn_pt_e[icorr][ic][ipt][i] = fcorr->GetParError(i);
          }
        } // ipt
      }

      //-- pT integrated

      dphi_corr[icorr][ic]->Fit(fcorr, "RQ0N");
      for (int i = 0; i < NPAR; i++)
      {
        cn[icorr][ic][i] = fcorr->GetParameter(i);
        cn_e[icorr][ic][i] = fcorr->GetParError(i);
      }

    } // ic
  } // icorr



  //==========================================================================//
  // PLOT OBJECTS
  //==========================================================================//
  cout << endl;
  cout << "--> plotting" << endl;

  char ctitle[500];

  TLatex ltitle;
  ltitle.SetTextAlign(22);
  ltitle.SetNDC();

  TLine l1;
  l1.SetLineStyle(2);

  TLatex le;
  le.SetNDC();

  TLatex lc;
  lc.SetNDC();
  lc.SetTextAlign(31);

  //==========================================================================//
  // PLOT
  //==========================================================================//

  TCanvas *cFGBGpt[NCORRPT][NPT];
  TCanvas *ccorrpt[NCORRPT][NPT];
  TCanvas *cFGBG[NCORR];
  TCanvas *ccorr[NCORR];
  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    //-- pT dependent
    if (icorr < NCORRPT)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
      {
        //FGBG
        sprintf(hname, "cFGBGpt_%i_%i", icorr, ipt);
        cFGBGpt[icorr][ipt] = new TCanvas(hname, hname, 1200, 1000);
        cFGBGpt[icorr][ipt]->Divide(3, 2, 0, 0);
        for (int ic = 0; ic < NC; ic++)
        {
          cFGBGpt[icorr][ipt]->GetPad(ic + 1)->SetTopMargin(0.1);
          cFGBGpt[icorr][ipt]->GetPad(ic + 1)->SetRightMargin(0.02);
          cFGBGpt[icorr][ipt]->GetPad(ic + 1)->SetBottomMargin(0.1);
          cFGBGpt[icorr][ipt]->GetPad(ic + 1)->SetLeftMargin(0.1);
          cFGBGpt[icorr][ipt]->GetPad(ic + 1)->SetTicks(1, 1);

          cFGBGpt[icorr][ipt]->cd(ic + 1);
          dphi_FGsum_pt[icorr][ic][ipt]->Scale(
            1. / dphi_FGsum_pt[icorr][ic][ipt]->Integral());
          dphi_BGsum_pt[icorr][ic][ipt]->Scale(
            1. / dphi_BGsum_pt[icorr][ic][ipt]->Integral());

          dphi_FGsum_pt[icorr][ic][ipt]->Draw();
          dphi_BGsum_pt[icorr][ic][ipt]->Draw("same");

          sprintf(ctitle, "%s FGBG %.2f<pT<%.2f %i-%i%%",
                  cCorr[icorr],
                  ptl[ipt], pth[ipt], cl[ic], ch[ic]);
          ltitle.DrawLatex(0.5, 0.95, ctitle);

          if (ic == 0)
            le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
        } // ic

        //correlation
        sprintf(hname, "ccorrpt_%i_%i", icorr, ipt);
        ccorrpt[icorr][ipt] = new TCanvas(hname, hname, 1200, 1000);
        ccorrpt[icorr][ipt]->Divide(3, 2, 0, 0);
        for (int ic = 0; ic < NC; ic++)
        {
          ccorrpt[icorr][ipt]->GetPad(ic + 1)->SetTopMargin(0.1);
          ccorrpt[icorr][ipt]->GetPad(ic + 1)->SetRightMargin(0.02);
          ccorrpt[icorr][ipt]->GetPad(ic + 1)->SetBottomMargin(0.1);
          ccorrpt[icorr][ipt]->GetPad(ic + 1)->SetLeftMargin(0.1);
          ccorrpt[icorr][ipt]->GetPad(ic + 1)->SetTicks(1, 1);

          ccorrpt[icorr][ipt]->cd(ic + 1);

          //set the scale
          for (int i = 0; i < NPAR; i++)
          {
            if (i == 0)
              fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i]);
            else
              fcorr->SetParameter(i, 0);
          }
          double min = fcorr->GetMinimum(-1, 1);
          double max = dphi_corr_pt[icorr][ic][ipt]->GetMaximum();

          dphi_corr_pt[icorr][ic][ipt]->GetYaxis()->SetRangeUser(
            (min - 0.1 * (1 - min)),
            (max + 0.1 * (max - 1)));


          // Draw
          dphi_corr_pt[icorr][ic][ipt]->Draw("P");

          sprintf(ctitle, "%s %.2f<pT<%.2f %i-%i%%",
                  cCorr[icorr],
                  ptl[ipt], pth[ipt],
                  cl[ic], ch[ic]);
          ltitle.DrawLatex(0.5, 0.95, ctitle);

          l1.DrawLine(-1 * TMath::Pi() / 2., 1, 3 * TMath::Pi() / 2., 1.);

          if (ic == 0)
            le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));

          // plot fit
          fcorr->SetLineColor(kBlack);
          for (int i = 0; i < NPAR; i++)
            fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i]);
          fcorr->DrawCopy("same");

          for (int ipar = 0; ipar < NPAR; ipar++)
          {
            fcorr->SetLineColor(fitColor[ipar]);
            for (int i = 0; i < NPAR; i++)
            {
              if (i == ipar)
                fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i]);
              else
                fcorr->SetParameter(i, 0);
            }
            fcorr->DrawCopy("same");
          }

        } // ic

      } // ipt
    }

    //-- pT integrated
    //FGBG
    sprintf(hname, "cFGBG_%i", icorr);
    cFGBG[icorr] = new TCanvas(hname, hname, 1200, 1000);
    cFGBG[icorr]->Divide(3, 2, 0, 0);
    for (int ic = 0; ic < NC; ic++)
    {
      cFGBG[icorr]->GetPad(ic + 1)->SetTopMargin(0.1);
      cFGBG[icorr]->GetPad(ic + 1)->SetRightMargin(0.02);
      cFGBG[icorr]->GetPad(ic + 1)->SetBottomMargin(0.1);
      cFGBG[icorr]->GetPad(ic + 1)->SetLeftMargin(0.1);
      cFGBG[icorr]->GetPad(ic + 1)->SetTicks(1, 1);

      cFGBG[icorr]->cd(ic + 1);
      dphi_FGsum[icorr][ic]->Scale(
        1. / dphi_FGsum[icorr][ic]->Integral());
      dphi_BGsum[icorr][ic]->Scale(
        1. / dphi_BGsum[icorr][ic]->Integral());

      dphi_FGsum[icorr][ic]->Draw();
      dphi_BGsum[icorr][ic]->Draw("same");

      sprintf(ctitle, "%s FGBG %i-%i%%",
              cCorr[icorr],
              cl[ic], ch[ic]);
      ltitle.DrawLatex(0.5, 0.95, ctitle);

      if (ic == 0)
        le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
    } // ic

    //correlation
    sprintf(hname, "ccorr_%i", icorr);
    ccorr[icorr] = new TCanvas(hname, hname, 1200, 1000);
    ccorr[icorr]->Divide(3, 2, 0, 0);
    for (int ic = 0; ic < NC; ic++)
    {
      ccorr[icorr]->GetPad(ic + 1)->SetTopMargin(0.1);
      ccorr[icorr]->GetPad(ic + 1)->SetRightMargin(0.02);
      ccorr[icorr]->GetPad(ic + 1)->SetBottomMargin(0.1);
      ccorr[icorr]->GetPad(ic + 1)->SetLeftMargin(0.1);
      ccorr[icorr]->GetPad(ic + 1)->SetTicks(1, 1);

      ccorr[icorr]->cd(ic + 1);

      //set the scale
      for (int i = 0; i < NPAR; i++)
      {
        if (i == 0)
          fcorr->SetParameter(i, cn[icorr][ic][i]);
        else
          fcorr->SetParameter(i, 0);
      }
      double min = fcorr->GetMinimum(-1, 1);
      double max = dphi_corr[icorr][ic]->GetMaximum();

      dphi_corr[icorr][ic]->GetYaxis()->SetRangeUser(
        (min - 0.1 * (1 - min)),
        (max + 0.1 * (max - 1)));


      // Draw
      dphi_corr[icorr][ic]->Draw("P");

      sprintf(ctitle, "%s %i-%i%%",
              cCorr[icorr],
              cl[ic], ch[ic]);
      ltitle.DrawLatex(0.5, 0.95, ctitle);

      l1.DrawLine(-1 * TMath::Pi() / 2., 1, 3 * TMath::Pi() / 2., 1.);

      if (ic == 0)
        le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));

      // plot fit
      fcorr->SetLineColor(kBlack);
      for (int i = 0; i < NPAR; i++)
        fcorr->SetParameter(i, cn[icorr][ic][i]);
      fcorr->DrawCopy("same");

      for (int ipar = 0; ipar < NPAR; ipar++)
      {
        fcorr->SetLineColor(fitColor[ipar]);
        for (int i = 0; i < NPAR; i++)
        {
          if (i == ipar)
            fcorr->SetParameter(i, cn[icorr][ic][i]);
          else
            fcorr->SetParameter(i, 0);
        }
        fcorr->DrawCopy("same");
      }

      lc.SetTextColor(fitColor[0]);
      lc.DrawLatex(0.45, 0.7, Form("C_{1} = % .4f",
                                   cn[icorr][ic][0]));
      lc.SetTextColor(fitColor[1]);
      lc.DrawLatex(0.45, 0.64, Form("C_{2} = % .4f",
                                   cn[icorr][ic][1]));
      lc.SetTextColor(kBlack);
      lc.DrawLatex(0.45, 0.58, Form("C_{2}/C_{1} = % .4f",
                                   cn[icorr][ic][1] / cn[icorr][ic][0]));


    } // ic

  } // icorr




  //==========================================================================//
  // PRINT PLOTS
  //==========================================================================//
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing plots" << endl;

    char cname[500];

    for (int icorr = 0; icorr < NCORR; icorr++)
    {
      //-- pT dependent
      if (icorr < NCORRPT)
      {
        for (int ipt = 0; ipt < NPT; ipt++)
        {
          sprintf(cname, "pdfs/dAu%i_%s_FGBG_pt%i.pdf", energy, cCorr[icorr], ipt);
          cFGBGpt[icorr][ipt]->Print(cname);

          sprintf(cname, "pdfs/dAu%i_%s_corr_pt%i.pdf", energy, cCorr[icorr], ipt);
          ccorrpt[icorr][ipt]->Print(cname);
        }
      }

      //-- pT integrated
      sprintf(cname, "pdfs/dAu%i_%s_FGBG.pdf", energy, cCorr[icorr]);
      cFGBG[icorr]->Print(cname);

      sprintf(cname, "pdfs/dAu%i_%s_corr.pdf", energy, cCorr[icorr]);
      ccorr[icorr]->Print(cname);

    } // icorr


  }

}