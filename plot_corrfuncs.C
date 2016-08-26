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
#include <vector>

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
  // const char *inFile = "rootfiles/CorrFuncdAuBES_dAu20.root";
  // int energy = 20;

  float ptl[] = {0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00};
  float pth[] = {0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 5.00};
  int cl[] = {0,  5, 10, 20, 40,  60};
  int ch[] = {5, 10, 20, 40, 60, 100};

  int ptsuml = 3;
  int ptsumh = 5;
  // int ptsuml = 0;
  // int ptsumh = NPT - 1;

  bool printPlots = true;

  // correlations between detectors
  const int NCORRPT = 4; // number of correlations with pT binning
  const int NCORR = 8;   // number of correlations
  enum CORR {
    CNTBBCS = 0,
    CNTBBCN,
    CNTFVTXS,
    CNTFVTXN,
    BBCNBBCS,
    BBCNFVTXS,
    FVTXNBBCS,
    FVTXNFVTXS,
  }; // Note: These better be in the same order as cCorr[]!!!
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
  int corrColor[] =  {kRed, kAzure, kBlue, kMagenta + 2,
                      kYellow + 2, kOrange + 2, kAzure + 5, kGreen + 2,
                     };
  int corrMarker[] = {kFullDiamond, kFullCross, kFullCircle, kFullSquare,
                      kOpenDiamond, kOpenCross,
                      kOpenSquare, kOpenCircle
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
  TGraphErrors *gc2c1[NCORR];


  //-- v2
  double vn_cent[NC][2];
  double vn_e_cent[NC][2];
  TGraphErrors *gvn_cent[2];

  double vn_pT[NC][2][NPT];
  double vn_e_pT[NC][2][NPT];
  TGraphErrors *gvn_pT[NC][2];

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
          dphi_FG12_pt[icorr][iz][ic][ipt]->Sumw2();
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
          dphi_BG12_pt[icorr][iz][ic][ipt]->Sumw2();
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
          if (iz == 0 && ipt == ptsuml)
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
          if (ipt > ptsuml && ipt <= ptsumh && iz > 0)
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
        dphi_FG12[icorr][iz][ic]->Sumw2();
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
        dphi_BG12[icorr][iz][ic]->Sumw2();
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
          // cout << " icorr:" << icorr << " ic:" << ic << " ipt:" << ipt << endl;

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
      // cout << " icorr:" << icorr << " ic:" << ic << endl;
      dphi_corr[icorr][ic] =
        (TH1D*) dphi_FGsum[icorr][ic]->Clone(
          Form("dphi_corr_%i_%i", icorr, ic));
      dphi_corr[icorr][ic]->SetDirectory(0);
      dphi_corr[icorr][ic]->Divide(dphi_BGsum[icorr][ic]);
      dphi_corr[icorr][ic]->Scale(
        dphi_BGsum[icorr][ic]->Integral() / dphi_FGsum[icorr][ic]->Integral());
      dphi_corr[icorr][ic]->SetMarkerStyle(kOpenCircle);
      dphi_corr[icorr][ic]->SetMarkerColor(kBlack);
    } // ic

  } // icorr


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
  // MAKE GRAPHS OF CORRELATION COEFFICIENTS
  //==========================================================================//
  cout << endl;
  cout << "--> Making Graphs of correlation coefficients" << endl;

  //-- pT dependent correlations (vs pT)
  for (int icorr = 0; icorr < NCORRPT; icorr++)
  {
    for (int ic = 0; ic < NC; ic++)
    {
      for (int ipar = 0; ipar < NPAR; ipar++)
      {
        gcn_pt[icorr][ic][ipar] = new TGraphErrors();
        gcn_pt[icorr][ic][ipar]->SetName(Form("gcn_pt_%i_%i_%i",
                                              icorr, ic, ipar));
        gcn_pt[icorr][ic][ipar]->SetMarkerStyle(corrMarker[icorr]);
        gcn_pt[icorr][ic][ipar]->SetMarkerColor(corrColor[icorr]);
        gcn_pt[icorr][ic][ipar]->SetLineColor(corrColor[icorr]);

        for (int ipt = 0; ipt < NPT; ipt++)
        {
          gcn_pt[icorr][ic][ipar]->SetPoint(ipt,
                                            0.5 * (ptl[ipt] + pth[ipt]),
                                            cn_pt[icorr][ic][ipt][ipar]);
          gcn_pt[icorr][ic][ipar]->SetPointError(ipt,
                                                 0,
                                                 cn_pt_e[icorr][ic][ipt][ipar]);
        } // ipt
      } // ipar
    } // ic
  } // icorr

  //-- pT integrated correlations (vs centrality)
  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    // C2/C1
    gc2c1[icorr] = new TGraphErrors();
    gc2c1[icorr]->SetName(Form("gc2c1_%i", icorr));
    gc2c1[icorr]->SetMarkerStyle(corrMarker[icorr]);
    gc2c1[icorr]->SetMarkerColor(corrColor[icorr]);
    gc2c1[icorr]->SetLineColor(corrColor[icorr]);
    gc2c1[icorr]->SetTitle(";centrality index;C_{2} / C_{1}");
    gc2c1[icorr]->GetYaxis()->CenterTitle();
    gc2c1[icorr]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {
      double c2c1 = cn[icorr][ic][1] / cn[icorr][ic][0];
      double unc = 0;
      unc += TMath::Power(cn_e[icorr][ic][0] / cn[icorr][ic][0], 2);
      unc += TMath::Power(cn_e[icorr][ic][1] / cn[icorr][ic][1], 2);
      unc = c2c1 * TMath::Sqrt(unc);

      gc2c1[icorr]->SetPoint(ic, ic + 0.5, c2c1);
      gc2c1[icorr]->SetPointError(ic, 0, unc);

      if (c2c1 > gc2c1[icorr]->GetMaximum())
        gc2c1[icorr]->SetMaximum(c2c1);
      if (c2c1 < gc2c1[icorr]->GetMinimum())
        gc2c1[icorr]->SetMinimum(c2c1);
    }

    // All Cn
    for (int ipar = 0; ipar < NPAR; ipar++)
    {
      gcn[icorr][ipar] = new TGraphErrors();
      gcn[icorr][ipar]->SetName(Form("gcn_%i_%i", icorr, ipar));
      gcn[icorr][ipar]->SetMarkerStyle(corrMarker[icorr]);
      gcn[icorr][ipar]->SetMarkerColor(corrColor[icorr]);
      gcn[icorr][ipar]->SetLineColor(corrColor[icorr]);
      gcn[icorr][ipar]->SetTitle(Form(";centrality index; C_{%i}", ipar));
      gcn[icorr][ipar]->SetMinimum(100);

      for (int ic = 0; ic < NC; ic++)
      {
        gcn[icorr][ipar]->SetPoint(ic, ic + 0.5, cn[icorr][ic][ipar]);
        gcn[icorr][ipar]->SetPointError(ic, 0, cn_e[icorr][ic][ipar]);

        if (cn[icorr][ic][ipar] > gcn[icorr][ipar]->GetMaximum())
          gcn[icorr][ipar]->SetMaximum(cn[icorr][ic][ipar]);
        if (cn[icorr][ic][ipar] < gcn[icorr][ipar]->GetMinimum())
          gcn[icorr][ipar]->SetMinimum(cn[icorr][ic][ipar]);
      }
    } // ipar
  } // icorr


  //==========================================================================//
  // CALCULATE V2 USING 3-SUB EVENT METHOD
  //==========================================================================//
  cout << endl;
  cout << "--> Calcualte v2 using 3 sub-event method" << endl;

  //--CNT-FVTXN-FVTXS
  for (int j = 0; j < 2; j++)
  {
    gvn_cent[j] = new TGraphErrors();
    sprintf(hname, "gvn_cent_n%i_dAu%i", j, energy);
    gvn_cent[j]->SetName(hname);
    gvn_cent[j]->SetMarkerStyle(20);
    gvn_cent[j]->SetMarkerColor(kBlue);
    gvn_cent[j]->SetLineColor(kBlue);
    gvn_cent[j]->SetMinimum(100);
    for (int ic = 0; ic < NC; ic++)
    {
      double c2_AB = cn[CNTFVTXS][ic][j + 1];
      double c2_AC = cn[CNTFVTXN][ic][j + 1];
      double c2_BC = cn[FVTXNFVTXS][ic][j + 1];
      vn_cent[ic][j] = TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);

      //uncertainty
      double dc2_AB = 0.5 * (c2_AC / c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);
      double dc2_AC = 0.5 * (c2_AB / c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);
      double dc2_BC = -0.5 * (c2_AB * c2_AC) / (c2_BC * c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);

      double sig_AB = cn_e[CNTFVTXS][ic][j + 1];
      double sig_AC = cn_e[CNTFVTXN][ic][j + 1];
      double sig_BC = cn_e[FVTXNFVTXS][ic][j + 1];

      double sig_vn = 0;
      sig_vn += TMath::Power(dc2_AB * sig_AB, 2);
      sig_vn += TMath::Power(dc2_AC * sig_AC, 2);
      sig_vn += TMath::Power(dc2_BC * sig_BC, 2);
      vn_e_cent[ic][j] = TMath::Sqrt(sig_vn);

      //fill tgraph
      gvn_cent[j]->SetPoint(ic, ic + 0.5, vn_cent[ic][j]);
      gvn_cent[j]->SetPointError(ic, 0, vn_e_cent[ic][j]);

      if (vn_cent[ic][j] > gvn_cent[j]->GetMaximum())
        gvn_cent[j]->SetMaximum(vn_cent[ic][j]);
      if (vn_cent[ic][j] < gvn_cent[j]->GetMinimum())
        gvn_cent[j]->SetMinimum(vn_cent[ic][j]);
    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gvn_pT_c%i_n%i_dAu%i", ic, j, energy);
      gvn_pT[ic][j]->SetName(hname);
      gvn_pT[ic][j]->SetMarkerStyle(20);
      gvn_pT[ic][j]->SetMarkerColor(kBlue);
      gvn_pT[ic][j]->SetLineColor(kBlue);
      gvn_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {
        double c2_AB = cn_pt[CNTFVTXS][ic][ipt][j + 1];
        double c2_AC = cn_pt[CNTFVTXN][ic][ipt][j + 1];
        double c2_BC = cn[FVTXNFVTXS][ic][j + 1];
        vn_pT[ic][j][ipt] = TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);

        //uncertainty
        double dc2_AB = 0.5 * (c2_AC / c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);
        double dc2_AC = 0.5 * (c2_AB / c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);
        double dc2_BC = -0.5 * (c2_AB * c2_AC) / (c2_BC * c2_BC) / TMath::Sqrt( (c2_AB * c2_AC) / c2_BC);

        double sig_AB = cn_pt_e[CNTFVTXS][ic][ipt][j + 1];
        double sig_AC = cn_pt_e[CNTFVTXN][ic][ipt][j + 1];
        double sig_BC = cn_e[FVTXNFVTXS][ic][j + 1];

        double sig_vn = 0;
        sig_vn += TMath::Power(dc2_AB * sig_AB, 2);
        sig_vn += TMath::Power(dc2_AC * sig_AC, 2);
        sig_vn += TMath::Power(dc2_BC * sig_BC, 2);
        vn_e_pT[ic][j][ipt] = TMath::Sqrt(sig_vn);

        //fill tgraph
        gvn_pT[ic][j]->SetPoint(ipt, 0.5 * (ptl[ipt] + pth[ipt]), vn_pT[ic][j][ipt]);
        gvn_pT[ic][j]->SetPointError(ipt, 0, vn_e_pT[ic][j][ipt]);

        if (vn_pT[ic][j][ipt] > gvn_pT[ic][j]->GetMaximum())
          gvn_pT[ic][j]->SetMaximum(vn_pT[ic][j][ipt]);
        if (vn_pT[ic][j][ipt] < gvn_pT[ic][j]->GetMinimum())
          gvn_pT[ic][j]->SetMinimum(vn_pT[ic][j][ipt]);

      }

    }
  } // j


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


  vector<int> corrdrawidx = {CNTBBCS, CNTFVTXS, CNTFVTXN, FVTXNFVTXS};
  TLegend *legcn = new TLegend(0.2, 0.2, 0.5, 0.45);
  legcn->SetFillStyle(0);
  legcn->SetBorderSize(0);
  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
    legcn->AddEntry(gcn[corrdrawidx.at(j)][0], cCorr[corrdrawidx.at(j)], "P");


  TH1F* haxis_cncent = new TH1F("haxis_cncent", "", NC, 0, NC);
  haxis_cncent->GetYaxis()->CenterTitle();
  haxis_cncent->GetYaxis()->SetTitleSize(0.08);
  haxis_cncent->GetYaxis()->SetTitleOffset(0.75);
  haxis_cncent->GetXaxis()->SetLabelSize(0.08);

  for (int ibin = 1; ibin <= haxis_cncent->GetNbinsX(); ibin++)
  {
    // if (ibin < 2 || ibin > NC)
    // haxis_cncent->GetXaxis()->SetBinLabel(ibin, "");
    // else
    // {
    sprintf(ctitle, "%i - %i%%", cl[ibin - 1], ch[ibin - 1]);
    haxis_cncent->GetXaxis()->SetBinLabel(ibin, ctitle);
    // }
  }

  float min, max;

  TH1D* hrat;

  //temporary bools to make some plotting quicker
  bool plot_CORRPT = false;
  bool plot_FGBG = true;
  bool plot_CORR = true;



  //==========================================================================//
  // PLOT
  //==========================================================================//



  TCanvas *cFGBGpt[NCORRPT][NPT];
  TCanvas *ccorrpt[NCORRPT][NPT];
  TCanvas *cFGBG[NCORR];
  TCanvas *cBGrat[NCORR];
  TCanvas *ccorr[NCORR];

  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    //-- pT dependent

    if (icorr < NCORRPT && plot_CORRPT)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
      {
        //FGBG
        if (plot_FGBG)
        {
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
        } // plot_FGBG


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
    if (plot_CORR)
    {
      //FGBG
      if (plot_FGBG)
      {
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

          dphi_FGsum[icorr][ic]->Draw("Hist,E,");
          dphi_BGsum[icorr][ic]->Draw("Hist,E,same");

          sprintf(ctitle, "%s FGBG %i-%i%%",
                  cCorr[icorr],
                  cl[ic], ch[ic]);
          ltitle.DrawLatex(0.5, 0.95, ctitle);

          if (ic == 0)
            le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
        } // ic


        sprintf(hname, "cBGrat_%i", icorr);
        cBGrat[icorr] = new TCanvas(hname, hname, 1200, 1000);
        cBGrat[icorr]->Divide(3, 2, 0, 0);
        for (int ic = 0; ic < NC; ic++)
        {
          cBGrat[icorr]->GetPad(ic + 1)->SetTopMargin(0.1);
          cBGrat[icorr]->GetPad(ic + 1)->SetRightMargin(0.02);
          cBGrat[icorr]->GetPad(ic + 1)->SetBottomMargin(0.1);
          cBGrat[icorr]->GetPad(ic + 1)->SetLeftMargin(0.1);
          cBGrat[icorr]->GetPad(ic + 1)->SetTicks(1, 1);

          cBGrat[icorr]->cd(ic + 1);

          hrat = (TH1D*) dphi_BGsum[icorr][ic]->Clone("hrat");
          hrat->Divide(dphi_BGsum[icorr][0]);
          hrat->DrawCopy("Hist,E");

          sprintf(ctitle, "%s BG ratio %i-%i%% / %i-%i%%",
                  cCorr[icorr],
                  cl[ic], ch[ic], cl[0], ch[0]);
          ltitle.DrawLatex(0.5, 0.95, ctitle);

          l1.DrawLine(-1 * TMath::Pi() / 2., 1, 3 * TMath::Pi() / 2., 1.);

          if (ic == 0)
            le.DrawLatex(0.2, 0.8, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));

          lc.SetTextColor(fitColor[0]);
          lc.DrawLatex(0.45, 0.7, Form("C_{1} = % .4f",
                                       cn[icorr][ic][0]));

          // lc.SetTextColor(kRed);
          // lc.DrawLatex(0.45, 0.64, Form("RMS = % .4f",
          //                               hrat->GetMean(2)));

        } // ic

      } // plot_FGBG

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

        sprintf(ctitle, "%s %i-%i%%",
                cCorr[icorr],
                cl[ic], ch[ic]);
        ltitle.DrawLatex(0.5, 0.95, ctitle);

        l1.DrawLine(-1 * TMath::Pi() / 2., 1, 3 * TMath::Pi() / 2., 1.);

        if (ic == 0)
          le.DrawLatex(0.2, 0.82, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
        if (icorr < NCORRPT)
          le.DrawLatex(0.15, 0.76, Form("%.1f<p_{T}^{trig} [GeV/c]<%.1f",
                                        ptl[ptsuml], pth[ptsumh]));


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
    } // plot_CORR

  } // icorr


  TCanvas *cncent = new TCanvas("cncent", "cncent", 800, 900);
  cncent->Divide(1, 3, 0, 0);

  cncent->GetPad(1)->SetTopMargin(0.1);
  cncent->GetPad(1)->SetRightMargin(0.02);
  cncent->GetPad(1)->SetBottomMargin(0.0);
  cncent->GetPad(1)->SetLeftMargin(0.12);
  cncent->GetPad(1)->SetTicks(1, 1);

  cncent->GetPad(2)->SetTopMargin(0.0);
  cncent->GetPad(2)->SetRightMargin(0.02);
  cncent->GetPad(2)->SetBottomMargin(0.0);
  cncent->GetPad(2)->SetLeftMargin(0.12);
  cncent->GetPad(2)->SetTicks(1, 1);

  cncent->GetPad(3)->SetTopMargin(0.0);
  cncent->GetPad(3)->SetRightMargin(0.02);
  cncent->GetPad(3)->SetBottomMargin(0.12);
  cncent->GetPad(3)->SetLeftMargin(0.12);
  cncent->GetPad(3)->SetTicks(1, 1);


  //--C_1
  cncent->cd(1);
  // get the min/max
  min = 100;
  max = -100;
  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
  {
    int idx = corrdrawidx.at(j);
    float hi = gcn[idx][0]->GetMaximum();
    float lo = gcn[idx][0]->GetMinimum();

    cout << "  C_1, " << cCorr[idx] << " lo:" << lo << " hi:" << hi << endl;

    if (hi > max)
      max = hi;
    if (lo < min)
      min = lo;
  }
  haxis_cncent->SetMinimum(1.1 * min);
  haxis_cncent->SetMaximum(0.9 * max);
  haxis_cncent->SetTitle(";;C_{1}");
  haxis_cncent->DrawCopy();

  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
    gcn[corrdrawidx.at(j)][0]->Draw("P");

  legcn->Draw("same");
  ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));

  //--C_2
  cncent->cd(2);
  // get the min/max
  min = 100;
  max = -100;
  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
  {
    int idx = corrdrawidx.at(j);
    float hi = gcn[idx][1]->GetMaximum();
    float lo = gcn[idx][1]->GetMinimum();

    cout << "  C_2, " << cCorr[idx] << " lo:" << lo << " hi:" << hi << endl;

    if (hi > max)
      max = hi;
    if (lo < min)
      min = lo;
  }
  haxis_cncent->SetMinimum(0.9 * min);
  haxis_cncent->SetMaximum(1.1 * max);
  haxis_cncent->SetTitle(";;C_{2}");
  haxis_cncent->DrawCopy();

  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
    gcn[corrdrawidx.at(j)][1]->Draw("P");

  //--C_2 / C_1
  cncent->cd(3);
  // get the min/max
  min = 100;
  max = -100;
  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
  {
    int idx = corrdrawidx.at(j);
    float hi = gc2c1[idx]->GetMaximum();
    float lo = gc2c1[idx]->GetMinimum();

    cout << "  C_2/C_1, " << cCorr[idx] << " lo:" << lo << " hi:" << hi << endl;

    if (hi > max)
      max = hi;
    if (lo < min)
      min = lo;
  }
  haxis_cncent->SetMinimum(1.1 * min);
  haxis_cncent->SetMaximum(0.9 * max);
  haxis_cncent->SetTitle(";;C_{2} / C_{1}");
  haxis_cncent->DrawCopy();

  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
    gc2c1[corrdrawidx.at(j)]->Draw("P");





  TCanvas *cv2cent = new TCanvas("cv2cent", "v2 cent", 800, 600);
  cv2cent->Divide(1, 2);

  for (int j = 0; j < 2; j++)
  {
    cv2cent->cd(j + 1);
    haxis_cncent->SetMinimum(0.0);
    haxis_cncent->SetMaximum(0.2);
    haxis_cncent->SetTitle(Form(";;v_{%i}", j + 2));
    haxis_cncent->GetYaxis()->SetTitleSize(0.06);
    haxis_cncent->GetXaxis()->SetLabelSize(0.06);
    haxis_cncent->DrawCopy();

    gvn_cent[j]->Draw("P");

    if (j == 0)
    {
      ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
      le.DrawLatex(0.15, 0.86, "3 sub-event #color[4]{CNT}-FVTXN-FVTXS");
      le.DrawLatex(0.15, 0.80, Form("%.1f<p_{T}^{trig} [GeV/c]<%.1f",
                                    ptl[ptsuml], pth[ptsumh]));
    }
  }


  TCanvas *cv2pt = new TCanvas("cv2pt", "v2 pt", 800, 600);
  cv2pt->Divide(1, 2);

  for (int j = 0; j < 2; j++)
  {

    cv2pt->cd(j + 1);
    gvn_pT[0][j]->SetMinimum(0);
    gvn_pT[0][j]->SetMaximum(0.2);
    gvn_pT[0][j]->GetXaxis()->SetRangeUser(0, 5);
    gvn_pT[0][j]->SetTitle(Form(";p_{T} [GeV/c];v_{%i}", j+2));
    gvn_pT[0][j]->Draw("AP");

    if (j == 0)
    {
      ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));
      le.DrawLatex(0.2, 0.8, "3 sub-event #color[4]{CNT}-FVTXN-FVTXS");
    }
  }


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
      if (icorr < NCORRPT && plot_CORRPT)
      {
        for (int ipt = 0; ipt < NPT; ipt++)
        {
          if (plot_FGBG)
          {
            sprintf(cname, "pdfs/dAu%i_%s_FGBG_pt%i.pdf", energy, cCorr[icorr], ipt);
            cFGBGpt[icorr][ipt]->Print(cname);
          }

          sprintf(cname, "pdfs/dAu%i_%s_corr_pt%i.pdf", energy, cCorr[icorr], ipt);
          ccorrpt[icorr][ipt]->Print(cname);
        }
      }

      //-- pT integrated
      if (plot_CORR)
      {
        if (plot_FGBG)
        {
          sprintf(cname, "pdfs/dAu%i_%s_FGBG.pdf", energy, cCorr[icorr]);
          cFGBG[icorr]->Print(cname);

          sprintf(cname, "pdfs/dAu%i_%s_BGrat.pdf", energy, cCorr[icorr]);
          cBGrat[icorr]->Print(cname);
        }

        sprintf(cname, "pdfs/dAu%i_%s_corr.pdf", energy, cCorr[icorr]);
        ccorr[icorr]->Print(cname);
      }
    } // icorr

    sprintf(cname, "pdfs/dAu%i_cncent.pdf", energy);
    cncent->Print(cname);

    sprintf(cname, "pdfs/dAu%i_v2cent.pdf", energy);
    cv2cent->Print(cname);

    sprintf(cname, "pdfs/dAu%i_v2pT.pdf", energy);
    cv2pt->Print(cname);
  } // printPlots

}