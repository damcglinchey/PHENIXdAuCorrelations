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
#include <utility>

#include "v2.h"

using namespace std;


typedef pair<double, double> ValErr;

ValErr vn_3sub(const ValErr &C_AB, const ValErr &C_AC, const ValErr &C_BC);
ValErr calc_cn(TH1D* hcorr, int order);


///-----------------------------------------------------------------------------
///
/// Main Function
///
/// The main switch below is "energy"
/// The names of the ROOT files containing the correlations and multiplicity
/// are generated automatically based on the energy:
///    correlation file: "rootfiles/CorrFuncdAuBES_dAu%i.root"
///    histogram file  : "rootfiles/Histos_CorrFuncdAuBES_dAu%i.root"
void plot_corrfuncs()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  //==========================================================================//
  // SET RUNNING CONDITIONS
  //==========================================================================//

  int energy = 20;

  bool printPlots = false;
  bool saveHistos = false;
  const char* outFile = "correlations.root";
  bool printTables = false;

  // plotting options
  bool plot_CORRPT = false;
  bool plot_FGBG = false;
  bool plot_CORR = true;


  const int NC  =  6; // number of centrality bins
  const int NZ  = 10; // number of z bins
  const int NPT =  7; // number of pT bins

  float ptl[] = {0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00};
  float pth[] = {0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 5.00};
  int cl[] = {0,  5, 10, 20, 40,  60};
  int ch[] = {5, 10, 20, 40, 60, 100};
  int centColor[] = {kBlack, kBlue, kRed, kGreen + 2, kMagenta + 2, kOrange + 2};

  int ptsuml = 3;
  int ptsumh = 5;
  // int ptsuml = 0;
  // int ptsumh = NPT - 1;



  // correlations between detectors
  const int NCORRPT = 4; // number of correlations with pT binning
  const int NCORR = 9;   // number of correlations
  enum CORR {
    CNTBBCS = 0,
    CNTBBCN,
    CNTFVTXS,
    CNTFVTXN,
    BBCNBBCS,
    BBCNFVTXS,
    FVTXNBBCS,
    FVTXNFVTXS,
    BBCSFVTXS,
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
    "BBCSFVTXS",
  };
  int corrColor[] =  {kRed, kAzure, kBlue, kMagenta + 2,
                      kYellow + 2, kOrange + 2, kAzure + 5, kGreen + 2,
                      kGray,
                     };
  int corrMarker[] = {kFullDiamond, kFullCross, kFullCircle, kFullSquare,
                      kOpenDiamond, kOpenCross,
                      kOpenSquare, kOpenCircle,
                      kFullStar,
                     };

  // Rebin values
  // Initialize all to 4, but make some 8
  vector<int> corrRebin(NCORR, 4);
  corrRebin[BBCNBBCS] = 8;
  corrRebin[CNTBBCN] = 8;
  if (energy == 20)
  {
    corrRebin[CNTBBCS] = 8;
    corrRebin[CNTFVTXS] = 8;
    corrRebin[CNTFVTXN] = 8;
    corrRebin[BBCSFVTXS] = 8;
  }

  const int NPAR = 4; // number of orders in the fourier fit
  int fitColor[] = {kBlue, kRed, kGreen + 2, kMagenta + 2, kYellow + 2};

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

  TH1D* dphi_BGsum_tot[NCORR];
  TH1D* dphi_corr_tot[NCORR][NC];

  // multiplicity
  TH2D* hnpc1centrality;
  TH1D* hnpc1_cent[NC];

  // deta
  TH1D* deta[NCORR][NC];
  double mdeta[NCORR][NC] = {0};


  //-- fitting

  // TF1* fcorr = new TF1("fcorr",
  //                      "1 + 2*[0]*TMath::Cos(x) + 2*[1]*TMath::Cos(2*x) + 2*[2]*TMath::Cos(3*x) + 2*[3]*TMath::Cos(4*x) + 2*[4]*TMath::Cos(5*x)",
  //                      -1 * TMath::Pi() / 2., 3 * TMath::Pi() / 2.);
  TF1* fcorr = new TF1("fcorr",
                       "1 + 2*[0]*TMath::Cos(x) + 2*[1]*TMath::Cos(2*x) + 2*[2]*TMath::Cos(3*x) + 2*[3]*TMath::Cos(4*x)",
                       -1 * TMath::Pi() / 2., 3 * TMath::Pi() / 2.);
  fcorr->SetLineColor(kBlack);
  fcorr->SetLineStyle(2);


  TF1* fpol = new TF1("fpol", "pol4", 0.3, 3);
  fpol->SetLineColor(kBlack);
  fpol->SetLineStyle(2);


  //-- cn's
  ValErr cn_pt[NCORRPT][NC][NPT][NPAR];
  ValErr cn[NCORR][NC][NPAR];

  // calculated numerically (i.e. not from fit)
  ValErr cn_num_pt[NCORRPT][NC][NPT][NPAR];
  ValErr cn_num[NCORR][NC][NPAR];

  //-- Multplicity
  ValErr mean_npc1[NC];

  //-- cn graphs
  TGraphErrors *gcn_pt[NCORRPT][NC][NPAR];
  TGraphErrors *gcn[NCORR][NPAR];
  TGraphErrors *gc2c1[NCORR];

  TGraphErrors *gc2c1_deta[NC];


  //-- vn (3 sub-event method)
  ValErr vn_CNTFVTXSFVTXN_cent[NC][NPAR];
  ValErr vn_CNTFVTXSFVTXN_pT[NC][NPAR][NPT];
  TGraphErrors *gvn_CNTFVTXSFVTXN_cent[NPAR];
  TGraphErrors *gvn_CNTFVTXSFVTXN_mult[NPAR];
  TGraphErrors *gvn_CNTFVTXSFVTXN_pT[NC][NPAR];

  ValErr vn_CNTBBCSFVTXS_cent[NC][NPAR];
  ValErr vn_CNTBBCSFVTXS_pT[NC][NPAR][NPT];
  TGraphErrors *gvn_CNTBBCSFVTXS_cent[NPAR];
  TGraphErrors *gvn_CNTBBCSFVTXS_mult[NPAR];
  TGraphErrors *gvn_CNTBBCSFVTXS_pT[NC][NPAR];

  ValErr vn_CNTBBCSFVTXN_cent[NC][NPAR];
  ValErr vn_CNTBBCSFVTXN_pT[NC][NPAR][NPT];
  TGraphErrors *gvn_CNTBBCSFVTXN_cent[NPAR];
  TGraphErrors *gvn_CNTBBCSFVTXN_mult[NPAR];
  TGraphErrors *gvn_CNTBBCSFVTXN_pT[NC][NPAR];

  ValErr vn_CNTFVTXSBBCN_cent[NC][NPAR];
  ValErr vn_CNTFVTXSBBCN_pT[NC][NPAR][NPT];
  TGraphErrors *gvn_CNTFVTXSBBCN_cent[NPAR];
  TGraphErrors *gvn_CNTFVTXSBBCN_mult[NPAR];
  TGraphErrors *gvn_CNTFVTXSBBCN_pT[NC][NPAR];

  ValErr vn_CNTBBCSBBCN_cent[NC][NPAR];
  ValErr vn_CNTBBCSBBCN_pT[NC][NPAR][NPT];
  TGraphErrors *gvn_CNTBBCSBBCN_cent[NPAR];
  TGraphErrors *gvn_CNTBBCSBBCN_mult[NPAR];
  TGraphErrors *gvn_CNTBBCSBBCN_pT[NC][NPAR];

  ValErr vn_num_CNTBBCSFVTXN_cent[NC][NPAR];
  TGraphErrors *gvn_num_CNTBBCSFVTXN_mult[NPAR];


  TGraphErrors* gv2_Ron;
  TGraphErrors* gv2fitrat_Ron;
  TGraphErrors* gv2fitrat_CNTFVTXSFVTXN;
  TGraphErrors* gv2fitrat_CNTBBCSFVTXS;
  TGraphErrors* gv2fitrat_CNTBBCSFVTXN;
  TGraphErrors* gv2fitrat_CNTFVTXSBBCN;
  TGraphErrors* gv2fitrat_CNTBBCSBBCN;


  //-- temp stuff
  char hname[500];
  char fname[500];

  TH2D* h2tmp;
  TH1D* h1tmp;

  //==========================================================================//
  // GET CORRELATION HISTOGRAMS
  //==========================================================================//
  cout << endl;

  sprintf(fname, "rootfiles/CorrFuncdAuBES_dAu%i.root", energy);
  TFile *fin = TFile::Open(fname);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << fname << endl;
    return;
  }
  cout << "--> For dAu " << energy << " GeV, "
       << "reading correlation histograms from " << fname << endl;

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
          dphi_FG12_pt[icorr][iz][ic][ipt]->Rebin(corrRebin[icorr]);


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
          dphi_BG12_pt[icorr][iz][ic][ipt]->Rebin(corrRebin[icorr]);


          // deta
          sprintf(hname,
                  "deta_dphi_%s_%03i_%03i_%03.0f_%03.0f_FG12",
                  cCorr[icorr],
                  cl[ic], ch[ic],
                  ptl[ipt] * 100., pth[ipt] * 100);

          h2tmp = (TH2D*) gROOT->FindObject(hname);
          if (!h2tmp)
          {
            cout << "ERROR!! Unable to find " << hname << " for "
                 << " ic:" << ic << " iz:" << iz << " ipt:" << ipt
                 << endl;
            return;
          }
          h1tmp = (TH1D*) h2tmp->ProjectionX("h1tmp");



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

            dphi_BGsum_tot[icorr] =
              (TH1D*) dphi_BG12_pt[icorr][iz][ic][ipt]->Clone(
                Form("dphi_BGsum_tot_%s", cCorr[icorr]));
            dphi_BGsum_tot[icorr]->SetDirectory(0);

            deta[icorr][ic] =
              (TH1D*) h1tmp->Clone(
                Form("deta_%s_%i", cCorr[icorr], ic));
            deta[icorr][ic]->SetDirectory(0);

          }
          if (ipt > ptsuml && ipt <= ptsumh && iz > 0)
          {
            dphi_FGsum[icorr][ic]->Add(dphi_FG12_pt[icorr][iz][ic][ipt]);
            dphi_BGsum[icorr][ic]->Add(dphi_BG12_pt[icorr][iz][ic][ipt]);

            dphi_BGsum_tot[icorr]->Add(dphi_BG12_pt[icorr][iz][ic][ipt]);

            deta[icorr][ic]->Add(h1tmp);

          }

          delete h1tmp;
          delete h2tmp;

        } // ipt
      } // icorr



      //-- Get the pT integrated correlations
      for (int icorr = NCORRPT; icorr < NCORR; icorr++)
      {

        // FG12
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
        dphi_FG12[icorr][iz][ic]->Rebin(corrRebin[icorr]);


        // BG12
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
        dphi_BG12[icorr][iz][ic]->Rebin(corrRebin[icorr]);

        // deta
        sprintf(hname,
                "deta_dphi_%s_%03i_%03i_FG12",
                cCorr[icorr],
                cl[ic], ch[ic]);

        h2tmp = (TH2D*) gROOT->FindObject(hname);
        if (!h2tmp)
        {
          cout << "ERROR!! Unable to find " << hname << " for "
               << " ic:" << ic << " iz:" << iz
               << endl;
          return;
        }
        h1tmp = (TH1D*) h2tmp->ProjectionX("h1tmp");

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

          dphi_BGsum_tot[icorr] =
            (TH1D*) dphi_BG12[icorr][iz][ic]->Clone(
              Form("dphi_BGsum_tot_%s", cCorr[icorr]));
          dphi_BGsum_tot[icorr]->SetDirectory(0);

          deta[icorr][ic] =
            (TH1D*) h1tmp->Clone(
              Form("deta_%s_%i", cCorr[icorr], ic));
          deta[icorr][ic]->SetDirectory(0);

        }
        else
        {
          dphi_FGsum[icorr][ic]->Add(dphi_FG12[icorr][iz][ic]);
          dphi_BGsum[icorr][ic]->Add(dphi_BG12[icorr][iz][ic]);

          dphi_BGsum_tot[icorr]->Add(dphi_BG12[icorr][iz][ic]);

          deta[icorr][ic]->Add(h1tmp);
        }

        delete h1tmp;
        delete h2tmp;

      } // icorr

    } // ic

    //-- change back to the main directory
    fin->cd();

  } // iz

  fin->Close();
  delete fin;

  //==========================================================================//
  // READ THE MULTIPLICITY HISTOGRAMS FROM FILE
  //==========================================================================//
  cout << endl;

  sprintf(fname, "rootfiles/Histos_CorrFuncdAuBES_dAu%i.root", energy);
  fin = TFile::Open(fname);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << fname << endl;
    return;
  }
  cout << "--> For dAu" << energy << " Getting histograms from:"
       << fname << endl;


  //-- Get PC1 multiplicity in each centrality bin
  hnpc1centrality = (TH2D*) fin->Get("hnpc1_centrality");
  if (!hnpc1centrality)
  {
    cout << "ERROR!! Unable to find hnpc1_centrality in " << fname << endl;
    return;
  }
  hnpc1centrality->SetName(Form("hnpc1centrality_dAu%i", energy));
  hnpc1centrality->SetDirectory(0);

  for (int ic = 0; ic < NC; ic++)
  {
    int bl = hnpc1centrality->GetXaxis()->FindBin(cl[ic]);
    int bh = hnpc1centrality->GetXaxis()->FindBin(ch[ic] - 1);

    hnpc1_cent[ic] = (TH1D*) hnpc1centrality->ProjectionY(
                       Form("hnpc1_cent_dAu%i_c%i", energy, ic),
                       bl, bh);
    hnpc1_cent[ic]->SetDirectory(0);
    hnpc1_cent[ic]->SetLineColor(centColor[ic]);
    hnpc1_cent[ic]->Scale(1. / ((double)(ch[ic] - cl[ic])));

    mean_npc1[ic].first = hnpc1_cent[ic]->GetMean();
    mean_npc1[ic].second = hnpc1_cent[ic]->GetRMS();

  }



  //==========================================================================//
  // CALCULTATE MEAN DELTA ETA
  //==========================================================================//
  cout << endl;
  cout << "--> Calculating the mean delta eta" << endl;

  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    cout << "  " << cCorr[icorr] << endl;
    for (int ic = 0; ic < NC; ic++)
    {
      mdeta[icorr][ic] = deta[icorr][ic]->GetMean();

      //-- Temporary until 2D histos get fixed!!
      if (mdeta[icorr][ic] == 0 && icorr == BBCNBBCS)
        mdeta[icorr][ic] = 7.0;
      if (mdeta[icorr][ic] == 0 && icorr == FVTXNFVTXS)
        mdeta[icorr][ic] = 4.0;
      if (mdeta[icorr][ic] == 0 && icorr == BBCNFVTXS)
        mdeta[icorr][ic] = 5.5;
      if (mdeta[icorr][ic] == 0 && icorr == FVTXNBBCS)
        mdeta[icorr][ic] = 5.5;



      cout << "    " << cl[ic] << "-" << ch[ic] << ": " << mdeta[icorr][ic] << endl;
    }
  }



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
          dphi_corr_pt[icorr][ic][ipt]->SetTitle(";#Delta#phi;C(#Delta#Phi)");
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
      // dphi_corr[icorr][ic]->SetTitle(";#Delta#phi;C(#Delta#Phi)");
      dphi_corr[icorr][ic]->SetTitle(";#Delta#phi;");
      dphi_corr[icorr][ic]->Divide(dphi_BGsum[icorr][ic]);
      dphi_corr[icorr][ic]->Scale(
        dphi_BGsum[icorr][ic]->Integral() / dphi_FGsum[icorr][ic]->Integral());
      dphi_corr[icorr][ic]->SetMarkerStyle(kOpenCircle);
      dphi_corr[icorr][ic]->SetMarkerColor(kBlack);
    } // ic

  } // icorr


  //==========================================================================//
  // CALCULATE THE CORRELATION FUNCTIONS WITH SUMMED BG
  //==========================================================================//
  cout << endl;
  cout << "--> Calculating the correlation functions" << endl;

  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    for (int ic = 0; ic < NC; ic++)
    {

      dphi_corr_tot[icorr][ic] =
        (TH1D*) dphi_FGsum[icorr][ic]->Clone(
          Form("dphi_corr_tot_%i", icorr));
      dphi_corr_tot[icorr][ic]->SetDirectory(0);
      dphi_corr_tot[icorr][ic]->Divide(dphi_BGsum_tot[icorr]);
      dphi_corr_tot[icorr][ic]->Scale(
        dphi_BGsum_tot[icorr]->Integral() / dphi_FGsum[icorr][ic]->Integral());
      dphi_corr_tot[icorr][ic]->SetMarkerStyle(kOpenCircle);
      dphi_corr_tot[icorr][ic]->SetMarkerColor(kBlack);
    } // ic

  } // icorr


  //==========================================================================//
  // CALCULTE C_N's NUMERICALLY FROM CORRELATION FUNCTION
  //==========================================================================//
  cout << endl;
  cout << "--> Calculating C_n's numerically" << endl;

  for (int icorr = 0; icorr < NCORR; icorr++)
  {
    for (int ic = 0; ic < NC; ic++)
    {
      //-- pT dependent
      if (icorr < NCORRPT)
      {
        for (int ipt = 0; ipt < NPT; ipt++)
        {
          for (int i = 0; i < NPAR; i++)
            cn_num_pt[icorr][ic][ipt][i] = calc_cn(dphi_corr_pt[icorr][ic][ipt], i + 1);
        } // ipt
      }

      //-- pT integrated
      for (int i = 0; i < NPAR; i++)
        cn_num[icorr][ic][i] = calc_cn(dphi_corr[icorr][ic], i + 1);

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
          dphi_corr_pt[icorr][ic][ipt]->Fit(fcorr, "RQ0");
          for (int i = 0; i < NPAR; i++)
          {
            cn_pt[icorr][ic][ipt][i].first = fcorr->GetParameter(i);
            cn_pt[icorr][ic][ipt][i].second = fcorr->GetParError(i);
          }
        } // ipt
      }

      //-- pT integrated

      dphi_corr[icorr][ic]->Fit(fcorr, "RQ0");
      for (int i = 0; i < NPAR; i++)
      {
        cn[icorr][ic][i].first = fcorr->GetParameter(i);
        cn[icorr][ic][i].second = fcorr->GetParError(i);

        if (icorr == 0 && i < 3)
        {
          ValErr test = calc_cn(dphi_corr[icorr][ic], i + 1);
          cout << " "
               << cCorr[icorr]
               << " n=" << i + 1
               << " Cn(fit):" << cn[icorr][ic][i].first << " +/- " << cn[icorr][ic][i].second
               << " Cn(num):" << test.first << " +/- " << test.second
               << endl;
        }

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
                                            cn_pt[icorr][ic][ipt][ipar].first);
          gcn_pt[icorr][ic][ipar]->SetPointError(ipt,
                                                 0,
                                                 cn_pt[icorr][ic][ipt][ipar].second);
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
      double c2c1 = -1 * cn[icorr][ic][1].first / cn[icorr][ic][0].first;
      double unc = 0;
      unc += TMath::Power(cn[icorr][ic][0].second / cn[icorr][ic][0].first, 2);
      unc += TMath::Power(cn[icorr][ic][1].second / cn[icorr][ic][1].first, 2);
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
        gcn[icorr][ipar]->SetPoint(ic, ic + 0.5, cn[icorr][ic][ipar].first);
        gcn[icorr][ipar]->SetPointError(ic, 0, cn[icorr][ic][ipar].second);

        if (cn[icorr][ic][ipar].first > gcn[icorr][ipar]->GetMaximum())
          gcn[icorr][ipar]->SetMaximum(cn[icorr][ic][ipar].first);
        if (cn[icorr][ic][ipar].first < gcn[icorr][ipar]->GetMinimum())
          gcn[icorr][ipar]->SetMinimum(cn[icorr][ic][ipar].first);
      }
    } // ipar
  } // icorr


  //-- C2/C1 vs deta
  for (int ic = 0; ic < NC; ic++)
  {
    gc2c1_deta[ic] = new TGraphErrors();
    gc2c1_deta[ic]->SetName(Form("gc2c1_deta_%i", ic));
    gc2c1_deta[ic]->SetMarkerStyle(kFullCircle);
    gc2c1_deta[ic]->SetMarkerColor(kBlue);
    gc2c1_deta[ic]->SetLineColor(kBlue);
    gc2c1_deta[ic]->SetTitle(";#Delta#eta;-C_{2}/C_{1}");

    int ip = 0;
    for (int icorr = 0; icorr < NCORRPT; icorr++)
    {
      // if (icorr == BBCNFVTXS ||
      //     icorr == CNTBBCN)
      //   continue;

      double c2c1 = -1 * cn[icorr][ic][1].first / cn[icorr][ic][0].first;
      double unc = 0;
      unc += TMath::Power(cn[icorr][ic][0].second / cn[icorr][ic][0].first, 2);
      unc += TMath::Power(cn[icorr][ic][1].second / cn[icorr][ic][1].first, 2);
      unc = c2c1 * TMath::Sqrt(unc);

      gc2c1_deta[ic]->SetPoint(ip, -1 * mdeta[icorr][ic], c2c1);
      gc2c1_deta[ic]->SetPointError(ip, 0, unc);
      ip++;

      if (c2c1 < gc2c1_deta[ic]->GetMinimum())
        gc2c1_deta[ic]->SetMinimum(c2c1);
      if (c2c1 > gc2c1_deta[ic]->GetMaximum())
        gc2c1_deta[ic]->SetMaximum(c2c1);
    }
  }


  //==========================================================================//
  // CALCULATE V2 USING 3-SUB EVENT METHOD
  //==========================================================================//
  cout << endl;
  cout << "--> Calcualte v2 using 3 sub-event method" << endl;

  double ptshift = 0.02;
  double cshift  = 0.05;

  //--CNT-FVTXN-FVTXS
  for (int j = 0; j < NPAR; j++)
  {
    // centrality dependence
    gvn_CNTFVTXSFVTXN_cent[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTFVTXSFVTXN_cent", j + 1, energy);
    gvn_CNTFVTXSFVTXN_cent[j]->SetName(hname);
    gvn_CNTFVTXSFVTXN_cent[j]->SetMarkerStyle(20);
    gvn_CNTFVTXSFVTXN_cent[j]->SetMarkerColor(kBlue);
    gvn_CNTFVTXSFVTXN_cent[j]->SetLineColor(kBlue);
    gvn_CNTFVTXSFVTXN_cent[j]->SetMinimum(100);

    // multiplicity dependence
    gvn_CNTFVTXSFVTXN_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTFVTXSFVTXN_mult", j + 1, energy);
    gvn_CNTFVTXSFVTXN_mult[j]->SetName(hname);
    gvn_CNTFVTXSFVTXN_mult[j]->SetMarkerStyle(20);
    gvn_CNTFVTXSFVTXN_mult[j]->SetMarkerColor(kBlue);
    gvn_CNTFVTXSFVTXN_mult[j]->SetLineColor(kBlue);
    gvn_CNTFVTXSFVTXN_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_CNTFVTXSFVTXN_cent[ic][j] = vn_3sub(cn[CNTFVTXS][ic][j],
                                             cn[CNTFVTXN][ic][j],
                                             cn[FVTXNFVTXS][ic][j] );

      //fill centrality tgraph
      gvn_CNTFVTXSFVTXN_cent[j]->SetPoint(ic, ic + 0.5, vn_CNTFVTXSFVTXN_cent[ic][j].first);
      gvn_CNTFVTXSFVTXN_cent[j]->SetPointError(ic, 0, vn_CNTFVTXSFVTXN_cent[ic][j].second);

      if (vn_CNTFVTXSFVTXN_cent[ic][j].first > gvn_CNTFVTXSFVTXN_cent[j]->GetMaximum())
        gvn_CNTFVTXSFVTXN_cent[j]->SetMaximum(vn_CNTFVTXSFVTXN_cent[ic][j].first);
      if (vn_CNTFVTXSFVTXN_cent[ic][j].first < gvn_CNTFVTXSFVTXN_cent[j]->GetMinimum())
        gvn_CNTFVTXSFVTXN_cent[j]->SetMinimum(vn_CNTFVTXSFVTXN_cent[ic][j].first);

      // fill multiplicity tgraph
      gvn_CNTFVTXSFVTXN_mult[j]->SetPoint(ic,
                                          mean_npc1[ic].first,
                                          vn_CNTFVTXSFVTXN_cent[ic][j].first);
      gvn_CNTFVTXSFVTXN_mult[j]->SetPointError(ic, 0, vn_CNTFVTXSFVTXN_cent[ic][j].second);

      if (vn_CNTFVTXSFVTXN_cent[ic][j].first > gvn_CNTFVTXSFVTXN_mult[j]->GetMaximum())
        gvn_CNTFVTXSFVTXN_mult[j]->SetMaximum(vn_CNTFVTXSFVTXN_cent[ic][j].first);
      if (vn_CNTFVTXSFVTXN_cent[ic][j].first < gvn_CNTFVTXSFVTXN_mult[j]->GetMinimum())
        gvn_CNTFVTXSFVTXN_mult[j]->SetMinimum(vn_CNTFVTXSFVTXN_cent[ic][j].first);
    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_CNTFVTXSFVTXN_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gv%i_dAu%i_CNTFVTXSFVTXN_pT_c%i", j + 1, energy, ic);
      gvn_CNTFVTXSFVTXN_pT[ic][j]->SetName(hname);
      gvn_CNTFVTXSFVTXN_pT[ic][j]->SetMarkerStyle(20);
      gvn_CNTFVTXSFVTXN_pT[ic][j]->SetMarkerColor(kBlue);
      gvn_CNTFVTXSFVTXN_pT[ic][j]->SetLineColor(kBlue);
      gvn_CNTFVTXSFVTXN_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {

        vn_CNTFVTXSFVTXN_pT[ic][j][ipt] = vn_3sub(cn_pt[CNTFVTXS][ic][ipt][j],
                                          cn_pt[CNTFVTXN][ic][ipt][j],
                                          cn[FVTXNFVTXS][ic][j]);

        //fill tgraph
        gvn_CNTFVTXSFVTXN_pT[ic][j]->SetPoint(ipt,
                                              0.5 * (ptl[ipt] + pth[ipt]),
                                              vn_CNTFVTXSFVTXN_pT[ic][j][ipt].first);
        gvn_CNTFVTXSFVTXN_pT[ic][j]->SetPointError(ipt, 0, vn_CNTFVTXSFVTXN_pT[ic][j][ipt].second);

        if (vn_CNTFVTXSFVTXN_pT[ic][j][ipt].first > gvn_CNTFVTXSFVTXN_pT[ic][j]->GetMaximum())
          gvn_CNTFVTXSFVTXN_pT[ic][j]->SetMaximum(vn_CNTFVTXSFVTXN_pT[ic][j][ipt].first);
        if (vn_CNTFVTXSFVTXN_pT[ic][j][ipt].first < gvn_CNTFVTXSFVTXN_pT[ic][j]->GetMinimum())
          gvn_CNTFVTXSFVTXN_pT[ic][j]->SetMinimum(vn_CNTFVTXSFVTXN_pT[ic][j][ipt].first);

      }

    }
  } // j


  //--CNT-BBCS-FVTXS
  for (int j = 0; j < NPAR; j++)
  {
    // centrality dependence
    gvn_CNTBBCSFVTXS_cent[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXS_cent", j + 1, energy);
    gvn_CNTBBCSFVTXS_cent[j]->SetName(hname);
    gvn_CNTBBCSFVTXS_cent[j]->SetMarkerStyle(kFullSquare);
    gvn_CNTBBCSFVTXS_cent[j]->SetMarkerColor(kRed);
    gvn_CNTBBCSFVTXS_cent[j]->SetLineColor(kRed);
    gvn_CNTBBCSFVTXS_cent[j]->SetMinimum(100);

    // multiplicity dependence
    gvn_CNTBBCSFVTXS_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXS_mult", j + 1, energy);
    gvn_CNTBBCSFVTXS_mult[j]->SetName(hname);
    gvn_CNTBBCSFVTXS_mult[j]->SetMarkerStyle(kFullSquare);
    gvn_CNTBBCSFVTXS_mult[j]->SetMarkerColor(kRed);
    gvn_CNTBBCSFVTXS_mult[j]->SetLineColor(kRed);
    gvn_CNTBBCSFVTXS_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_CNTBBCSFVTXS_cent[ic][j] = vn_3sub(
                                      cn[CNTBBCS][ic][j],
                                      cn[CNTFVTXS][ic][j],
                                      cn[BBCSFVTXS][ic][j] );

      //fill centrality tgraph
      gvn_CNTBBCSFVTXS_cent[j]->SetPoint(ic, ic + 0.5 - cshift, vn_CNTBBCSFVTXS_cent[ic][j].first);
      gvn_CNTBBCSFVTXS_cent[j]->SetPointError(ic, 0, vn_CNTBBCSFVTXS_cent[ic][j].second);

      if (vn_CNTBBCSFVTXS_cent[ic][j].first > gvn_CNTBBCSFVTXS_cent[j]->GetMaximum())
        gvn_CNTBBCSFVTXS_cent[j]->SetMaximum(vn_CNTBBCSFVTXS_cent[ic][j].first);
      if (vn_CNTBBCSFVTXS_cent[ic][j].first < gvn_CNTBBCSFVTXS_cent[j]->GetMinimum())
        gvn_CNTBBCSFVTXS_cent[j]->SetMinimum(vn_CNTBBCSFVTXS_cent[ic][j].first);

      //fill multiplicity tgraph
      gvn_CNTBBCSFVTXS_mult[j]->SetPoint(ic,
                                         mean_npc1[ic].first,
                                         vn_CNTBBCSFVTXS_cent[ic][j].first);
      gvn_CNTBBCSFVTXS_mult[j]->SetPointError(ic, 0, vn_CNTBBCSFVTXS_cent[ic][j].second);

      if (vn_CNTBBCSFVTXS_cent[ic][j].first > gvn_CNTBBCSFVTXS_mult[j]->GetMaximum())
        gvn_CNTBBCSFVTXS_mult[j]->SetMaximum(vn_CNTBBCSFVTXS_cent[ic][j].first);
      if (vn_CNTBBCSFVTXS_cent[ic][j].first < gvn_CNTBBCSFVTXS_mult[j]->GetMinimum())
        gvn_CNTBBCSFVTXS_mult[j]->SetMinimum(vn_CNTBBCSFVTXS_cent[ic][j].first);
    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_CNTBBCSFVTXS_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXS_pT_c%i", j + 1, energy, ic);
      gvn_CNTBBCSFVTXS_pT[ic][j]->SetName(hname);
      gvn_CNTBBCSFVTXS_pT[ic][j]->SetMarkerStyle(kFullSquare);
      gvn_CNTBBCSFVTXS_pT[ic][j]->SetMarkerColor(kRed);
      gvn_CNTBBCSFVTXS_pT[ic][j]->SetLineColor(kRed);
      gvn_CNTBBCSFVTXS_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {

        vn_CNTBBCSFVTXS_pT[ic][j][ipt] = vn_3sub(
                                           cn_pt[CNTBBCS][ic][ipt][j],
                                           cn_pt[CNTFVTXS][ic][ipt][j],
                                           cn[BBCSFVTXS][ic][j] );

        //fill tgraph
        gvn_CNTBBCSFVTXS_pT[ic][j]->SetPoint(ipt,
                                             0.5 * (ptl[ipt] + pth[ipt]) - ptshift,
                                             vn_CNTBBCSFVTXS_pT[ic][j][ipt].first);
        gvn_CNTBBCSFVTXS_pT[ic][j]->SetPointError(ipt, 0, vn_CNTBBCSFVTXS_pT[ic][j][ipt].second);

        if (vn_CNTBBCSFVTXS_pT[ic][j][ipt].first > gvn_CNTBBCSFVTXS_pT[ic][j]->GetMaximum())
          gvn_CNTBBCSFVTXS_pT[ic][j]->SetMaximum(vn_CNTBBCSFVTXS_pT[ic][j][ipt].first);
        if (vn_CNTBBCSFVTXS_pT[ic][j][ipt].first < gvn_CNTBBCSFVTXS_pT[ic][j]->GetMinimum())
          gvn_CNTBBCSFVTXS_pT[ic][j]->SetMinimum(vn_CNTBBCSFVTXS_pT[ic][j][ipt].first);

      }

    }
  } // j


  //--CNT-BBCS-FVTXN
  for (int j = 0; j < NPAR; j++)
  {
    // centrality dependence
    gvn_CNTBBCSFVTXN_cent[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXN_cent", j + 1, energy);
    gvn_CNTBBCSFVTXN_cent[j]->SetName(hname);
    gvn_CNTBBCSFVTXN_cent[j]->SetMarkerStyle(kFullDiamond);
    gvn_CNTBBCSFVTXN_cent[j]->SetMarkerColor(kGreen + 2);
    gvn_CNTBBCSFVTXN_cent[j]->SetLineColor(kGreen + 2);
    gvn_CNTBBCSFVTXN_cent[j]->SetMinimum(100);

    // multiplicity dependence
    gvn_CNTBBCSFVTXN_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXN_mult", j + 1, energy);
    gvn_CNTBBCSFVTXN_mult[j]->SetName(hname);
    gvn_CNTBBCSFVTXN_mult[j]->SetMarkerStyle(kFullDiamond);
    gvn_CNTBBCSFVTXN_mult[j]->SetMarkerColor(kGreen + 2);
    gvn_CNTBBCSFVTXN_mult[j]->SetLineColor(kGreen + 2);
    gvn_CNTBBCSFVTXN_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_CNTBBCSFVTXN_cent[ic][j] = vn_3sub(
                                      cn[CNTBBCS][ic][j],
                                      cn[CNTFVTXN][ic][j],
                                      cn[FVTXNBBCS][ic][j] );

      //fill centrality tgraph
      gvn_CNTBBCSFVTXN_cent[j]->SetPoint(ic, ic + 0.5 + cshift, vn_CNTBBCSFVTXN_cent[ic][j].first);
      gvn_CNTBBCSFVTXN_cent[j]->SetPointError(ic, 0, vn_CNTBBCSFVTXN_cent[ic][j].second);

      if (vn_CNTBBCSFVTXN_cent[ic][j].first > gvn_CNTBBCSFVTXN_cent[j]->GetMaximum())
        gvn_CNTBBCSFVTXN_cent[j]->SetMaximum(vn_CNTBBCSFVTXN_cent[ic][j].first);
      if (vn_CNTBBCSFVTXN_cent[ic][j].first < gvn_CNTBBCSFVTXN_cent[j]->GetMinimum())
        gvn_CNTBBCSFVTXN_cent[j]->SetMinimum(vn_CNTBBCSFVTXN_cent[ic][j].first);

      //fill multiplicity tgraph
      gvn_CNTBBCSFVTXN_mult[j]->SetPoint(ic,
                                         mean_npc1[ic].first,
                                         vn_CNTBBCSFVTXN_cent[ic][j].first);
      gvn_CNTBBCSFVTXN_mult[j]->SetPointError(ic, 0, vn_CNTBBCSFVTXN_cent[ic][j].second);

      if (vn_CNTBBCSFVTXN_cent[ic][j].first > gvn_CNTBBCSFVTXN_mult[j]->GetMaximum())
        gvn_CNTBBCSFVTXN_mult[j]->SetMaximum(vn_CNTBBCSFVTXN_cent[ic][j].first);
      if (vn_CNTBBCSFVTXN_cent[ic][j].first < gvn_CNTBBCSFVTXN_mult[j]->GetMinimum())
        gvn_CNTBBCSFVTXN_mult[j]->SetMinimum(vn_CNTBBCSFVTXN_cent[ic][j].first);

    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_CNTBBCSFVTXN_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gv%i_dAu%i_CNTBBCSFVTXN_pT_c%i", j + 1, energy, ic);
      gvn_CNTBBCSFVTXN_pT[ic][j]->SetName(hname);
      gvn_CNTBBCSFVTXN_pT[ic][j]->SetMarkerStyle(kFullDiamond);
      gvn_CNTBBCSFVTXN_pT[ic][j]->SetMarkerColor(kGreen + 2);
      gvn_CNTBBCSFVTXN_pT[ic][j]->SetLineColor(kGreen + 2);
      gvn_CNTBBCSFVTXN_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {
        vn_CNTBBCSFVTXN_pT[ic][j][ipt] = vn_3sub(
                                           cn_pt[CNTBBCS][ic][ipt][j],
                                           cn_pt[CNTFVTXN][ic][ipt][j],
                                           cn[FVTXNBBCS][ic][j] );

        //fill tgraph
        gvn_CNTBBCSFVTXN_pT[ic][j]->SetPoint(ipt,
                                             0.5 * (ptl[ipt] + pth[ipt]) + ptshift,
                                             vn_CNTBBCSFVTXN_pT[ic][j][ipt].first);
        gvn_CNTBBCSFVTXN_pT[ic][j]->SetPointError(ipt, 0, vn_CNTBBCSFVTXN_pT[ic][j][ipt].second);

        if (vn_CNTBBCSFVTXN_pT[ic][j][ipt].first > gvn_CNTBBCSFVTXN_pT[ic][j]->GetMaximum())
          gvn_CNTBBCSFVTXN_pT[ic][j]->SetMaximum(vn_CNTBBCSFVTXN_pT[ic][j][ipt].first);
        if (vn_CNTBBCSFVTXN_pT[ic][j][ipt].first < gvn_CNTBBCSFVTXN_pT[ic][j]->GetMinimum())
          gvn_CNTBBCSFVTXN_pT[ic][j]->SetMinimum(vn_CNTBBCSFVTXN_pT[ic][j][ipt].first);

      }

    }
  } // j




  //--CNT-FVTXS-BBCN
  for (int j = 0; j < NPAR; j++)
  {
    // centrality dependence
    gvn_CNTFVTXSBBCN_cent[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTFVTXSBBCN_cent", j + 1, energy);
    gvn_CNTFVTXSBBCN_cent[j]->SetName(hname);
    gvn_CNTFVTXSBBCN_cent[j]->SetMarkerStyle(kFullCross);
    gvn_CNTFVTXSBBCN_cent[j]->SetMarkerColor(kMagenta + 2);
    gvn_CNTFVTXSBBCN_cent[j]->SetLineColor(kMagenta + 2);
    gvn_CNTFVTXSBBCN_cent[j]->SetMinimum(100);

    // multiplicity dependence
    gvn_CNTFVTXSBBCN_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTFVTXSBBCN_mult", j + 1, energy);
    gvn_CNTFVTXSBBCN_mult[j]->SetName(hname);
    gvn_CNTFVTXSBBCN_mult[j]->SetMarkerStyle(kFullCross);
    gvn_CNTFVTXSBBCN_mult[j]->SetMarkerColor(kMagenta + 2);
    gvn_CNTFVTXSBBCN_mult[j]->SetLineColor(kMagenta + 2);
    gvn_CNTFVTXSBBCN_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_CNTFVTXSBBCN_cent[ic][j] = vn_3sub(
                                      cn[CNTFVTXS][ic][j],
                                      cn[CNTBBCN][ic][j],
                                      cn[BBCNFVTXS][ic][j] );

      //fill centrality tgraph
      gvn_CNTFVTXSBBCN_cent[j]->SetPoint(ic, ic + 0.5 - 2*cshift, vn_CNTFVTXSBBCN_cent[ic][j].first);
      gvn_CNTFVTXSBBCN_cent[j]->SetPointError(ic, 0, vn_CNTFVTXSBBCN_cent[ic][j].second);

      if (vn_CNTFVTXSBBCN_cent[ic][j].first > gvn_CNTFVTXSBBCN_cent[j]->GetMaximum())
        gvn_CNTFVTXSBBCN_cent[j]->SetMaximum(vn_CNTFVTXSBBCN_cent[ic][j].first);
      if (vn_CNTFVTXSBBCN_cent[ic][j].first < gvn_CNTFVTXSBBCN_cent[j]->GetMinimum())
        gvn_CNTFVTXSBBCN_cent[j]->SetMinimum(vn_CNTFVTXSBBCN_cent[ic][j].first);

      //fill multiplicity tgraph
      gvn_CNTFVTXSBBCN_mult[j]->SetPoint(ic,
                                         mean_npc1[ic].first,
                                         vn_CNTFVTXSBBCN_cent[ic][j].first);
      gvn_CNTFVTXSBBCN_mult[j]->SetPointError(ic, 0, vn_CNTFVTXSBBCN_cent[ic][j].second);

      if (vn_CNTFVTXSBBCN_cent[ic][j].first > gvn_CNTFVTXSBBCN_mult[j]->GetMaximum())
        gvn_CNTFVTXSBBCN_mult[j]->SetMaximum(vn_CNTFVTXSBBCN_cent[ic][j].first);
      if (vn_CNTFVTXSBBCN_cent[ic][j].first < gvn_CNTFVTXSBBCN_mult[j]->GetMinimum())
        gvn_CNTFVTXSBBCN_mult[j]->SetMinimum(vn_CNTFVTXSBBCN_cent[ic][j].first);

    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_CNTFVTXSBBCN_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gv%i_dAu%i_CNTFVTXSBBCN_pT_c%i", j + 1, energy, ic);
      gvn_CNTFVTXSBBCN_pT[ic][j]->SetName(hname);
      gvn_CNTFVTXSBBCN_pT[ic][j]->SetMarkerStyle(kFullCross);
      gvn_CNTFVTXSBBCN_pT[ic][j]->SetMarkerColor(kMagenta + 2);
      gvn_CNTFVTXSBBCN_pT[ic][j]->SetLineColor(kMagenta + 2);
      gvn_CNTFVTXSBBCN_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {
        vn_CNTFVTXSBBCN_pT[ic][j][ipt] = vn_3sub(
                                           cn_pt[CNTFVTXS][ic][ipt][j],
                                           cn_pt[CNTBBCN][ic][ipt][j],
                                           cn[BBCNFVTXS][ic][j] );

        //fill tgraph
        gvn_CNTFVTXSBBCN_pT[ic][j]->SetPoint(ipt,
                                             0.5 * (ptl[ipt] + pth[ipt]) - 2 * ptshift,
                                             vn_CNTFVTXSBBCN_pT[ic][j][ipt].first);
        gvn_CNTFVTXSBBCN_pT[ic][j]->SetPointError(ipt, 0, vn_CNTFVTXSBBCN_pT[ic][j][ipt].second);

        if (vn_CNTFVTXSBBCN_pT[ic][j][ipt].first > gvn_CNTFVTXSBBCN_pT[ic][j]->GetMaximum())
          gvn_CNTFVTXSBBCN_pT[ic][j]->SetMaximum(vn_CNTFVTXSBBCN_pT[ic][j][ipt].first);
        if (vn_CNTFVTXSBBCN_pT[ic][j][ipt].first < gvn_CNTFVTXSBBCN_pT[ic][j]->GetMinimum())
          gvn_CNTFVTXSBBCN_pT[ic][j]->SetMinimum(vn_CNTFVTXSBBCN_pT[ic][j][ipt].first);

      }

    }
  } // j




  //--CNT-BBCS-FVTXN
  for (int j = 0; j < NPAR; j++)
  {
    // centrality dependence
    gvn_CNTBBCSBBCN_cent[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSBBCN_cent", j + 1, energy);
    gvn_CNTBBCSBBCN_cent[j]->SetName(hname);
    gvn_CNTBBCSBBCN_cent[j]->SetMarkerStyle(kOpenCircle);
    gvn_CNTBBCSBBCN_cent[j]->SetMarkerColor(kOrange + 2);
    gvn_CNTBBCSBBCN_cent[j]->SetLineColor(kOrange + 2);
    gvn_CNTBBCSBBCN_cent[j]->SetMinimum(100);

    // multiplicity dependence
    gvn_CNTBBCSBBCN_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_dAu%i_CNTBBCSBBCN_mult", j + 1, energy);
    gvn_CNTBBCSBBCN_mult[j]->SetName(hname);
    gvn_CNTBBCSBBCN_mult[j]->SetMarkerStyle(kOpenCircle);
    gvn_CNTBBCSBBCN_mult[j]->SetMarkerColor(kOrange + 2);
    gvn_CNTBBCSBBCN_mult[j]->SetLineColor(kOrange + 2);
    gvn_CNTBBCSBBCN_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_CNTBBCSBBCN_cent[ic][j] = vn_3sub(
                                     cn[CNTBBCS][ic][j],
                                     cn[CNTBBCN][ic][j],
                                     cn[BBCNBBCS][ic][j] );

      //fill centrality tgraph
      gvn_CNTBBCSBBCN_cent[j]->SetPoint(ic, ic + 0.5 + 2*cshift, vn_CNTBBCSBBCN_cent[ic][j].first);
      gvn_CNTBBCSBBCN_cent[j]->SetPointError(ic, 0, vn_CNTBBCSBBCN_cent[ic][j].second);

      if (vn_CNTBBCSBBCN_cent[ic][j].first > gvn_CNTBBCSBBCN_cent[j]->GetMaximum())
        gvn_CNTBBCSBBCN_cent[j]->SetMaximum(vn_CNTBBCSBBCN_cent[ic][j].first);
      if (vn_CNTBBCSBBCN_cent[ic][j].first < gvn_CNTBBCSBBCN_cent[j]->GetMinimum())
        gvn_CNTBBCSBBCN_cent[j]->SetMinimum(vn_CNTBBCSBBCN_cent[ic][j].first);

      //fill multiplicity tgraph
      gvn_CNTBBCSBBCN_mult[j]->SetPoint(ic,
                                        mean_npc1[ic].first,
                                        vn_CNTBBCSBBCN_cent[ic][j].first);
      gvn_CNTBBCSBBCN_mult[j]->SetPointError(ic, 0, vn_CNTBBCSBBCN_cent[ic][j].second);

      if (vn_CNTBBCSBBCN_cent[ic][j].first > gvn_CNTBBCSBBCN_mult[j]->GetMaximum())
        gvn_CNTBBCSBBCN_mult[j]->SetMaximum(vn_CNTBBCSBBCN_cent[ic][j].first);
      if (vn_CNTBBCSBBCN_cent[ic][j].first < gvn_CNTBBCSBBCN_mult[j]->GetMinimum())
        gvn_CNTBBCSBBCN_mult[j]->SetMinimum(vn_CNTBBCSBBCN_cent[ic][j].first);

    }

    //pT dependence
    for (int ic = 0; ic < NC; ic++)
    {
      gvn_CNTBBCSBBCN_pT[ic][j] = new TGraphErrors();
      sprintf(hname, "gv%i_dAu%i_CNTBBCSBBCN_pT_c%i", j + 1, energy, ic);
      gvn_CNTBBCSBBCN_pT[ic][j]->SetName(hname);
      gvn_CNTBBCSBBCN_pT[ic][j]->SetMarkerStyle(kOpenCircle);
      gvn_CNTBBCSBBCN_pT[ic][j]->SetMarkerColor(kOrange + 2);
      gvn_CNTBBCSBBCN_pT[ic][j]->SetLineColor(kOrange + 2);
      gvn_CNTBBCSBBCN_pT[ic][j]->SetMinimum(100);

      for (int ipt = 0; ipt < NPT; ipt++)
      {
        vn_CNTBBCSBBCN_pT[ic][j][ipt] = vn_3sub(
                                          cn_pt[CNTBBCS][ic][ipt][j],
                                          cn_pt[CNTBBCN][ic][ipt][j],
                                          cn[BBCNBBCS][ic][j] );

        //fill tgraph
        gvn_CNTBBCSBBCN_pT[ic][j]->SetPoint(ipt,
                                            0.5 * (ptl[ipt] + pth[ipt]) + 2* ptshift,
                                            vn_CNTBBCSBBCN_pT[ic][j][ipt].first);
        gvn_CNTBBCSBBCN_pT[ic][j]->SetPointError(ipt, 0, vn_CNTBBCSBBCN_pT[ic][j][ipt].second);

        if (vn_CNTBBCSBBCN_pT[ic][j][ipt].first > gvn_CNTBBCSBBCN_pT[ic][j]->GetMaximum())
          gvn_CNTBBCSBBCN_pT[ic][j]->SetMaximum(vn_CNTBBCSBBCN_pT[ic][j][ipt].first);
        if (vn_CNTBBCSBBCN_pT[ic][j][ipt].first < gvn_CNTBBCSBBCN_pT[ic][j]->GetMinimum())
          gvn_CNTBBCSBBCN_pT[ic][j]->SetMinimum(vn_CNTBBCSBBCN_pT[ic][j][ipt].first);

      }

    }
  } // j




  //--CNT-BBCS-FVTXN (numerical C_n's)
  for (int j = 0; j < NPAR; j++)
  {
    // multiplicity dependence
    gvn_num_CNTBBCSFVTXN_mult[j] = new TGraphErrors();
    sprintf(hname, "gv%i_num_dAu%i_CNTBBCSFVTXN_mult", j + 1, energy);
    gvn_num_CNTBBCSFVTXN_mult[j]->SetName(hname);
    gvn_num_CNTBBCSFVTXN_mult[j]->SetMarkerStyle(kOpenDiamond);
    gvn_num_CNTBBCSFVTXN_mult[j]->SetMarkerColor(kGreen + 2);
    gvn_num_CNTBBCSFVTXN_mult[j]->SetLineColor(kGreen + 2);
    gvn_num_CNTBBCSFVTXN_mult[j]->SetMinimum(100);

    for (int ic = 0; ic < NC; ic++)
    {

      vn_num_CNTBBCSFVTXN_cent[ic][j] = vn_3sub(
                                          cn_num[CNTBBCS][ic][j],
                                          cn_num[CNTFVTXN][ic][j],
                                          cn_num[FVTXNBBCS][ic][j] );

      //fill multiplicity tgraph
      gvn_num_CNTBBCSFVTXN_mult[j]->SetPoint(ic,
                                             mean_npc1[ic].first + 0.05,
                                             vn_num_CNTBBCSFVTXN_cent[ic][j].first);
      gvn_num_CNTBBCSFVTXN_mult[j]->SetPointError(ic, 0, vn_num_CNTBBCSFVTXN_cent[ic][j].second);

      if (vn_num_CNTBBCSFVTXN_cent[ic][j].first > gvn_num_CNTBBCSFVTXN_mult[j]->GetMaximum())
        gvn_num_CNTBBCSFVTXN_mult[j]->SetMaximum(vn_num_CNTBBCSFVTXN_cent[ic][j].first);
      if (vn_num_CNTBBCSFVTXN_cent[ic][j].first < gvn_num_CNTBBCSFVTXN_mult[j]->GetMinimum())
        gvn_num_CNTBBCSFVTXN_mult[j]->SetMinimum(vn_num_CNTBBCSFVTXN_cent[ic][j].first);

    } // ic
  } // j

  //==========================================================================//
  // GET RON'S EP METHOD, FIT IT, AND TAKE RATIOS TO FIT
  //==========================================================================//
  if (energy > 20)
  {
    //-- Difine TGraph
    gv2_Ron = new TGraphErrors();
    gv2_Ron->SetMarkerStyle(kFullCircle);
    gv2_Ron->SetMarkerColor(kBlack);
    gv2_Ron->SetLineColor(kBlack);


    if ( energy == 200 )
    {
      for (int i = 0; i < Ron_v2_dAu200_N; i++)
      {
        gv2_Ron->SetPoint(i, Ron_v2_dAu200_pT[i], Ron_v2_dAu200_v2[i]);
        gv2_Ron->SetPointError(i, 0, Ron_v2_dAu200_v2e[i]);
      }
    }
    if ( energy == 62 )
    {
      for (int i = 0; i < Ron_v2_dAu62_N; i++)
      {
        gv2_Ron->SetPoint(i, Ron_v2_dAu62_pT[i], Ron_v2_dAu62_v2[i]);
        gv2_Ron->SetPointError(i, 0, Ron_v2_dAu62_v2e[i]);
      }
    }
    if ( energy == 39 )
    {
      for (int i = 0; i < Ron_v2_dAu39_N; i++)
      {
        gv2_Ron->SetPoint(i, Ron_v2_dAu39_pT[i], Ron_v2_dAu39_v2[i]);
        gv2_Ron->SetPointError(i, 0, Ron_v2_dAu39_v2e[i]);
      }
    }

    //-- Fit
    gv2_Ron->Fit(fpol, "RQ0N");


    //-- Ratios
    double x, y;

    gv2fitrat_Ron = (TGraphErrors*) gv2_Ron->Clone("gv2fitrat_Ron");
    for (int i = 0; i < gv2fitrat_Ron->GetN(); i++)
    {
      gv2fitrat_Ron->GetPoint(i, x, y);
      gv2fitrat_Ron->SetPoint(i, x,
                              y / fpol->Eval(x));
      gv2fitrat_Ron->SetPointError(i, 0,
                                   gv2fitrat_Ron->GetErrorY(i) / fpol->Eval(x));
    }


    gv2fitrat_CNTFVTXSFVTXN =
      (TGraphErrors*) gvn_CNTFVTXSFVTXN_pT[0][1]->Clone(
        "gv2fitrat_CNTFVTXSFVTXN");
    for (int i = 0; i < gv2fitrat_CNTFVTXSFVTXN->GetN(); i++)
    {
      gv2fitrat_CNTFVTXSFVTXN->GetPoint(i, x, y);
      y = y / fpol->Eval(x);
      double ey = gv2fitrat_CNTFVTXSFVTXN->GetErrorY(i) / fpol->Eval(x);

      gv2fitrat_CNTFVTXSFVTXN->SetPoint(i, x, y);
      gv2fitrat_CNTFVTXSFVTXN->SetPointError(i, 0, ey);
    }


    gv2fitrat_CNTBBCSFVTXS =
      (TGraphErrors*) gvn_CNTBBCSFVTXS_pT[0][1]->Clone(
        "gv2fitrat_CNTBBCSFVTXS");
    for (int i = 0; i < gv2fitrat_CNTBBCSFVTXS->GetN(); i++)
    {
      gv2fitrat_CNTBBCSFVTXS->GetPoint(i, x, y);
      y = y / fpol->Eval(x);
      double ey = gv2fitrat_CNTBBCSFVTXS->GetErrorY(i) / fpol->Eval(x);

      gv2fitrat_CNTBBCSFVTXS->SetPoint(i, x, y);
      gv2fitrat_CNTBBCSFVTXS->SetPointError(i, 0, ey);
    }


    gv2fitrat_CNTFVTXSBBCN =
      (TGraphErrors*) gvn_CNTFVTXSBBCN_pT[0][1]->Clone(
        "gv2fitrat_CNTFVTXSBBCN");
    for (int i = 0; i < gv2fitrat_CNTFVTXSBBCN->GetN(); i++)
    {
      gv2fitrat_CNTFVTXSBBCN->GetPoint(i, x, y);
      y = y / fpol->Eval(x);
      double ey = gv2fitrat_CNTFVTXSBBCN->GetErrorY(i) / fpol->Eval(x);

      gv2fitrat_CNTFVTXSBBCN->SetPoint(i, x, y);
      gv2fitrat_CNTFVTXSBBCN->SetPointError(i, 0, ey);
    }

    gv2fitrat_CNTBBCSFVTXN =
      (TGraphErrors*) gvn_CNTBBCSFVTXN_pT[0][1]->Clone(
        "gv2fitrat_CNTBBCSFVTXN");
    for (int i = 0; i < gv2fitrat_CNTBBCSFVTXN->GetN(); i++)
    {
      gv2fitrat_CNTBBCSFVTXN->GetPoint(i, x, y);
      y = y / fpol->Eval(x);
      double ey = gv2fitrat_CNTBBCSFVTXN->GetErrorY(i) / fpol->Eval(x);

      gv2fitrat_CNTBBCSFVTXN->SetPoint(i, x, y);
      gv2fitrat_CNTBBCSFVTXN->SetPointError(i, 0, ey);
    }

    gv2fitrat_CNTBBCSBBCN =
      (TGraphErrors*) gvn_CNTBBCSBBCN_pT[0][1]->Clone(
        "gv2fitrat_CNTBBCSBBCN");
    for (int i = 0; i < gv2fitrat_CNTBBCSBBCN->GetN(); i++)
    {
      gv2fitrat_CNTBBCSBBCN->GetPoint(i, x, y);
      y = y / fpol->Eval(x);
      double ey = gv2fitrat_CNTBBCSBBCN->GetErrorY(i) / fpol->Eval(x);

      gv2fitrat_CNTBBCSBBCN->SetPoint(i, x, y);
      gv2fitrat_CNTBBCSBBCN->SetPointError(i, 0, ey);
    }






  }


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

  TLegend *leg3sub = new TLegend(0.15, 0.6, 0.6, 0.88);
  leg3sub->SetFillStyle(0);
  leg3sub->SetBorderSize(0);
  leg3sub->AddEntry(gvn_CNTFVTXSFVTXN_pT[0][0], "CNT--FVTXN--FVTXS", "P");
  leg3sub->AddEntry(gvn_CNTBBCSFVTXN_pT[0][0], "CNT--FVTXN--BBCS", "P");
  leg3sub->AddEntry(gvn_CNTFVTXSBBCN_pT[0][0], "CNT--BBCN--FVTXS", "P");
  leg3sub->AddEntry(gvn_CNTBBCSBBCN_pT[0][0], "CNT--BBCN--BBCS", "P");
  leg3sub->AddEntry(gvn_CNTBBCSFVTXS_pT[0][0], "CNT--BBCS--FVTXS", "P");


  TLegend *legv2 = new TLegend(0.15, 0.6, 0.6, 0.88);
  legv2->SetFillStyle(0);
  legv2->SetBorderSize(0);
  legv2->AddEntry(gvn_CNTFVTXSFVTXN_pT[0][0], "3 sub-event CNT-FVTXN-FVTXS", "P");
  legv2->AddEntry(gvn_CNTBBCSFVTXN_pT[0][0], "3 sub-event CNT-FVTXN-BBCS", "P");
  if (energy == 200)
  {
    legv2->AddEntry(gvn_CNTFVTXSBBCN_pT[0][0], "3 sub-event CNT-BBCN-FVTXS", "P");
    legv2->AddEntry(gvn_CNTBBCSBBCN_pT[0][0], "3 sub-event CNT-BBCN-BBCS", "P");
  }
  legv2->AddEntry(gvn_CNTBBCSFVTXS_pT[0][0], "3 sub-event CNT-BBCS-FVTXS", "P");
  if (energy > 20)
    legv2->AddEntry(gv2_Ron, "FVTXS EP (Ron)", "P");


  TH1F* haxis_cncent = new TH1F("haxis_cncent", "", NC, 0, NC);
  haxis_cncent->GetYaxis()->CenterTitle();
  haxis_cncent->GetYaxis()->SetTitleFont(63);
  haxis_cncent->GetYaxis()->SetTitleSize(20);
  haxis_cncent->GetYaxis()->SetTitleOffset(1.7);
  haxis_cncent->GetYaxis()->SetLabelFont(63);
  haxis_cncent->GetYaxis()->SetLabelSize(14);
  haxis_cncent->GetXaxis()->SetTitleFont(63);
  haxis_cncent->GetXaxis()->SetTitleSize(20);
  haxis_cncent->GetXaxis()->SetTitleOffset(1.7);
  haxis_cncent->GetXaxis()->SetLabelFont(63);
  haxis_cncent->GetXaxis()->SetLabelSize(14);

  for (int ibin = 1; ibin <= haxis_cncent->GetNbinsX(); ibin++)
  {
    sprintf(ctitle, "%i - %i%%", cl[ibin - 1], ch[ibin - 1]);
    haxis_cncent->GetXaxis()->SetBinLabel(ibin, ctitle);
  }


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

  TH1F* haxis_ratpt = new TH1F("haxis_ratpt", ";p_{T} [GeV/c];v_{2} / fit", 100, 0, 5);
  haxis_ratpt->SetMinimum(0.01);
  haxis_ratpt->SetMaximum(0.29);
  haxis_ratpt->GetYaxis()->CenterTitle();
  haxis_ratpt->GetYaxis()->SetTitleFont(63);
  haxis_ratpt->GetYaxis()->SetTitleSize(20);
  haxis_ratpt->GetYaxis()->SetTitleOffset(1.7);
  haxis_ratpt->GetYaxis()->SetLabelFont(63);
  haxis_ratpt->GetYaxis()->SetLabelSize(14);
  haxis_ratpt->GetXaxis()->SetTitleFont(63);
  haxis_ratpt->GetXaxis()->SetTitleSize(20);
  haxis_ratpt->GetXaxis()->SetTitleOffset(1.7);
  haxis_ratpt->GetXaxis()->SetLabelFont(63);
  haxis_ratpt->GetXaxis()->SetLabelSize(14);




  float min, max;

  TH1D* hrat;



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
              fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i].first);
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
            fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i].first);
          fcorr->DrawCopy("same");

          for (int ipar = 0; ipar < NPAR; ipar++)
          {
            fcorr->SetLineColor(fitColor[ipar]);
            for (int i = 0; i < NPAR; i++)
            {
              if (i == ipar)
                fcorr->SetParameter(i, cn_pt[icorr][ic][ipt][i].first);
              else
                fcorr->SetParameter(i, 0);
            }
            fcorr->DrawCopy("same");
          }


        lc.SetTextColor(fitColor[0]);
        lc.DrawLatex(0.45, 0.7, Form("C_{1} = % .4f",
                                     cn_pt[icorr][ic][ipt][0].first));
        lc.SetTextColor(fitColor[1]);
        lc.DrawLatex(0.45, 0.64, Form("C_{2} = % .4f",
                                      cn_pt[icorr][ic][ipt][1].first));
        lc.SetTextColor(kBlack);
        lc.DrawLatex(0.45, 0.58, Form("C_{2}/C_{1} = % .4f",
                                      cn_pt[icorr][ic][ipt][1].first / cn_pt[icorr][ic][ipt][0].first));


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
                                       cn[icorr][ic][0].first));

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
            fcorr->SetParameter(i, cn[icorr][ic][i].first);
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
          fcorr->SetParameter(i, cn[icorr][ic][i].first);
        fcorr->DrawCopy("same");

        for (int ipar = 0; ipar < NPAR; ipar++)
        {
          fcorr->SetLineColor(fitColor[ipar]);
          for (int i = 0; i < NPAR; i++)
          {
            if (i == ipar)
              fcorr->SetParameter(i, cn[icorr][ic][i].first);
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
                                     cn[icorr][ic][0].first));
        lc.SetTextColor(fitColor[1]);
        lc.DrawLatex(0.45, 0.64, Form("C_{2} = % .4f",
                                      cn[icorr][ic][1].first));
        lc.SetTextColor(kBlack);
        lc.DrawLatex(0.45, 0.58, Form("C_{2}/C_{1} = % .4f",
                                      cn[icorr][ic][1].first / cn[icorr][ic][0].first));


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

    if (hi > max)
      max = hi;
    if (lo < min)
      min = lo;
  }
  min = (min < 0) ? min : 0; // set the min to 0 unless it's negative
  haxis_cncent->SetMinimum(1.1 * min);
  haxis_cncent->SetMaximum(0.9 * max);
  haxis_cncent->SetTitle(";;-C_{2} / C_{1}");
  haxis_cncent->DrawCopy();

  for (vector<int>::size_type j = 0; j < corrdrawidx.size(); j++)
    gc2c1[corrdrawidx.at(j)]->Draw("P");





  TCanvas *cvncent = new TCanvas("cvncent", "vn cent", 800, 600);
  cvncent->Divide(1, 2);

  for (int j = 1; j < 3; j++)
  {
    cvncent->cd(j);
    haxis_cncent->SetMinimum(0.0);
    haxis_cncent->SetMaximum(0.24);
    haxis_cncent->SetTitle(Form(";;v_{%i}", j + 1));
    haxis_cncent->DrawCopy();

    gvn_CNTFVTXSFVTXN_cent[j]->Draw("P");
    gvn_CNTBBCSFVTXS_cent[j]->Draw("P");
    gvn_CNTBBCSFVTXN_cent[j]->Draw("P");

    if (j == 1)
    {
      ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV", energy));
      le.DrawLatex(0.15, 0.86, "3 sub-event #color[4]{CNT}-FVTXN-FVTXS");
      le.DrawLatex(0.15, 0.80, "3 sub-event #color[2]{CNT}-BBCS-FVTXS");
      le.DrawLatex(0.15, 0.74, "3 sub-event #color[418]{CNT}-BBCS-FVTXN");
      le.DrawLatex(0.15, 0.68, Form("%.1f<p_{T}^{trig} [GeV/c]<%.1f",
                                    ptl[ptsuml], pth[ptsumh]));
    }
  }


  TCanvas *cvnpt = new TCanvas("cvnpt", "vn pt", 800, 600);
  cvnpt->Divide(1, 2);

  for (int j = 1; j < 3; j++)
  {

    cvnpt->cd(j);
    gvn_CNTFVTXSFVTXN_pT[0][j]->SetMinimum(0);
    gvn_CNTFVTXSFVTXN_pT[0][j]->SetMaximum(0.24);
    gvn_CNTFVTXSFVTXN_pT[0][j]->GetXaxis()->SetRangeUser(0, 5);
    gvn_CNTFVTXSFVTXN_pT[0][j]->SetTitle(Form(";p_{T} [GeV/c];v_{%i}", j + 1));
    gvn_CNTFVTXSFVTXN_pT[0][j]->Draw("AP");

    gvn_CNTBBCSFVTXN_pT[0][j]->Draw("P");
    gvn_CNTBBCSFVTXS_pT[0][j]->Draw("P");

    if (j == 1)
    {
      ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));
      le.DrawLatex(0.2, 0.86, "3 sub-event #color[4]{CNT}-FVTXN-FVTXS");
      le.DrawLatex(0.2, 0.80, "3 sub-event #color[2]{CNT}-BBCS-FVTXS");
      le.DrawLatex(0.2, 0.74, "3 sub-event #color[418]{CNT}-BBCS-FVTXN");

    }
  }


  TCanvas *cv2 = new TCanvas("cv2", "v2", 600, 800);
  cv2->Divide(1, 2);

  cv2->cd(1);
  haxis_vnpt->GetYaxis()->SetRangeUser(0, 0.29);
  haxis_vnpt->GetXaxis()->SetRangeUser(0, 5);
  haxis_vnpt->DrawCopy();

  gvn_CNTFVTXSFVTXN_pT[0][1]->Draw("P");
  gvn_CNTBBCSFVTXN_pT[0][1]->Draw("P");
  gvn_CNTBBCSFVTXS_pT[0][1]->Draw("P");
  gvn_CNTFVTXSBBCN_pT[0][1]->Draw("P");
  gvn_CNTBBCSBBCN_pT[0][1]->Draw("P");

  ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));
  leg3sub->Draw("same");



  cv2->cd(2);
  cv2->GetPad(2)->SetTopMargin(0.02);
  haxis_cncent->SetMinimum(0.0);
  haxis_cncent->SetMaximum(0.26);
  haxis_cncent->SetTitle(";;v_{2}");
  haxis_cncent->DrawCopy();

  gvn_CNTFVTXSFVTXN_cent[1]->Draw("P");
  gvn_CNTBBCSFVTXS_cent[1]->Draw("P");
  gvn_CNTBBCSFVTXN_cent[1]->Draw("P");
  gvn_CNTFVTXSBBCN_cent[1]->Draw("P");
  gvn_CNTBBCSBBCN_cent[1]->Draw("P");

  le.DrawLatex(0.15, 0.86, Form("%.1f<p_{T}^{trig} [GeV/c]<%.1f",
                                ptl[ptsuml], pth[ptsumh]));





  TCanvas *cv2multnum = new TCanvas("cv2multnum", "v2 mult num", 800, 600);

  cv2multnum->cd(1);
  gvn_CNTBBCSFVTXN_mult[1]->GetYaxis()->SetRangeUser(0, 0.30);
  gvn_CNTBBCSFVTXN_mult[1]->Draw("AP");
  gvn_num_CNTBBCSFVTXN_mult[1]->Draw("P");


  TCanvas* cdeta = new TCanvas("cdeta", "deta", 800, 600);

  cdeta->cd(1);
  if (gc2c1_deta[0]->GetMinimum() > 0)
    gc2c1_deta[0]->SetMinimum(0);
  gc2c1_deta[0]->SetMaximum(1.1 * gc2c1_deta[0]->GetMaximum());

  gc2c1_deta[0]->Draw("AP");

  ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));





  TCanvas *cv2ptfitrat;
  if (energy > 20)
  {
    cv2ptfitrat = new TCanvas("cv2ptfitrat", "v2 pT fit rat", 600, 700);
    cv2ptfitrat->SetTopMargin(0);
    cv2ptfitrat->SetRightMargin(0);
    cv2ptfitrat->SetBottomMargin(0);
    cv2ptfitrat->SetLeftMargin(0);

    //-- v2
    // cv2ptfitrat->cd(1);
    TPad *pv2 = new TPad("pv2", "v2", 0, 0.4, 1, 1);
    pv2->SetTopMargin(0.10);
    pv2->SetRightMargin(0.02);
    pv2->SetBottomMargin(0.0);
    pv2->SetLeftMargin(0.12);
    pv2->SetTicks(1, 1);
    pv2->Draw();
    pv2->cd();

    haxis_vnpt->GetXaxis()->SetRangeUser(0, 3);
    haxis_vnpt->DrawCopy();

    gvn_CNTFVTXSFVTXN_pT[0][1]->Draw("P");
    gvn_CNTBBCSFVTXN_pT[0][1]->Draw("P");
    gvn_CNTBBCSFVTXS_pT[0][1]->Draw("P");
    if (energy == 200)
    {
      gvn_CNTFVTXSBBCN_pT[0][1]->Draw("P");
      gvn_CNTBBCSBBCN_pT[0][1]->Draw("P");
    }

    gv2_Ron->Draw("P");
    fpol->Draw("same");

    ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));
    legv2->Draw("same");

    //-- Ratio to fit
    // cv2ptfitrat->cd(2);
    cv2ptfitrat->cd();
    TPad *prat = new TPad("prat", "rat", 0, 0, 1, 0.4);
    prat->SetTopMargin(0.0);
    prat->SetRightMargin(0.02);
    prat->SetBottomMargin(0.18);
    prat->SetLeftMargin(0.12);
    prat->SetTicks(1, 1);
    prat->Draw();
    prat->cd();

    haxis_ratpt->GetYaxis()->SetRangeUser(0.51, 1.49);
    haxis_ratpt->GetXaxis()->SetRangeUser(0, 3);
    haxis_ratpt->Draw();

    gv2fitrat_CNTFVTXSFVTXN->Draw("P");
    gv2fitrat_CNTBBCSFVTXN->Draw("P");
    gv2fitrat_CNTBBCSFVTXS->Draw("P");
    if (energy == 200)
    {
      gv2fitrat_CNTFVTXSBBCN->Draw("P");
      gv2fitrat_CNTBBCSBBCN->Draw("P");
    }

    gv2fitrat_Ron->Draw("P");

    l1.DrawLine(0, 1, 3, 1);

  }


  //-- v2 vs pT in centrality bins
  TCanvas* cv2ptcent = new TCanvas("cv2ptcent", "v2 pt centrality", 1200, 500);
  cv2ptcent->SetMargin(0, 0, 0, 0);
  cv2ptcent->Divide(3, 1, 0, 0);

  cv2ptcent->GetPad(1)->SetMargin(0.12, 0.00, 0.12, 0.08);
  cv2ptcent->GetPad(2)->SetMargin(0.00, 0.00, 0.12, 0.08);
  cv2ptcent->GetPad(3)->SetMargin(0.00, 0.02, 0.12, 0.08);

  cv2ptcent->GetPad(1)->SetTicks(1, 1);
  cv2ptcent->GetPad(2)->SetTicks(1, 1);
  cv2ptcent->GetPad(3)->SetTicks(1, 1);

  cv2ptcent->cd(1);
  haxis_vnpt->GetYaxis()->SetRangeUser(0.0, 0.34);
  haxis_vnpt->GetXaxis()->SetRangeUser(0.1, 2.9);
  haxis_vnpt->GetXaxis()->SetTitleOffset(1.2);
  haxis_vnpt->GetYaxis()->SetTitleOffset(1.5);
  haxis_vnpt->DrawCopy();

  TGraphErrors *gtmp_CNTFVTXSFVTXN =
    (TGraphErrors*) gvn_CNTFVTXSFVTXN_pT[0][1]->Clone("gtmp_CNTFVTXSFVTXN");
  gtmp_CNTFVTXSFVTXN->SetMarkerColor(centColor[0]);
  gtmp_CNTFVTXSFVTXN->SetLineColor(centColor[0]);
  gtmp_CNTFVTXSFVTXN->Draw("P");


  for (int ic = 1; ic < NC; ic++)
  {
    gvn_CNTFVTXSFVTXN_pT[ic][1]->SetMarkerColor(centColor[ic]);
    gvn_CNTFVTXSFVTXN_pT[ic][1]->SetLineColor(centColor[ic]);
    gvn_CNTFVTXSFVTXN_pT[ic][1]->Draw("P");
  }

  le.DrawLatex(0.2, 0.84, "3 sub-event CNT-FVTXN-FVTXS");

  TLegend *legc = new TLegend(0.7, 0.15, 0.98, 0.45);
  legc->SetFillStyle(0);
  legc->SetBorderSize(0);
  legc->AddEntry(gtmp_CNTFVTXSFVTXN, Form("%i--%i%%", cl[0], ch[0]), "P");
  for (int ic = 1; ic < NC; ic++)
    legc->AddEntry(gvn_CNTFVTXSFVTXN_pT[ic][1], Form("%i--%i%%", cl[ic], ch[ic]), "P");
  legc->Draw("same");

  cv2ptcent->cd(2);
  haxis_vnpt->DrawCopy();

  TGraphErrors *gtmp_CNTBBCSFVTXN =
    (TGraphErrors*) gvn_CNTBBCSFVTXN_pT[0][1]->Clone("gtmp_CNTBBCSFVTXN");
  gtmp_CNTBBCSFVTXN->SetMarkerColor(centColor[0]);
  gtmp_CNTBBCSFVTXN->SetLineColor(centColor[0]);
  gtmp_CNTBBCSFVTXN->Draw("P");

  for (int ic = 1; ic < NC; ic++)
  {
    gvn_CNTBBCSFVTXN_pT[ic][1]->SetMarkerColor(centColor[ic]);
    gvn_CNTBBCSFVTXN_pT[ic][1]->SetLineColor(centColor[ic]);
    gvn_CNTBBCSFVTXN_pT[ic][1]->Draw("P");
  }

  ltitle.DrawLatex(0.5, 0.95, Form("d+Au #sqrt{s_{_{NN}}}=%i GeV %i-%i%%", energy, cl[0], ch[0]));
  le.DrawLatex(0.2, 0.84, "3 sub-event CNT-FVTXN-BBCS");


  cv2ptcent->cd(3);
  haxis_vnpt->DrawCopy();

  TGraphErrors *gtmp_CNTBBCSFVTXS =
    (TGraphErrors*) gvn_CNTBBCSFVTXS_pT[0][1]->Clone("gtmp_CNTBBCSFVTXS");
  gtmp_CNTBBCSFVTXS->SetMarkerColor(centColor[0]);
  gtmp_CNTBBCSFVTXS->SetLineColor(centColor[0]);
  gtmp_CNTBBCSFVTXS->Draw("P");

  for (int ic = 1; ic < NC; ic++)
  {
    gvn_CNTBBCSFVTXS_pT[ic][1]->SetMarkerColor(centColor[ic]);
    gvn_CNTBBCSFVTXS_pT[ic][1]->SetLineColor(centColor[ic]);
    gvn_CNTBBCSFVTXS_pT[ic][1]->Draw("P");
  }

  le.DrawLatex(0.2, 0.84, "3 sub-event CNT-FVTXS-BBCS");


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

    sprintf(cname, "pdfs/dAu%i_vncent.pdf", energy);
    cvncent->Print(cname);

    sprintf(cname, "pdfs/dAu%i_vnpT.pdf", energy);
    cvnpt->Print(cname);

    sprintf(cname, "pdfs/dAu%i_v2.pdf", energy);
    cv2->Print(cname);

    sprintf(cname, "pdfs/dAu%i_v2pT_fitrat.pdf", energy);
    cv2ptfitrat->Print(cname);

    sprintf(cname, "pdfs/dAu%i_v2pT_cent.pdf", energy);
    cv2ptcent->Print(cname);

  } // printPlots

  //==========================================================================//
  // SAVE RESULTS TO FILE
  //==========================================================================//
  if (saveHistos)
  {
    cout << endl;
    cout << "--> Saving results to " << outFile << endl;

    TFile *fout = new TFile(outFile, "UPDATE");



    for (int icorr = 0; icorr < NCORR; icorr++)
    {
      for (int ic = 0; ic < NC; ic++)
      {
        //-- pT dependent
        if (icorr < NCORRPT)
        {
          for (int ipt = 0; ipt < NPT; ipt++)
          {
            sprintf(hname, "dphi_corr_dAu%i_%s_c%i_pt%i", energy, cCorr[icorr], ic, ipt);
            dphi_corr_pt[icorr][ic][ipt]->Write(hname);

            for (int ipar = 0; ipar < NPAR; ipar++)
            {
              sprintf(hname, "c%i_dAu%i_%s_c%i_pt%i", ipar, energy, cCorr[icorr], ic, ipt);
              gcn_pt[icorr][ic][ipar]->Write(hname);

            }
          } // ipt

        }

        sprintf(hname, "dphi_corr_dAu%i_%s_c%i", energy, cCorr[icorr], ic);
        dphi_corr[icorr][ic]->Write(hname);

        for (int ipar = 0; ipar < NPAR; ipar++)
        {
          sprintf(hname, "c%i_dAu%i_%s_c%i", ipar, energy, cCorr[icorr], ic);
          gcn[icorr][ipar]->Write(hname);
        }

      } // ic
    } // icorr

    //-- v_n for 3 sub-event method
    for (int ipar = 0; ipar < NPAR; ipar++)
    {
      for (int ic = 0; ic < NC; ic++)
      {
        gvn_CNTFVTXSFVTXN_pT[ic][ipar]->Write();
        gvn_CNTBBCSFVTXS_pT[ic][ipar]->Write();
        gvn_CNTBBCSFVTXN_pT[ic][ipar]->Write();
      }

      gvn_CNTFVTXSFVTXN_cent[ipar]->Write();
      gvn_CNTBBCSFVTXS_cent[ipar]->Write();
      gvn_CNTBBCSFVTXN_cent[ipar]->Write();

      gvn_CNTFVTXSFVTXN_mult[ipar]->Write();
      gvn_CNTBBCSFVTXS_mult[ipar]->Write();
      gvn_CNTBBCSFVTXN_mult[ipar]->Write();
    }

    //-- c2/c1
    sprintf(hname, "gc2c1_deta_dAu%i_c0", energy);
    gc2c1_deta[0]->SetName(hname);
    gc2c1_deta[0]->Write();


    fout->Close();
    delete fout;

  }

  //==========================================================================//
  // PRINT TABLES
  //==========================================================================//
  if (printTables)
  {
    cout << endl;
    cout << "--> Printing Tables" << endl;

    //-- pT dependent central
    cout << endl;
    cout << setw(5) << "$\\sqrt{s_{_{NN}}}$" << " & "
         << setw(10) << "Correlation" << " & "
         << setw(10) << "Centrality" << " & "
         << setw(12) << "$p_T$ [GeV/c]" << " & "
         << setw(8) << "$C_1$" << " & "
         << setw(8) << "$C_2$" << " & "
         << setw(8) << "$C_3$" << " \\\\"
         << endl;
    cout << "\\hline" << endl;
    for (int icorr = 0; icorr < NCORRPT; icorr++)
    {

      for (int ipt = 0; ipt < NPT; ipt++)
      {
        cout << setw(5) << energy << " & "
             << setw(10) << cCorr[icorr] << " & "
             << setw(10) << Form("%i--%i\\%%", cl[0], ch[0]) << " & "
             << setw(12) << Form("[%.2f, %.2f]", ptl[ipt], pth[ipt]);
        for (int ipar = 0; ipar < 3; ipar++)
          cout << " & " << setw(8) << Form("%.5f$\\pm$%.5f", cn_pt[icorr][0][ipt][ipar].first, cn_pt[icorr][0][ipt][ipar].second);
        // cout << " & " << setw(8) << Form("%.5f$\\pm$%.5f (%.5f$\\pm$%.5f)", cn_pt[icorr][0][ipt][ipar].first, cn_pt[icorr][0][ipt][ipar].second, cn_num_pt[icorr][0][ipt][ipar].first, cn_num_pt[icorr][0][ipt][ipar].second);
        cout << " \\\\" << endl;
      } // ipt
      cout << "\\hline" << endl;
    }// icorr


    //-- Centrality dependent (pT integrated)
    cout << endl;
    cout << setw(5) << "$\\sqrt{s_{_{NN}}}$" << " & "
         << setw(10) << "Correlation" << " & "
         << setw(10) << "Centrality" << " & "
         << setw(12) << "$p_T$ [GeV/c]" << " & "
         << setw(8) << "$C_1$" << " & "
         << setw(8) << "$C_2$" << " & "
         << setw(8) << "$C_3$" << " \\\\"
         << endl;
    cout << "\\hline" << endl;

    for (int icorr = 0; icorr < NCORR; icorr++)
    {
      for (int ic = 0; ic < NC; ic++)
      {
        cout << setw(5) << energy << " & "
             << setw(10) << cCorr[icorr] << " & "
             << setw(10) << Form("%i--%i\\%%", cl[ic], ch[ic]) << " & "
             << setw(12) << Form("[%.2f, %.2f]", ptl[ptsuml], pth[ptsumh]);
        for (int ipar = 0; ipar < 3; ipar++)
          cout << " & " << setw(8) << Form("%.5f$\\pm$%.5f", cn[icorr][ic][ipar].first, cn[icorr][ic][ipar].second);
        // cout << " & " << setw(8) << Form("%.5f$\\pm$%.5f (%.5f$\\pm$%.5f)", cn[icorr][ic][ipar].first, cn[icorr][ic][ipar].second, cn_num[icorr][ic][ipar].first, cn_num[icorr][ic][ipar].second);
        cout << " \\\\" << endl;

      } // ic
      cout << "\\hline" << endl;

    } // icorr


    //-- Raw v2 from 3-sub event methods, 0-5% vs pT
    cout << endl;
    cout << endl;
    cout << setw(5) << "$\\sqrt{s_{_{NN}}}$" << " & "
         << setw(10) << "Centrality" << " & "
         << setw(10) << "$p_T$ [GeV/c]" << " & "
         << setw(10) << "$v_{2}$(CNT--FVTXN--FVTXS)" << " & "
         << setw(10) << "$v_{2}$(CNT--FVTXN--BBCS)" << " & "
         << setw(10) << "$v_{2}$(CNT--BBCN--FVTXS)" << " & "
         << setw(10) << "$v_{2}$(CNT--BBCS--FVTXS)" << " \\\\"
         << endl;
    cout << "\\hline" << endl;
    for (int ic = 0; ic < NC; ic++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
      {
        cout << setw(5) << energy << " & "
             << setw(10) << Form("%i--%i\\%%", cl[ic], ch[ic]) << " & "
             << setw(10) << Form("[%.2f, %.2f]", ptl[ipt], pth[ipt]) << " & "
             << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTFVTXSFVTXN_pT[ic][1][ipt].first, vn_CNTFVTXSFVTXN_pT[ic][1][ipt].second) << " & "
             << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTBBCSFVTXN_pT[ic][1][ipt].first, vn_CNTBBCSFVTXN_pT[ic][1][ipt].second) << " & "
             << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTFVTXSBBCN_pT[ic][1][ipt].first, vn_CNTFVTXSBBCN_pT[ic][1][ipt].second) << " & "
             << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTBBCSFVTXS_pT[ic][1][ipt].first, vn_CNTBBCSFVTXS_pT[ic][1][ipt].second) << " \\\\ "
             << endl;
      }
      cout << "\\hline" << endl;
    }


    //-- Raw v2 from 3-sub event methods, pT integrated vs centrality
    cout << endl;
    cout << endl;
    cout << setw(5) << "$\\sqrt{s_{_{NN}}}$" << " & "
         << setw(10) << "Centrality" << " & "
         << setw(10) << "$p_T$ [GeV/c]" << " & "
         << setw(10) << "$v_{2}$(CNT--FVTXN--FVTXS)" << " & "
         << setw(10) << "$v_{2}$(CNT--FVTXN--BBCS)" << " & "
         << setw(10) << "$v_{2}$(CNT--BBCN--FVTXS)" << " & "
         << setw(10) << "$v_{2}$(CNT--BBCS--FVTXS)" << " \\\\"
         << endl;
    cout << "\\hline" << endl;
    for (int ic = 0; ic < NC; ic++)
    {
      cout << setw(5) << energy << " & "
           << setw(10) << Form("%i--%i\\%%", cl[ic], ch[ic]) << " & "
           << setw(10) << Form("[%.2f, %.2f]", ptl[ptsuml], pth[ptsumh]) << " & "
           << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTFVTXSFVTXN_cent[ic][1].first, vn_CNTFVTXSFVTXN_cent[ic][1].second) << " & "
           << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTBBCSFVTXN_cent[ic][1].first, vn_CNTBBCSFVTXN_cent[ic][1].second) << " & "
           << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTFVTXSBBCN_cent[ic][1].first, vn_CNTFVTXSBBCN_cent[ic][1].second) << " & "
           << setw(10) << Form("%.4f$\\pm$%.4f", vn_CNTBBCSFVTXS_cent[ic][1].first, vn_CNTBBCSFVTXS_cent[ic][1].second) << " \\\\ "
           << endl;
    }
    cout << "\\hline" << endl;

  } // if (printTables)


}
///-----------------------------------------------------------------------------



///-----------------------------------------------------------------------------
///
/// Calculate v_n using the 3 sub-event method
///
ValErr vn_3sub(const ValErr &C_AB, const ValErr &C_AC, const ValErr &C_BC)
{
  // return values;
  double vn = 0;
  double err = 0;

  // Calculate the value
  vn = TMath::Sqrt( (C_AB.first * C_AC.first) / C_BC.first );

  // Calculate the uncertainty using error propogation
  double dC_AB = +0.5 * (C_AC.first / C_BC.first) / vn;
  double dC_AC = +0.5 * (C_AB.first / C_BC.first) / vn;
  double dC_BC = -0.5 * (C_AB.first * C_AC.first) / (C_BC.first * C_BC.first) / vn;

  err += TMath::Power(dC_AB * C_AB.second, 2);
  err += TMath::Power(dC_AC * C_AC.second, 2);
  err += TMath::Power(dC_BC * C_BC.second, 2);
  err = TMath::Sqrt(err);

  // check for nan's (imaginary numbers)
  if (vn != vn)
  {
    vn = 0;
    err = 0;
  }

  return make_pair(vn, err);

}
///-----------------------------------------------------------------------------



///-----------------------------------------------------------------------------
///
/// Calculate C_n coefficients from the histogram numerically
///
ValErr calc_cn(TH1D* hcorr, int order)
{
  if (order < 0 || order > 5)
  {
    cout << "WARNING!! calc_cn() - order parameter " << order
         << " is out of range, returnning 0's" << endl;
    return make_pair(0, 0);
  }
  if (!hcorr)
  {
    cout << "ERROR!! calc_cn() - hcorr not a valid pointer! Returning 0's!" << endl;
    return make_pair(0, 0);
  }

  double val = 0;
  double err = 0;
  double tot = 0;
  for (int ibin = 1; ibin <= hcorr->GetNbinsX(); ibin++)
  {
    double dphi = hcorr->GetBinCenter(ibin);
    double bc = hcorr->GetBinContent(ibin);
    double e = hcorr->GetBinError(ibin);

    val += bc * TMath::Cos((float)order * dphi);
    err += TMath::Power(e * TMath::Cos((float)order * dphi), 2);
    tot += bc;

  }
  val = val / tot;
  err = TMath::Sqrt(err) / tot;

  return make_pair(val, err);
}
///-----------------------------------------------------------------------------

