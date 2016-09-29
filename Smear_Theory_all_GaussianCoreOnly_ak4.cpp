#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <math.h>
#include <TF1.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
 
#include <Riostream.h>
#include <Strlen.h>
#include <TDirectory.h>
#include <TClassTable.h>
#include <TInterpreter.h>
#include <THashList.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TError.h>
#include <TClass.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TVirtualMutex.h>
#include <TMath.h>

//#include <TThreadSlots.h>

#include <iostream>
#include <cmath>
#include "TSpline.h"
#include <TH1D.h>

#include "RooUnfoldResponse.h"


using namespace std;

void Column_Normalization(TH2D *h1,TH2D *h2);
void Raw_Normalization(TH2D *h1,TH2D *h2);

TH1D *cross_section(TH1D *Spectra_);

double my_crystal_ball(double x, double N, double mu, double sig, double a1, double p1, double a2, double p2){
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}



////////// NPs from Paolo 17/9/2015
//Dear all,
//please find the fits for NP corrections (ak4 jets). The function is:

/// y = [0] +[1]*pow(pT,[2]) with the three parameters below
///                      [0]            [1]             [2]
/// rapidity bin 1:   0.271307      0.771988       -0.0108214
/// rapidity bin 2:   0.27094       0.771557       -0.0106443
/// rapidity bin 3:   0.268122      0.768218       -0.00895077
/// rapidity bin 4:   0.263524      0.762893       -0.00660873
/// rapidity bin 5:   0.258164      0.756676       -0.00359792
/// rapidity bin 6:   0.277025      0.778438       -0.0140284
/// rapidity bin 7:   -65.5857      66.7221        -0.000413302


double NP_y0(double x){return (0.271307  + 0.771988  * pow (x ,-0.0108214));}
double NP_y1(double x){return (0.27094   + 0.771557  * pow (x ,-0.0106443));}
double NP_y2(double x){return (0.268122  + 0.768218  * pow (x ,-0.00895077));}
double NP_y3(double x){return (0.263524  + 0.762893  * pow (x ,-0.00660873));}
double NP_y4(double x){return (0.258164  + 0.756676  * pow (x ,-0.00359792));}
double NP_y5(double x){return (0.277025  + 0.778438  * pow (x ,-0.0140284));}
double NP_y6(double x){return (-65.5857  + 66.7221   * pow (x ,-0.000413302));}



////////// Gaussian Core Distributions from Panos 19/4/2016   25ns P8CUETM1
double gres_y0(double x){return ( 9.72078e-03 + 5.21606e-01/(pow(x , 3.81938e-01)+( 1.93633e-03*x)) + ( 2.84308e-06*x) );}
double gres_y1(double x){return ( 1.43784e-02 + 6.76487e-01/(pow(x , 4.60585e-01)+(-3.70305e-03*x)) + ( 5.15904e-07*x) );}
double gres_y2(double x){return ( 6.04647e-02 + 6.41042e-03/(pow(x ,-1.04867e+00)+( 1.05028e-03*x)) + (-5.98738e-06*x) );}
double gres_y3(double x){return ( 4.59239e-02 + 9.18303e-02/(pow(x ,-4.22337e-01)+( 1.36487e-02*x)) + (-8.35534e-06*x) );}
double gres_y4(double x){return ( 3.63903e-02 + 1.70734e+00/(pow(x , 8.71882e-01)+(-2.53203e-01*x)) + (-1.06428e-05*x) );}
double gres_y5(double x){return ( 4.77114e-02 + 8.85395e-01/(pow(x , 8.86030e-01)+(-4.45290e-01*x)) + (-7.69834e-05*x) );}
double gres_y6(double x){return (-7.59909e-02 + 5.81719e-01/(pow(x , 2.17956e-01)+( 1.73311e-03*x)) );}




void Smear_Theory_all_GaussianCoreOnly_ak4()
{

  TFile *outf    = new TFile("Smear_Theory_CT14_GaussianCoreOnly.root","RECREATE");

 // Theory spectrum
 TFile *f1 = new TFile("PredictionsFastNLO-ak4-CT14nlo.root");

 TH1D *theory_y0 = (TH1D*)f1->Get("InclusiveJet_1bin;1");
 TH1D *theory_y1 = (TH1D*)f1->Get("InclusiveJet_2bin;1");
 TH1D *theory_y2 = (TH1D*)f1->Get("InclusiveJet_3bin;1");
 TH1D *theory_y3 = (TH1D*)f1->Get("InclusiveJet_4bin;1");
 TH1D *theory_y4 = (TH1D*)f1->Get("InclusiveJet_5bin;1");
 TH1D *theory_y5 = (TH1D*)f1->Get("InclusiveJet_6bin;1");
 TH1D *theory_y6 = (TH1D*)f1->Get("InclusiveJet_7bin;1");


/////////////// Inclusive jets Binning ////////////////////////////////////////////////////////
     double newBins_incl[81]={0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
                              97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
                              507, 548, 592, 638, 686, 737, 790, 846, 905, 967,                    
                              1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                              2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
                              4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

     int n_newbins_incl=80;


     ////////////////////////////
     double xbins_incl[80]={1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
                            97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
                            507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
                            1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                            2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
	                    4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
     int nbins_incl=79;


     ////////////////////////////
     double reduced_xbins_incl[65]={56, 64, 74, 84,
                            97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
                            507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
                            1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                            2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
	                    4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
     int reduced_nbins_incl=64;






/////////////// NPs to theory spectra ////////////////////////////////////////////////////////
/////////////// Create functions to check NPs

  TF1 *fNP_y0 = new TF1("fNP_y0","NP_y0(x)",32,6000);
  TF1 *fNP_y1 = new TF1("fNP_y1","NP_y1(x)",32,5200);
  TF1 *fNP_y2 = new TF1("fNP_y2","NP_y2(x)",32,4000);
  TF1 *fNP_y3 = new TF1("fNP_y3","NP_y3(x)",32,2500);
  TF1 *fNP_y4 = new TF1("fNP_y4","NP_y4(x)",32,1588);
  TF1 *fNP_y5 = new TF1("fNP_y5","NP_y5(x)",32,967);
  TF1 *fNP_y6 = new TF1("fNP_y6","NP_y6(x)",32,468);

  /// check NPs
  TCanvas *plotNPs = new TCanvas("plotNPs", "plotNPs",1200,800);
  plotNPs->Divide(4,2);
  plotNPs->cd(1);
  fNP_y0->Draw();
  plotNPs->cd(2);
  fNP_y1->Draw();
  plotNPs->cd(3);
  fNP_y2->Draw();
  plotNPs->cd(4);
  fNP_y3->Draw();
  plotNPs->cd(5);
  fNP_y4->Draw();
  plotNPs->cd(6);
  fNP_y5->Draw();
  plotNPs->cd(7);
  fNP_y6->Draw();
  plotNPs->Update();


/////////////// Core p_T Resolution ////////////////////////////////////////////////////////
/////////////// Create functions

  TF1 *fres_y0 = new TF1("fres_y0","gres_y0(x)",56,3000);
  TF1 *fres_y1 = new TF1("fres_y1","gres_y1(x)",56,3000);
  TF1 *fres_y2 = new TF1("fres_y2","gres_y2(x)",56,2500);
  TF1 *fres_y3 = new TF1("fres_y3","gres_y3(x)",56,1600);
  TF1 *fres_y4 = new TF1("fres_y4","gres_y4(x)",56,1100);
  TF1 *fres_y5 = new TF1("fres_y5","gres_y5(x)",56,600);
  TF1 *fres_y6 = new TF1("fres_y6","gres_y6(x)",56,180);

  fres_y0->SetLineColor(2);
  fres_y0->SetMinimum(0);fres_y0->SetMaximum(0.25);
  fres_y1->SetLineColor(3);
  fres_y2->SetLineColor(4);
  fres_y3->SetLineColor(6);
  fres_y4->SetLineColor(kOrange+1);
  fres_y5->SetLineColor(kCyan-3);
  fres_y6->SetLineColor(kYellow-5);




  fres_y0->SetTitle("Gaussian Core Resolution (#sigma);jet p_{T}     ; #sigma         ");

  /// check Core p_T Resolution
  TCanvas *plotres = new TCanvas("plotres", "plotres",1200,800);
  plotres->cd();
  fres_y0->Draw();
  fres_y1->Draw("same");
  fres_y2->Draw("same");
  fres_y3->Draw("same");
  fres_y4->Draw("same");
  fres_y5->Draw("same");
  fres_y6->Draw("same");

  TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
  leg->AddEntry(fres_y0,"PYTHIA8 CUETM1 @ 13TeV", "");
  leg->AddEntry(fres_y0,"AK4 and PFchs Jets", "");
  leg->AddEntry(fres_y0,"0<|y|<0.5", "l"); 
  leg->AddEntry(fres_y1,"0.5<|y|<1.0", "l"); 
  leg->AddEntry(fres_y2,"1.0<|y|<1.5", "l"); 
  leg->AddEntry(fres_y3,"1.5<|y|<2.0", "l"); 
  leg->AddEntry(fres_y4,"2.0<|y|<2.5", "l"); 
  leg->AddEntry(fres_y5,"2.5<|y|<3.0", "l"); 
  leg->AddEntry(fres_y6,"3.2<|y|<4.7", "l"); 
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->Draw();

  plotres->Update();




  //////////// Apply NPs to theory
  /// Cross section y0
  TH1D *Xsection_theory_y0=(TH1D*)theory_y0->Clone();
  for(int i=0; i<theory_y0->GetNbinsX(); ++i){
     double bin_value  = theory_y0->GetBinContent(1+i);
     double bin_center = theory_y0->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y0(bin_center);
     Xsection_theory_y0->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y1
  TH1D *Xsection_theory_y1=(TH1D*)theory_y1->Clone();
  for(int i=0; i<theory_y1->GetNbinsX(); ++i){
     double bin_value  = theory_y1->GetBinContent(1+i);
     double bin_center = theory_y1->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y1(bin_center);
     Xsection_theory_y1->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y2
  TH1D *Xsection_theory_y2=(TH1D*)theory_y2->Clone();
  for(int i=0; i<theory_y2->GetNbinsX(); ++i){
     double bin_value  = theory_y2->GetBinContent(1+i);
     double bin_center = theory_y2->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y2(bin_center);
     Xsection_theory_y2->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y3
  TH1D *Xsection_theory_y3=(TH1D*)theory_y3->Clone();
  for(int i=0; i<theory_y3->GetNbinsX(); ++i){
     double bin_value  = theory_y3->GetBinContent(1+i);
     double bin_center = theory_y3->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y3(bin_center);
     Xsection_theory_y3->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y4
  TH1D *Xsection_theory_y4=(TH1D*)theory_y4->Clone();
  for(int i=0; i<theory_y4->GetNbinsX(); ++i){
     double bin_value  = theory_y4->GetBinContent(1+i);
     double bin_center = theory_y4->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y4(bin_center);
     Xsection_theory_y4->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y5
  TH1D *Xsection_theory_y5=(TH1D*)theory_y5->Clone();
  for(int i=0; i<theory_y5->GetNbinsX(); ++i){
     double bin_value  = theory_y5->GetBinContent(1+i);
     double bin_center = theory_y5->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y5(bin_center);
     Xsection_theory_y5->SetBinContent(1+i,new_bin_value);
  }

  /// Cross section y6
  TH1D *Xsection_theory_y6=(TH1D*)theory_y6->Clone();
  for(int i=0; i<theory_y6->GetNbinsX(); ++i){
     double bin_value  = theory_y6->GetBinContent(1+i);
     double bin_center = theory_y6->GetBinCenter(1+i);
     double new_bin_value = bin_value * NP_y6(bin_center);
     Xsection_theory_y6->SetBinContent(1+i,new_bin_value);
  }



//  TH1D *ratio_y6=(TH1D*)theory_y6->Clone(); ratio_y6->SetTitle("ratio_y6; p_{T} (GeV); xsection/theory y6   ");
//  ratio_y6->Divide(Xsection_theory_y6,theory_y6,1.,1.,"B");
//  TCanvas *cc = new TCanvas("cc", "cc",800,600);
//  cc->cd();gPad->SetLogy();
//  ratio_y6->Draw();
//  fNP_y6->Draw("same");
//  cc->Update();




///Create Cubic Splines using Cross sections
     TSpline3 *spline3_y0 = new TSpline3(Xsection_theory_y0);
     TSpline3 *spline3_y1 = new TSpline3(Xsection_theory_y1);
     TSpline3 *spline3_y2 = new TSpline3(Xsection_theory_y2);
     TSpline3 *spline3_y3 = new TSpline3(Xsection_theory_y3);
     TSpline3 *spline3_y4 = new TSpline3(Xsection_theory_y4);
     TSpline3 *spline3_y5 = new TSpline3(Xsection_theory_y5);
     TSpline3 *spline3_y6 = new TSpline3(Xsection_theory_y6);


     spline3_y0->SetLineColor(6);
     spline3_y1->SetLineColor(6);
     spline3_y2->SetLineColor(6);
     spline3_y3->SetLineColor(6);
     spline3_y4->SetLineColor(6);
     spline3_y5->SetLineColor(6);
     spline3_y6->SetLineColor(6);


  //Xsection_theory_y0->SetAxisRange(56,5777,"X"); 
  //Xsection_theory_y1->SetAxisRange(56,5220,"X"); 
  //Xsection_theory_y2->SetAxisRange(56,3832,"X"); 
  //Xsection_theory_y3->SetAxisRange(56,2500,"X"); 
  //Xsection_theory_y4->SetAxisRange(56,1588,"X"); 
  //Xsection_theory_y5->SetAxisRange(56,967,"X"); 
  //Xsection_theory_y6->SetAxisRange(56,448,"X"); 


  Xsection_theory_y0->SetAxisRange(56,5700,"X"); 
  Xsection_theory_y1->SetAxisRange(56,5000,"X"); 
  Xsection_theory_y2->SetAxisRange(56,3700,"X"); 
  Xsection_theory_y3->SetAxisRange(56,2400,"X"); 
  Xsection_theory_y4->SetAxisRange(56,1588,"X"); 
  Xsection_theory_y5->SetAxisRange(56,900,"X"); 
  Xsection_theory_y6->SetAxisRange(56,448,"X"); 



  /// check spectra
  TCanvas *TheorySpectra = new TCanvas("TheorySpectra", "Theory Spectra",1200,800);
  TheorySpectra->Divide(4,2);
  TheorySpectra->cd(1);gPad->SetLogy();
  Xsection_theory_y0->Draw();
  spline3_y0->Draw("same");
  TheorySpectra->cd(2);gPad->SetLogy();
  Xsection_theory_y1->Draw();
  spline3_y1->Draw("same");
  TheorySpectra->cd(3);gPad->SetLogy();
  Xsection_theory_y2->Draw();
  spline3_y2->Draw("same");  
  TheorySpectra->cd(4);gPad->SetLogy();
  Xsection_theory_y3->Draw();
  spline3_y3->Draw("same");  
  TheorySpectra->cd(5);gPad->SetLogy();
  Xsection_theory_y4->Draw();
  spline3_y4->Draw("same");  
  TheorySpectra->cd(6);gPad->SetLogy();
  Xsection_theory_y5->Draw();
  spline3_y5->Draw("same");  
  TheorySpectra->cd(7);gPad->SetLogy();
  Xsection_theory_y6->Draw();
  spline3_y6->Draw("same");  


  TCanvas *TheorySplineSpectra = new TCanvas("TheorySplineSpectra", "Theory Spline Spectra",1200,600);
  TheorySplineSpectra->Divide(4,2);
  TheorySplineSpectra->cd(1);gPad->SetLogy();
  spline3_y0->Draw();
  TheorySplineSpectra->cd(2);gPad->SetLogy();
  spline3_y1->Draw();
  TheorySplineSpectra->cd(3);gPad->SetLogy();
  spline3_y2->Draw();
  TheorySplineSpectra->cd(4);gPad->SetLogy();
  spline3_y3->Draw();
  TheorySplineSpectra->cd(5);gPad->SetLogy();
  spline3_y4->Draw();
  TheorySplineSpectra->cd(6);gPad->SetLogy();
  spline3_y5->Draw();
  TheorySplineSpectra->cd(7);gPad->SetLogy();
  spline3_y6->Draw();






  TRandom3 *rnd = new TRandom3();
  int nEvents=100000000;

  TH1D *True_incl_jet_y0    = new TH1D("True_incl_jet_y0","True_incl_jet_y0", n_newbins_incl, newBins_incl); True_incl_jet_y0->Sumw2();
  TH1D *Smeared_incl_jet_y0 = new TH1D("Smeared_incl_jet_y0","Smeared_incl_jet_y0", n_newbins_incl, newBins_incl); Smeared_incl_jet_y0->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y0    = new TH1D("NormBinWidth_True_incl_jet_y0","NormBinWidth_True_incl_jet_y0", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y0->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y0 = new TH1D("NormBinWidth_Smeared_incl_jet_y0","NormBinWidth_Smeared_incl_jet_y0", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y0->Sumw2();
  RooUnfoldResponse response_y0(Smeared_incl_jet_y0,True_incl_jet_y0);

  TH1D *True_incl_jet_y1    = new TH1D("True_incl_jet_y1","True_incl_jet_y1", n_newbins_incl, newBins_incl); True_incl_jet_y1->Sumw2();
  TH1D *Smeared_incl_jet_y1 = new TH1D("Smeared_incl_jet_y1","Smeared_incl_jet_y1", n_newbins_incl, newBins_incl); Smeared_incl_jet_y1->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y1    = new TH1D("NormBinWidth_True_incl_jet_y1","NormBinWidth_True_incl_jet_y1", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y1->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y1 = new TH1D("NormBinWidth_Smeared_incl_jet_y1","NormBinWidth_Smeared_incl_jet_y1", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y1->Sumw2();
  RooUnfoldResponse response_y1(Smeared_incl_jet_y1,True_incl_jet_y1);

  TH1D *True_incl_jet_y2    = new TH1D("True_incl_jet_y2","True_incl_jet_y2", n_newbins_incl, newBins_incl); True_incl_jet_y2->Sumw2();
  TH1D *Smeared_incl_jet_y2 = new TH1D("Smeared_incl_jet_y2","Smeared_incl_jet_y2", n_newbins_incl, newBins_incl); Smeared_incl_jet_y2->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y2    = new TH1D("NormBinWidth_True_incl_jet_y2","NormBinWidth_True_incl_jet_y2", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y2->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y2 = new TH1D("NormBinWidth_Smeared_incl_jet_y2","NormBinWidth_Smeared_incl_jet_y2", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y2->Sumw2();
  RooUnfoldResponse response_y2(Smeared_incl_jet_y2,True_incl_jet_y2);

  TH1D *True_incl_jet_y3    = new TH1D("True_incl_jet_y3","True_incl_jet_y3", n_newbins_incl, newBins_incl); True_incl_jet_y3->Sumw2();
  TH1D *Smeared_incl_jet_y3 = new TH1D("Smeared_incl_jet_y3","Smeared_incl_jet_y3", n_newbins_incl, newBins_incl); Smeared_incl_jet_y3->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y3    = new TH1D("NormBinWidth_True_incl_jet_y3","NormBinWidth_True_incl_jet_y3", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y3->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y3 = new TH1D("NormBinWidth_Smeared_incl_jet_y3","NormBinWidth_Smeared_incl_jet_y3", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y3->Sumw2();
  RooUnfoldResponse response_y3(Smeared_incl_jet_y3,True_incl_jet_y3);

  TH1D *True_incl_jet_y4    = new TH1D("True_incl_jet_y4","True_incl_jet_y4", n_newbins_incl, newBins_incl); True_incl_jet_y4->Sumw2();
  TH1D *Smeared_incl_jet_y4 = new TH1D("Smeared_incl_jet_y4","Smeared_incl_jet_y4", n_newbins_incl, newBins_incl); Smeared_incl_jet_y4->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y4    = new TH1D("NormBinWidth_True_incl_jet_y4","NormBinWidth_True_incl_jet_y4", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y4->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y4 = new TH1D("NormBinWidth_Smeared_incl_jet_y4","NormBinWidth_Smeared_incl_jet_y4", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y4->Sumw2();
  RooUnfoldResponse response_y4(Smeared_incl_jet_y4,True_incl_jet_y4);

  TH1D *True_incl_jet_y5    = new TH1D("True_incl_jet_y5","True_incl_jet_y5", n_newbins_incl, newBins_incl); True_incl_jet_y5->Sumw2();
  TH1D *Smeared_incl_jet_y5 = new TH1D("Smeared_incl_jet_y5","Smeared_incl_jet_y5", n_newbins_incl, newBins_incl); Smeared_incl_jet_y5->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y5    = new TH1D("NormBinWidth_True_incl_jet_y5","NormBinWidth_True_incl_jet_y5", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y5->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y5 = new TH1D("NormBinWidth_Smeared_incl_jet_y5","NormBinWidth_Smeared_incl_jet_y5", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y5->Sumw2();
  RooUnfoldResponse response_y5(Smeared_incl_jet_y5,True_incl_jet_y5);

  TH1D *True_incl_jet_y6    = new TH1D("True_incl_jet_y6","True_incl_jet_y6", n_newbins_incl, newBins_incl); True_incl_jet_y6->Sumw2();
  TH1D *Smeared_incl_jet_y6 = new TH1D("Smeared_incl_jet_y6","Smeared_incl_jet_y6", n_newbins_incl, newBins_incl); Smeared_incl_jet_y6->Sumw2();
  TH1D *NormBinWidth_True_incl_jet_y6    = new TH1D("NormBinWidth_True_incl_jet_y6","NormBinWidth_True_incl_jet_y6", n_newbins_incl, newBins_incl); NormBinWidth_True_incl_jet_y6->Sumw2();
  TH1D *NormBinWidth_Smeared_incl_jet_y6 = new TH1D("NormBinWidth_Smeared_incl_jet_y6","NormBinWidth_Smeared_incl_jet_y6", n_newbins_incl, newBins_incl); NormBinWidth_Smeared_incl_jet_y6->Sumw2();
  RooUnfoldResponse response_y6(Smeared_incl_jet_y6,True_incl_jet_y6);



  printf("Evaluating y0\n");

  //// y0
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 5777;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y0(ptTrue);

     sigma=sigma*1.086;   /// October 2015

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y0->Eval(ptTrue);

     True_incl_jet_y0    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y0 -> Fill(ptSmeared,pt_w);
     
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);

     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y0.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y0.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y0.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y0    = cross_section(True_incl_jet_y0);
   NormBinWidth_Smeared_incl_jet_y0 = cross_section(Smeared_incl_jet_y0);
   NormBinWidth_True_incl_jet_y0->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y0->SetLineColor(2);


  
  printf("Evaluating y1\n");

  //// y1
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 5220;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y1(ptTrue);

     sigma=sigma*1.128;   /// October 2015

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y1->Eval(ptTrue);

     True_incl_jet_y1    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y1 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
     
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y1.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y1.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y1.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y1    = cross_section(True_incl_jet_y1);
   NormBinWidth_Smeared_incl_jet_y1 = cross_section(Smeared_incl_jet_y1);
   NormBinWidth_True_incl_jet_y1->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y1->SetLineColor(2);


  printf("Evaluating y2\n");

  //// y2
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 3832;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y2(ptTrue);

     sigma=sigma*1.143;    /// October 2015

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y2->Eval(ptTrue);

     True_incl_jet_y2    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y2 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
     
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y2.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y2.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y2.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y2    = cross_section(True_incl_jet_y2);
   NormBinWidth_Smeared_incl_jet_y2 = cross_section(Smeared_incl_jet_y2);
   NormBinWidth_True_incl_jet_y2->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y2->SetLineColor(2);



  printf("Evaluating y3\n");

  //// y3
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 2500;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y3(ptTrue);

     sigma=sigma*1.109;   /// October 2015

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y3->Eval(ptTrue);

     True_incl_jet_y3    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y3 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
 
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y3.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y3.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y3.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y3    = cross_section(True_incl_jet_y3);
   NormBinWidth_Smeared_incl_jet_y3 = cross_section(Smeared_incl_jet_y3);
   NormBinWidth_True_incl_jet_y3->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y3->SetLineColor(2);


  printf("Evaluating y4\n");


  //// y4
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 1588;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y4(ptTrue);     
          
     sigma=sigma*1.109;   /// October 2015

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y4->Eval(ptTrue);

     True_incl_jet_y4    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y4 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
 
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y4.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y4.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y4.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y4    = cross_section(True_incl_jet_y4);
   NormBinWidth_Smeared_incl_jet_y4 = cross_section(Smeared_incl_jet_y4);
   NormBinWidth_True_incl_jet_y4->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y4->SetLineColor(2);


  printf("Evaluating y5\n");

  //// y5
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 967;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y5(ptTrue);     
  
     sigma=sigma*1.109;     ///1.254 From 8 TeV actually

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y5->Eval(ptTrue);

     True_incl_jet_y5    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y5 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);

     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y5.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y5.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y5.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y5    = cross_section(True_incl_jet_y5);
   NormBinWidth_Smeared_incl_jet_y5 = cross_section(Smeared_incl_jet_y5);
   NormBinWidth_True_incl_jet_y5->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y5->SetLineColor(2);


  printf("Evaluating y6\n");

  //// y6
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 468;   //468;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y6(ptTrue);     

     sigma=sigma*1.056;    /// From 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y6->Eval(ptTrue);

     True_incl_jet_y6    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y6 -> Fill(ptSmeared,pt_w);

     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
     
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y6.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y6.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y6.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y6    = cross_section(True_incl_jet_y6);
   NormBinWidth_Smeared_incl_jet_y6 = cross_section(Smeared_incl_jet_y6);
   NormBinWidth_True_incl_jet_y6->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y6->SetLineColor(2);



   TH1F *rat_y0_true_smeared=(TH1F*)True_incl_jet_y0->Clone(); rat_y0_true_smeared->SetTitle("y0 True/Smeared;p_T (GeV); y0 True/Smeared          ");
   rat_y0_true_smeared->Divide(True_incl_jet_y0,Smeared_incl_jet_y0,1.,1.,"B");

   TH1F *rat_y1_true_smeared=(TH1F*)True_incl_jet_y1->Clone(); rat_y1_true_smeared->SetTitle("y1 True/Smeared;p_T (GeV); y1 True/Smeared          ");
   rat_y1_true_smeared->Divide(True_incl_jet_y1,Smeared_incl_jet_y1,1.,1.,"B");

   TH1F *rat_y2_true_smeared=(TH1F*)True_incl_jet_y2->Clone(); rat_y2_true_smeared->SetTitle("y2 True/Smeared;p_T (GeV); y2 True/Smeared          ");
   rat_y2_true_smeared->Divide(True_incl_jet_y2,Smeared_incl_jet_y2,1.,1.,"B");

   TH1F *rat_y3_true_smeared=(TH1F*)True_incl_jet_y3->Clone(); rat_y3_true_smeared->SetTitle("y3 True/Smeared;p_T (GeV); y3 True/Smeared          ");
   rat_y3_true_smeared->Divide(True_incl_jet_y3,Smeared_incl_jet_y3,1.,1.,"B");

   TH1F *rat_y4_true_smeared=(TH1F*)True_incl_jet_y4->Clone(); rat_y4_true_smeared->SetTitle("y4 True/Smeared;p_T (GeV); y4 True/Smeared          ");
   rat_y4_true_smeared->Divide(True_incl_jet_y4,Smeared_incl_jet_y4,1.,1.,"B");

   TH1F *rat_y5_true_smeared=(TH1F*)True_incl_jet_y5->Clone(); rat_y5_true_smeared->SetTitle("y5 True/Smeared;p_T (GeV); y5 True/Smeared          ");
   rat_y5_true_smeared->Divide(True_incl_jet_y5,Smeared_incl_jet_y5,1.,1.,"B");

   TH1F *rat_y6_true_smeared=(TH1F*)True_incl_jet_y6->Clone(); rat_y6_true_smeared->SetTitle("y6 True/Smeared;p_T (GeV); y6 True/Smeared          ");
   rat_y6_true_smeared->Divide(True_incl_jet_y6,Smeared_incl_jet_y6,1.,1.,"B");

   
   
   rat_y0_true_smeared->SetMaximum(1.1); rat_y0_true_smeared->SetMinimum(0); rat_y0_true_smeared->SetAxisRange(97,6000,"X"); 
   rat_y1_true_smeared->SetMaximum(1.1); rat_y1_true_smeared->SetMinimum(0);  rat_y1_true_smeared->SetAxisRange(97,5200,"X"); 
   rat_y2_true_smeared->SetMaximum(1.1); rat_y2_true_smeared->SetMinimum(0);  rat_y2_true_smeared->SetAxisRange(97,4000,"X"); 
   rat_y3_true_smeared->SetMaximum(1.1); rat_y3_true_smeared->SetMinimum(0);  rat_y3_true_smeared->SetAxisRange(97,2500,"X"); 
   rat_y4_true_smeared->SetMaximum(1.1); rat_y4_true_smeared->SetMinimum(0);  rat_y4_true_smeared->SetAxisRange(97,1588,"X"); 
   rat_y5_true_smeared->SetMaximum(1.1); rat_y5_true_smeared->SetMinimum(0);  rat_y5_true_smeared->SetAxisRange(97,967,"X"); 
   rat_y6_true_smeared->SetMaximum(1.1); rat_y6_true_smeared->SetMinimum(0);  rat_y6_true_smeared->SetAxisRange(97,448,"X"); 
   
   
   
   TCanvas *trueOverSmeared = new TCanvas("trueOverSmeared", "trueOverSmeared",1600,800);
   trueOverSmeared->Divide(4,2);
   trueOverSmeared->cd(1);gPad->SetGrid(); 
   rat_y0_true_smeared->Draw();     
   trueOverSmeared->cd(2);gPad->SetGrid(); 
   rat_y1_true_smeared->Draw();     
   trueOverSmeared->cd(3);gPad->SetGrid(); 
   rat_y2_true_smeared->Draw();     
   trueOverSmeared->cd(4);gPad->SetGrid(); 
   rat_y3_true_smeared->Draw();     
   trueOverSmeared->cd(5);gPad->SetGrid(); 
   rat_y4_true_smeared->Draw();     
   trueOverSmeared->cd(6); gPad->SetGrid();
   rat_y5_true_smeared->Draw();     
   trueOverSmeared->cd(7);gPad->SetGrid(); 
   rat_y6_true_smeared->Draw();     
   trueOverSmeared->Update();











   /// some checks Only in the begining Then swich off
  int do_checks=0;
  if(do_checks==1){
  
      double integ_Xsection_theory_y0=0,integ_True_incl_jet_y0_normBinWidth=0;
      for(int i=0;i<Xsection_theory_y0->GetNbinsX();++i){integ_Xsection_theory_y0=integ_Xsection_theory_y0+Xsection_theory_y0->GetBinContent(i+1);}
      for(int i=15;i<NormBinWidth_True_incl_jet_y0->GetNbinsX();++i){integ_True_incl_jet_y0_normBinWidth=integ_True_incl_jet_y0_normBinWidth+NormBinWidth_True_incl_jet_y0->GetBinContent(i+1);}
      NormBinWidth_True_incl_jet_y0->Scale(integ_Xsection_theory_y0/integ_True_incl_jet_y0_normBinWidth);


      TH1D *rat_y0=(TH1D*)True_incl_jet_y0->Clone();
      for(int i=0;i<Xsection_theory_y0->GetNbinsX();++i){ 
         if(True_incl_jet_y0->GetBinContent(i+16)!=0){
            rat_y0->SetBinContent(i+16,Xsection_theory_y0->GetBinContent(i+1)/NormBinWidth_True_incl_jet_y0->GetBinContent(i+16));
	    rat_y0->SetBinError(i+16,0);
         }
      }
      rat_y0->SetMaximum(1.2); rat_y0->SetMinimum(0.8);  rat_y0->GetXaxis()->SetRange(16,79);


     TCanvas *SmearedSpectra = new TCanvas("SmearedSpectra", "Smeared Spectra",1200,600);
     SmearedSpectra->Divide(4,1);
  
     SmearedSpectra->cd(1); gPad->SetLogy();
     Xsection_theory_y0->Draw();     
     NormBinWidth_True_incl_jet_y0->Draw("same"); 
     SmearedSpectra->cd(2);gPad->SetLogy();
     NormBinWidth_True_incl_jet_y0->Draw();     
     SmearedSpectra->cd(3);gPad->SetLogy();
     NormBinWidth_Smeared_incl_jet_y0->Draw();     
     SmearedSpectra->cd(4);
     rat_y0->Draw();


  }////End do checks


  ////// Write out histos
  outf->cd();
  outf->Write();
  Xsection_theory_y0->Write();
  True_incl_jet_y0->Write();
  NormBinWidth_True_incl_jet_y0->Write("NormBinWidth_True_incl_jet_y0");    
  Smeared_incl_jet_y0->Write(); 
  NormBinWidth_Smeared_incl_jet_y0->Write("NormBinWidth_Smeared_incl_jet_y0"); 
  response_y0.Write("response_y0");

  Xsection_theory_y1->Write();
  True_incl_jet_y1->Write();
  NormBinWidth_True_incl_jet_y1->Write("NormBinWidth_True_incl_jet_y1");    
  Smeared_incl_jet_y1->Write(); 
  NormBinWidth_Smeared_incl_jet_y1->Write("NormBinWidth_Smeared_incl_jet_y1"); 
  response_y1.Write("response_y1");

  Xsection_theory_y2->Write();
  True_incl_jet_y2->Write();
  NormBinWidth_True_incl_jet_y2->Write("NormBinWidth_True_incl_jet_y2");    
  Smeared_incl_jet_y2->Write(); 
  NormBinWidth_Smeared_incl_jet_y2->Write("NormBinWidth_Smeared_incl_jet_y2"); 
  response_y2.Write("response_y2");

  Xsection_theory_y3->Write();
  True_incl_jet_y3->Write();
  NormBinWidth_True_incl_jet_y3->Write("NormBinWidth_True_incl_jet_y3");    
  Smeared_incl_jet_y3->Write(); 
  NormBinWidth_Smeared_incl_jet_y3->Write("NormBinWidth_Smeared_incl_jet_y3"); 
  response_y3.Write("response_y3");

  Xsection_theory_y4->Write();
  True_incl_jet_y4->Write();
  NormBinWidth_True_incl_jet_y4->Write("NormBinWidth_True_incl_jet_y4");    
  Smeared_incl_jet_y4->Write(); 
  NormBinWidth_Smeared_incl_jet_y4->Write("NormBinWidth_Smeared_incl_jet_y4"); 
  response_y4.Write("response_y4");

  Xsection_theory_y5->Write();
  True_incl_jet_y5->Write();
  NormBinWidth_True_incl_jet_y5->Write("NormBinWidth_True_incl_jet_y5");    
  Smeared_incl_jet_y5->Write(); 
  NormBinWidth_Smeared_incl_jet_y5->Write("NormBinWidth_Smeared_incl_jet_y5"); 
  response_y5.Write("response_y5");

  Xsection_theory_y6->Write();
  True_incl_jet_y6->Write();
  NormBinWidth_True_incl_jet_y6->Write("NormBinWidth_True_incl_jet_y6");    
  Smeared_incl_jet_y6->Write(); 
  NormBinWidth_Smeared_incl_jet_y6->Write("NormBinWidth_Smeared_incl_jet_y6"); 
  response_y6.Write("response_y6");






return;




}



////////////////////////////////////////////////////////////////////
void Column_Normalization(TH2D *h1,TH2D *h2)
{
   for(int i=0;i<h1->GetNbinsX();i++) {
      double a=0;
      for(int j=0;j<h1->GetNbinsY();j++) {
         a=a+h1->GetBinContent(i+1,j+1);
      }  
      for(int j=0;j<h1->GetNbinsY();j++) {
         double x = h1->GetXaxis()->GetBinCenter(i+1);
	 double y = h1->GetYaxis()->GetBinCenter(j+1); 
	 h2->SetBinContent(i+1,j+1,0.);
	 if(a>0)h2->Fill(x,y,h1->GetBinContent(i+1,j+1)/a);
      }
   }
}

////////////////////////////////////////////////////////////////////
void Raw_Normalization(TH2D *h1,TH2D *h2)
{
   for(int j=0;j<h1->GetNbinsY();j++) {
      double a=0;
      for(int i=0;i<h1->GetNbinsX();i++) {
         a=a+h1->GetBinContent(i+1,j+1);
      }  
      for(int i=0;i<h1->GetNbinsX();i++) {
         double x = h1->GetXaxis()->GetBinCenter(i+1);
	 double y = h1->GetYaxis()->GetBinCenter(j+1);
	 h2->SetBinContent(i+1,j+1,0.);
	 if(a>0)h2->Fill(x,y,h1->GetBinContent(i+1,j+1)/a);
      }
   }
}



TH1D *cross_section(TH1D *Spectra_){

  /// Evaluate xSections
  TH1D *Xsection_=(TH1D*)Spectra_->Clone();
  for(int i=0; i<Spectra_->GetNbinsX(); ++i){
     float bin_value = Spectra_->GetBinContent(1+i);
     float bin_error = Spectra_->GetBinError(1+i);
     float bin_width = Spectra_->GetBinWidth(1+i);
     float new_bin_value = bin_value/bin_width;
     float new_bin_error = bin_error/bin_width;
     Xsection_->SetBinContent(1+i,new_bin_value);
     Xsection_->SetBinError(1+i,new_bin_error);
  }
  //   std::cout<<Integral_<<std::endl;

  return Xsection_;
}



