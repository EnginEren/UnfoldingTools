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
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"

using namespace std;

void unfoldAK4(){

    // Get measured spectra
    //TFile *h = new TFile("../measuredSData2016B-ak4.root");
    TFile *h = new TFile("measuredS-ak4Rebinned.root");
   


    TH1F *y1 = (TH1F*)gDirectory->Get("pt_DETInclJet_1bin");
    TH1F *y2 = (TH1F*)gDirectory->Get("pt_DETInclJet_2bin");
    TH1F *y3 = (TH1F*)gDirectory->Get("pt_DETInclJet_3bin");
    TH1F *y4 = (TH1F*)gDirectory->Get("pt_DETInclJet_4bin");
    TH1F *y5 = (TH1F*)gDirectory->Get("pt_DETInclJet_5bin");
    TH1F *y6 = (TH1F*)gDirectory->Get("pt_DETInclJet_6bin");
    TH1F *y7 = (TH1F*)gDirectory->Get("pt_DETInclJet_7bin");
   
    
    // Get response matrix
    TFile *r = new TFile("Smear_Theory_CT14_GaussianCoreOnly_AK4_80x.root");
    RooUnfoldResponse *responseY1 = (RooUnfoldResponse*)r->Get("response_y0");
    RooUnfoldResponse *responseY2 = (RooUnfoldResponse*)r->Get("response_y1");
    RooUnfoldResponse *responseY3 = (RooUnfoldResponse*)r->Get("response_y2");
    RooUnfoldResponse *responseY4 = (RooUnfoldResponse*)r->Get("response_y3");
    RooUnfoldResponse *responseY5 = (RooUnfoldResponse*)r->Get("response_y4");
    RooUnfoldResponse *responseY6 = (RooUnfoldResponse*)r->Get("response_y5");
    RooUnfoldResponse *responseY7 = (RooUnfoldResponse*)r->Get("response_y678");
    //std::cout << ((TH1*) responseY7->Htruth())->Sizeof()<< std::endl;
   
    int iteration = 5;

    RooUnfoldBayes unfoldY1 (responseY1, y1, iteration);
    RooUnfoldBayes unfoldY2 (responseY2, y2, iteration);
    RooUnfoldBayes unfoldY3 (responseY3, y3, iteration);
    RooUnfoldBayes unfoldY4 (responseY4, y4, iteration);
    RooUnfoldBayes unfoldY5 (responseY5, y5, iteration);
    RooUnfoldBayes unfoldY6 (responseY6, y6, iteration);
    RooUnfoldBayes unfoldY7 (responseY7, y7, iteration);


    RooUnfold::ErrorTreatment f = RooUnfold::kCovariance;
    TH1F *recy1= (TH1F*) unfoldY1.Hreco(f);
    TMatrixD erry1 = (TMatrixD) unfoldY1.Ereco(f);
    recy1->SetName("y1");
    
    TH1F *recy2= (TH1F*) unfoldY2.Hreco(f);
    TMatrixD erry2 = (TMatrixD) unfoldY2.Ereco(f);
    recy2->SetName("y2");

    TH1F *recy3= (TH1F*) unfoldY3.Hreco(f);
    TMatrixD erry3 = (TMatrixD) unfoldY3.Ereco(f);
    recy3->SetName("y3");
    
    TH1F *recy4= (TH1F*) unfoldY4.Hreco(f);
    TMatrixD erry4 = (TMatrixD) unfoldY4.Ereco(f);
    recy4->SetName("y4");
    
    TH1F *recy5= (TH1F*) unfoldY5.Hreco(f);
    TMatrixD erry5 = (TMatrixD) unfoldY5.Ereco(f);
    recy5->SetName("y5");
    
    TH1F *recy6= (TH1F*) unfoldY6.Hreco(f);
    TMatrixD erry6 = (TMatrixD) unfoldY6.Ereco(f);
    recy6->SetName("y6");
    
    TH1F *recy7= (TH1F*) unfoldY7.Hreco(f);
    TMatrixD erry7 = (TMatrixD) unfoldY7.Ereco(f);
    recy7->SetName("y7");
    //unfoldY1.Print();
  

    TFile *out = new TFile("Unfolded-TOYMC-2016C-ak4.root","RECREATE");
    recy1->Write();
    erry1.Write();
    recy2->Write();
    erry2.Write();
    recy3->Write();
    erry3.Write();
    recy4->Write();
    erry4.Write();
    recy5->Write();
    erry5.Write();
    recy6->Write();
    erry6.Write();
    recy7->Write();
    erry7.Write();

    out->Close();
    
    /*
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->SetLogy();
    c1->SetLogx();
    recy1->Draw();
    y1->SetLineColor(kRed);
    y1->Draw("same");
    */
}








