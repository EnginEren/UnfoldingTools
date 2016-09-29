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
#include <string>
#include <sstream>

//#include <TThreadSlots.h>

#include <iostream>
#include <cmath>
#include "TSpline.h"
#include <TH1D.h>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"

using namespace std;

void unfoldAK7_P8toDATA (){

    // Get measured spectra
    
    TFile *m = new TFile("../measuredSData2016BC-ak7.root"); 
    //TFile *m = new TFile("../measuredS-ak7DAS.root"); 
   
    TH1F *y1 = (TH1F*)gDirectory->Get("pt_DETInclJet_1bin;1");
    TH1F *y2 = (TH1F*)gDirectory->Get("pt_DETInclJet_2bin;1");
    TH1F *y3 = (TH1F*)gDirectory->Get("pt_DETInclJet_3bin;1");
    TH1F *y4 = (TH1F*)gDirectory->Get("pt_DETInclJet_4bin;1");
    TH1F *y5 = (TH1F*)gDirectory->Get("pt_DETInclJet_5bin;1");
    TH1F *y6 = (TH1F*)gDirectory->Get("pt_DETInclJet_6bin;1");
    TH1F *y7 = (TH1F*)gDirectory->Get("pt_DETInclJet_7bin;1");
    
    
    //TFile *h = new TFile("../MC.ak4.cmsJobs.P8Nuples.13:11:17-06-17-2016/P8-TuneCUETP8M_13TeV_pythia85_25nsV7-ak4.root");
    //TFile *h = new TFile("../MC.ak4R02.cmsJobs.P8Nuples.10:33:36-07-05-2016/P8-TuneCUETP8M_13TeV_pythia85_25nsV7-ak4.root");
    TFile *h = new TFile("../MC.ak7.cmsJobs.P8Nuples.18:26:28-08-15-2016/P8-TuneCUETP8M_13TeV_pythia85_25nsV7-withJECs-AK7-pthatCut.root");
    //TFile *h = new TFile("../MC.ak7ForDas.cmsJobs.P8.12:18:44-09-12-2016/P8-TuneCUETP8M_13TeV_pythia85_25nsV7-withJECs-AK7-pthatCut.root");
    //TFile *h = new TFile("../MC.ak4.cmsJobs.P8Nuples.18:48:29-06-14-2016/P8-TuneCUETP8M_13TeV_pythia85_25nsV7-ak4.root");
    //int Ybins[8] = {0,1,2,3,4,5,6,7};
    //TH1F *y[7];
    TH1F *hGen[7];
    TH1F *hMiss[7];
    TH1F *hReco[7];
    TH1F *hFake[7];
    TH2F *hMatch[7];
    RooUnfoldResponse response[7];
    for (int k=0; k < 7; k++){
        std::string Gen ("efficiency/pt_GENInclJet_");
        std::string GenMiss ("efficiency/Gen_MissInclusiveJets_");
        std::string Reco("efficiency/pt_DETInclJet_");
        std::string RecoFake ("efficiency/PF_FakeInclusiveJets_");
        std::string Match ("efficiency/TwoD_MatchedInclusiveJets_");
        std::string binN ("bin");
        std::string N ("");
        ostringstream convert;
        convert << k+1;
        N = convert.str();

         
        std::string totM, totG, totGM, totR, totRF, totMat;
        //totM = Mea + N + binN;
        totG = Gen + N + binN;
        totGM = GenMiss + N + binN;
        totR = Reco + N + binN;
        totRF = RecoFake + N + binN;
        totMat = Match + N + binN;


        hGen[k] = (TH1F*)gDirectory->Get(totG.c_str());
        hMiss[k] = (TH1F*)gDirectory->Get(totGM.c_str());
        hReco[k] = (TH1F*)gDirectory->Get(totR.c_str());
        hFake[k] = (TH1F*)gDirectory->Get(totRF.c_str());
        hMatch[k] = (TH2F*)gDirectory->Get(totMat.c_str());
        
        
        // Fill response matrix from MC

        response[k].Setup(hReco[k],hGen[k],hMatch[k]);

        /*
        for (int i=1; i < hReco[k]->GetXaxis()->GetNbins() + 1; i++)
            response[k].Fill(hReco[k]->GetBinContent(i),hGen[k]->GetBinContent(i));
        

         
        for (int j=1; j < hMiss[k]->GetXaxis()->GetNbins() + 1; j++)
            response[k].Miss(hMiss[k]->GetBinContent(j));
         
         
        for (int l=1; l < hFake[k]->GetXaxis()->GetNbins() + 1; l++)
            response[k].Fake(hFake[k]->GetBinContent(l));
        */
    }

    
    int iteration = 4;

    RooUnfoldBayes unfoldY1 (&response[0], y1, iteration);
    RooUnfoldBayes unfoldY2 (&response[1], y2, iteration);
    RooUnfoldBayes unfoldY3 (&response[2], y3, iteration);
    RooUnfoldBayes unfoldY4 (&response[3], y4, iteration);
    RooUnfoldBayes unfoldY5 (&response[4], y5, iteration);
    RooUnfoldBayes unfoldY6 (&response[5], y6, iteration);
    RooUnfoldBayes unfoldY7 (&response[6], y7, iteration);
   
    RooUnfold::ErrorTreatment f = RooUnfold::kCovariance;
    TH1D *recy1= (TH1D*) unfoldY1.Hreco(f);
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
      
    
    //TFile *out = new TFile("Unfolded.DATA2016BC-ak7-withP8ResponseAK8Triggers-iter5.root","RECREATE");
    TFile *out = new TFile("Unfolded.DATA2016BC-ak7-withP8ResponseAK8Triggers-iter5.root","RECREATE");
    //TFile *out = new TFile("Unfolded.DATA2015-ak7-withP8Response-iter5.root","RECREATE");
    //TFile *out = new TFile("Unfolded.DATAforDAS-ak7-withP8Response-iter5.root","RECREATE");
    recy1->Write();
    //response.Write();
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
    //m->Close();
    /*
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->SetLogy();
    c1->SetLogx();
    recy1->Draw();
    y1->SetLineColor(kRed);
    y1->Draw("same");
    */
}








