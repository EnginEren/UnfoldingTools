#include <vector>
// Root
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TColor.h>
// RooUnfold
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfold.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldBinByBin.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldErrors.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldParms.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldSvd.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldBayes.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldDagostini.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldInvert.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldResponse.h"
#include "/afs/desy.de/user/c/connorpa/Libraries/RooUnfold/src/RooUnfoldTUnfold.h"

TCanvas * unfolding (const double mu = 0, const double sigma = 0.01, const unsigned nevents = 1e4)
{
    TColor::CreateColorWheel();
    cout << "TOY UNFOLDING - WELCOME\n"
         << "\tmu=" << mu
         << "\tsigma=" << sigma 
         << "\tnevents=" << nevents << endl;

    cout << "=== Defining resolution function (used only if sigma different from 0)" << endl;
    TF1 * f_resolution = new TF1 ("resolution", "gaus(0)", -1, 1);
    f_resolution->SetParameters(1, mu, sigma);

    cout << "=== Defining histograms" << endl;
    vector<double> binning = {0,18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032/*,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3450,3637,3832,4037,4252,4477,4713,4961,5220,5492,5777,6076,6389,6717,7000*/};
    // --> standard pt binning
    const unsigned minpt = binning.front(), maxpt = binning.back();
    TH1D * h_gen = new TH1D ("gen", "truth", binning.size()-1, &binning[0]),
         * h_rec = new TH1D ("rec", "measurement", binning.size()-1, &binning[0]);
    TH2D * h_RM = new TH2D ("RM", "RM", binning.size()-1, &binning[0], binning.size()-1, &binning[0]);
    TH2D * h_resolution = new TH2D("resolution", "resolution", 41, -1, 1, binning.size()-1, &binning[0]);

    cout << "=== Filling histograms" << endl;
    TRandom3 r;
    cout << ">> RM and resolution" << endl;
    for (unsigned i = 0 ; i < nevents ; i++)
    {
        const float pt_gen = r.Rndm()*(maxpt-minpt),
                    xsec = pow(pt_gen, -4);
        for (unsigned j = 0 ; j < nevents ; j++)
        {
            float pt_rec, x, resolution;
            if (sigma == 0) // then we want a diagonal matrix
            {
                x = mu;
                pt_rec = pt_gen*(1-x);
                resolution = 1;
            }
            else // otherwise we use the resolution to get the right smearing
            {
                pt_rec = r.Rndm()*(maxpt-minpt);
                x = (pt_gen-pt_rec)/pt_gen;
                resolution = f_resolution->Eval(x);
            }
            h_RM        ->Fill(pt_rec, pt_gen, resolution*xsec);
            h_resolution->Fill(x     , pt_gen, resolution*xsec);
        }
    }

    cout << ">> truth and measurement" << endl;
    for (unsigned i = 0 ; i < nevents ; i++)
    {
        // truth
        float pt_gen = r.Rndm()*(maxpt-minpt),
              xsec = pow(pt_gen,-4);
        h_gen->Fill(pt_gen, xsec);
        // measurement
        float x, pt_rec;
        if (sigma == 0)
            x = mu;
        else
            x = f_resolution->GetRandom(-1,1);
        pt_rec = pt_gen*(1-x);
        h_rec->Fill(pt_rec, xsec);
    }

    cout << "=== Plotting RM" << endl;
    TCanvas * c = new TCanvas (TString::Format("unfolding%f%f", mu, sigma));

    c->Divide(2,2);

    c->cd(1);
    c->GetPad(1)->SetLogx();
    c->GetPad(1)->SetLogy();
    c->GetPad(1)->SetLogz();
    // get nonnegative minimum
    double minimum = h_RM->GetMaximum();
    for (unsigned short ibin = 1 ; ibin <= h_RM->GetNbinsX() ; ibin++)
        for (unsigned short jbin = 1 ; jbin <= h_RM->GetNbinsY() ; jbin++)
        {
            const double current_content = h_RM->GetBinContent(ibin,jbin);
            if (current_content > 0 && current_content < minimum)
                minimum = current_content;
        }
    h_RM->SetMinimum(minimum);
    h_RM->GetXaxis()->SetNoExponent();
    h_RM->GetXaxis()->SetMoreLogLabels();
    h_RM->GetYaxis()->SetNoExponent();
    h_RM->GetYaxis()->SetMoreLogLabels();
    h_RM->SetStats(0);
    h_RM->DrawCopy("colz");

    cout << "=== Plotting resolution" << endl;
    c->cd(3);
    c->GetPad(3)->SetLogy();
    c->GetPad(3)->SetLogz();
    // get nonnegative minimum
    minimum = h_resolution->GetMaximum();
    for (unsigned short ibin = 1 ; ibin <= h_resolution->GetNbinsX() ; ibin++)
        for (unsigned short jbin = 1 ; jbin <= h_resolution->GetNbinsY() ; jbin++)
        {
            const double current_content = h_resolution->GetBinContent(ibin,jbin);
            if (current_content > 0 && current_content < minimum)
                minimum = current_content;
        }
    h_resolution->GetYaxis()->SetNoExponent();
    h_resolution->GetYaxis()->SetMoreLogLabels();
    h_resolution->SetStats(0);
    h_resolution->DrawCopy("colz");

    cout << "=== Plotting spectra" << endl;
    c->cd(2);
    c->GetPad(2)->SetLogx();
    c->GetPad(2)->SetLogy();

    // the following vector will contain all the spectra at hadron level: the generated spectrum and the unfolded spectra
    h_gen->Rebin(2);
    h_RM->Rebin2D(1,2);
    vector<TH1 *> numerators;
    numerators.push_back(h_gen);

    RooUnfoldResponse * RU_response;
    RooUnfold * RU_unfolding;
    // bin by bin
    cout << ">> Bin-by-bin unfolding" << endl;
    RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("recbinbybin")), static_cast<TH1D*>(h_gen->Clone("genbinbybin")), static_cast<TH2D*>(h_RM->Clone("RMbinbybin")));
    RU_unfolding = new RooUnfoldBinByBin (RU_response, h_rec, "binbybin");
    numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("BinByBin")));
    numerators.back()->SetTitle("BinByBin");
    // Bayes
    cout << ">> Bayes unfolding" << endl;
    RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("recbayes")), static_cast<TH1D*>(h_gen->Clone("genbayes")), static_cast<TH2D*>(h_RM->Clone("RMbayes")));
    RU_unfolding = new RooUnfoldBayes (RU_response, h_rec, 3, "bayes");
    numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("Bayes")));
    numerators.back()->SetTitle("Bayes");
    //// D'Agostini
    //cout << ">> D'Agostini unfolding" << endl;
    //RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("recdagostini")), static_cast<TH1D*>(h_gen->Clone("gendagostini")), static_cast<TH2D*>(h_RM->Clone("RMdagostini")));
    //RU_unfolding = new RooUnfoldDagostini (RU_response, h_rec, 3, "dagostini");
    //numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("DAgostini")));
    //numerators.back()->SetTitle("DAgostini");
    // SVD
    cout << ">> SVD unfolding" << endl;
    RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("recsvd")), static_cast<TH1D*>(h_gen->Clone("gensvd")), static_cast<TH2D*>(h_RM->Clone("RMsvd")));
    RU_unfolding = new RooUnfoldSvd (RU_response, h_rec, 20, "svd");
    numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("SVD")));
    numerators.back()->SetTitle("SVD");
    // inversion
    cout << ">> matrix-inversion unfolding" << endl;
    RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("recinversion")), static_cast<TH1D*>(h_gen->Clone("geninversion")), static_cast<TH2D*>(h_RM->Clone("RMinversion")));
    RU_unfolding = new RooUnfoldInvert (RU_response, h_rec, "inv");
    numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("Inverse")));
    numerators.back()->SetTitle("Inverse");
    // TUnfold 
    cout << ">> TUnfold unfolding" << endl;
    RU_response = new RooUnfoldResponse (static_cast<TH1D*>(h_rec->Clone("rectu")), static_cast<TH1D*>(h_gen->Clone("gentu")), static_cast<TH2D*>(h_RM->Clone("RMtu")));
    RU_unfolding = new RooUnfoldTUnfold (RU_response, h_rec, TUnfold::kRegModeDerivative, "tu");
    numerators.push_back(static_cast<TH1D*>(RU_unfolding->Hreco()->Clone("TU")));
    numerators.back()->SetTitle("TUnfold");
    
    h_rec->SetLineColor(kBlack);
    //h_rec->SetBinContent(0, 0);
    //h_rec->SetBinContent(1+h_rec->GetNbinsX(), 0);
    h_rec->SetStats(0);
    h_rec->SetTitle("measurement");
    h_rec->GetYaxis()->SetRangeUser(1e-10,1e-1);
    h_rec->SetMarkerStyle(20);
    h_rec->SetMarkerSize(0.5);
    h_rec->Rebin(2);
    h_rec->GetXaxis()->SetNoExponent();
    h_rec->GetXaxis()->SetMoreLogLabels();
    h_rec->DrawCopy();

    for (unsigned short i = 0 ; i < numerators.size() ; i++)
    {
        numerators[i]->SetLineColor(2+i);
        //numerators[i]->SetBinContent(0, 0);
        //numerators[i]->SetBinContent(1+numerators[i]->GetNbinsX(), 0);
        numerators[i]->SetStats(0);
        numerators[i]->DrawCopy("same");
    }
    c->GetPad(2)->BuildLegend();

    cout << "=== Plotting ratios" << endl;
    c->cd(4);
    c->GetPad(4)->Divide(1,numerators.size(), 0,0);
    for (unsigned short i = 0 ; i < numerators.size() ; i++)
    {
        c->GetPad(4)->cd(i+1);
        c->GetPad(4)->GetPad(i+1)->SetLogx();
        c->GetPad(4)->GetPad(i+1)->SetGridx();
        c->GetPad(4)->GetPad(i+1)->SetGridy();
        TH1D * ratio = static_cast<TH1D*>(numerators[i]->Clone(TString::Format("%s/rec", numerators[i]->GetName())));
        ratio->Divide(h_rec);
        ratio->SetTitle(TString::Format(";;%s", ratio->GetName()));
        ratio->SetTitleSize(2);
        ratio->SetStats(0);
        ratio->GetYaxis()->SetRangeUser(0.4,1.6);
        ratio->DrawCopy();
    }

    // esthetics
    cout << "=== Redrawing axes" << endl;
    for (unsigned short i = 1 ; i <= 3 ; i++)
        c->GetPad(i)->RedrawAxis();
    for (unsigned short i = 1 ; i <= numerators.size() ; i++)
        c->GetPad(4)->GetPad(i)->RedrawAxis();

    cout << "=== The end" << endl;
    return c;
}

