#include <iostream>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <dirent.h>
using namespace std;

#include "CommandLine.h"
#include "PlotHelper4.h"
#include "SetStyle.h"
#include "uti.h"
// #include "LabelHelper.h"

// g++ getTopicModelingCsv.cpp $(root-config --cflags --libs) -O2 -o "getTopicModelingCsv.exe"
// ./getTopicModelingCsv.exe -input "pp_pt80100_trkpt0:1:4:10_dijet.root"
// ./getTopicModelingCsv.exe -input /data/kying/lda/TMVAObsApplication/lda_PbPbWide_pt100120_photonjet_quark_truth_cwola.root -labels pbpb80_0_10_wide_photonjet_quark_truth -output /data/kying/lda/TMVAObsCsvs/lda_pt100120_cwola.csv

void writeToCsv(TH1D *multiplicityHist, string label, ofstream &outFile)
{
    vector<double> v, vError;

    for (int i = 0; i < multiplicityHist->GetNbinsX(); i++)
    {
        v.push_back(multiplicityHist->GetBinContent(i + 1));
        vError.push_back(multiplicityHist->GetBinError(i + 1));
    }
    outFile << label << ",";
    for (const auto &e : v)
        outFile << e << ",";
    outFile << endl;

    outFile << label << "_error"
            << ",";
    for (const auto &e : vError)
        outFile << e << ",";
    outFile << endl;
}

int main(int argc, char *argv[])
{
    SetThesisStyle();

    CommandLine CL(argc, argv);

    vector<string> InputFileNames = CL.GetStringVector("input");
    vector<string> InputFileLabels = CL.GetStringVector("labels");
    string outputName = CL.Get("output");

    // stuff for the output csv
    ofstream outFile(outputName.c_str(), ios_base::app);

    for (int i = 0; i < InputFileNames.size(); i++)
    {
        string fileName = InputFileNames[i];
        string label = InputFileLabels[i];
        TFile File(fileName.c_str());

        TTree *Tree = (TTree *)File.Get("tree");

        if (Tree == nullptr)
        {
            File.Close();
            continue;
        }

        // TODO LDA change here
        // TH1D Histogram("Hist", "Title;X;Y", 350, -2.5, 2.5);
        TH1D Hist("Hist", "Title;X;Y", 200, 0, 1);
        // Tree->Draw("LDA_response>>Hist");
        Tree->Draw("MLP_response>>Hist");

        TCanvas *c = new TCanvas("c", "canvas", 800, 800);
        Hist.Draw();
        c->Update();
        c->SaveAs((TString)("plots/" + InputFileLabels[i] + ".jpg"));

        writeToCsv(&Hist, label, outFile);

        File.Close();
    }

    return 0;
}
