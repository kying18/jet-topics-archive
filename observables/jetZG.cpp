// #include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <dirent.h>
#include <math.h>
#include <set>
using namespace std;

#include "ObsHelper.h"
#include "CommandLine.h"
#include "PlotHelper4.h"
#include "SetStyle.h"
#include "uti.h"

// NOTE: this is export for jupyter notebook
// 0 is dijet, 1 is photonjet

// jet ZG definition according to https://reader.elsevier.com/reader/sd/pii/S037026931730936X?token=776D19B17AE19EAF3D411B4E1B1842651E46111417A3A2AD9DD99B64F97D5B7DDC97B8EC1ADB9129D5ED6FFBB1AF6662&originRegion=us-east-1&originCreation=20211213182758
// m = sqrt(E^2 - pT^2 - pZ^2)

// g++ jetZG.cpp $(root-config --cflags --libs) -O2 -o "jetZG.exe"
// ./jetZG.exe -input "/data/kying/final80/pp80,/data/kying/final80/ppPhotonJet80" -minpt 80 -maxpt 100
// ./jetZG.exe -input "/data/kying/final80/PbPbWide80_0_10,/data/kying/final80/PbPbWidePhotonJet80_0_10" -minpt 80 -maxpt 100

int main(int argc, char *argv[])
{
    SetThesisStyle();

    CommandLine CL(argc, argv);

    vector<string> InputFileFolders = CL.GetStringVector("input");
    double MinPT = CL.GetDouble("minpt", 0);
    double MaxPT = CL.GetDouble("maxpt", 500);
    double trackPtCut = 0;

    // int NBins = 10; // let's just do 1-5
    const Int_t NTopics = 2;
    TH1D *ZGHists[NTopics * 2];
    double yMin = INFINITY;
    double yMax = -INFINITY;

    for (int i = 0; i < NTopics * 2; i++)
    {
        // TODO make sure axes labels are right
        ZGHists[i] = new TH1D("ZG", "Jet ZG;M;1/N dN/dM", 8, 0.1, 0.5);
        ZGHists[i]->SetStats(0);
    }

    bool pp = ispp(InputFileFolders[0]);
    ofstream csvFile;
    string name = "data_output_path/" + to_string((int)MinPT) + "_trackpt" + to_string((int)trackPtCut) + (pp ? "_pp" : "_pbpb") + ".csv";
    csvFile.open(name);

    for (int f = 0; f < NTopics; f++)
    {
        string folder = InputFileFolders[f];

        bool dijet = isDijet(folder);

        struct dirent *de; // to get files in path
        DIR *dir = opendir(folder.c_str());
        if (dir == NULL) // opendir returns NULL if couldn't open directory
        {
            cout << "Could not open directory: " << folder << endl;
            return 0;
        }

        int fileCount = 0;
        // int totalJets = 0;
        while ((de = readdir(dir)) != NULL)
        {
            if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, ".."))
                continue;
            // if (fileCount >= 5) continue; // REMOVE
            string fileName = de->d_name;
            TFile File((folder + "/" + fileName).c_str());

            TTree *Tree = (TTree *)File.Get("JetTree");

            if (Tree == nullptr)
            {
                File.Close();
                continue;
            }

            fileCount += 1; // REMOVE

            vector<double> *SignalJetNConstituent = nullptr;
            vector<int> *SignalJetMatrixElem = nullptr;
            vector<double> *SignalJetMatrixElemPhi = nullptr;
            vector<double> *SignalJetMatrixElemEta = nullptr;
            vector<double> *SignalJetPt = nullptr;
            vector<double> *Weight = nullptr;
            vector<double> *SignalJetPhi = nullptr;
            vector<double> *SignalJetEta = nullptr;
            vector<double> *SignalJetZG = nullptr;
            vector<double> *LeadingZPhi = nullptr;
            vector<double> *LeadingPhotonPhi = nullptr;
            vector<vector<double>> *SignalJetConstituentPt = nullptr;
            vector<vector<double>> *SignalJetConstituentEta = nullptr;
            vector<vector<double>> *SignalJetConstituentPhi = nullptr;
            vector<vector<double>> *SignalJetConstituentPid = nullptr;

            Tree->SetBranchAddress("SignalJetNConstituent", &SignalJetNConstituent);
            Tree->SetBranchAddress("SignalJetMatrixElem", &SignalJetMatrixElem);
            Tree->SetBranchAddress("SignalJetMatrixElemPhi", &SignalJetMatrixElemPhi);
            Tree->SetBranchAddress("SignalJetMatrixElemEta", &SignalJetMatrixElemEta);
            Tree->SetBranchAddress("SignalJetPt", &SignalJetPt);
            Tree->SetBranchAddress("SignalJetPhi", &SignalJetPhi);
            Tree->SetBranchAddress("SignalJetEta", &SignalJetEta);
            Tree->SetBranchAddress("SignalJetSD1ZG", &SignalJetZG);
            Tree->SetBranchAddress("LeadingZPhi", &LeadingZPhi);
            Tree->SetBranchAddress("LeadingPhotonPhi", &LeadingPhotonPhi);
            Tree->SetBranchAddress("EventWeight", &Weight);
            Tree->SetBranchAddress("SignalJetConstituentPt", &SignalJetConstituentPt);
            Tree->SetBranchAddress("SignalJetConstituentEta", &SignalJetConstituentEta);
            Tree->SetBranchAddress("SignalJetConstituentPhi", &SignalJetConstituentPhi);
            Tree->SetBranchAddress("SignalJetConstituentPid", &SignalJetConstituentPid);

            int EntryCount = Tree->GetEntries();

            for (int iE = 0; iE < EntryCount; iE++)
            {
                Tree->GetEntry(iE);

                if (SignalJetPt == nullptr)
                    continue;

                // get max pt jet indices for the event
                int NJet = SignalJetPt->size();
                pair<int, int> MaxJetPtIndices = getMaxPtJetIndices(NJet, dijet, (*LeadingPhotonPhi)[0], SignalJetPhi, SignalJetPt, SignalJetConstituentPt, SignalJetConstituentPid);
                pair<double, double> MaxJetPts((*SignalJetPt)[MaxJetPtIndices.first], (*SignalJetPt)[MaxJetPtIndices.second]);

                int jetIdxs[2] = {MaxJetPtIndices.first, MaxJetPtIndices.second};
                for (int iJ : jetIdxs)
                {
                    if (iJ == -1)
                        continue;
                    if ((*SignalJetPt)[iJ] < MinPT || (*SignalJetPt)[iJ] > MaxPT)
                        continue;
                    if (hasLargePhotonPtFrac((*SignalJetConstituentPt)[iJ], (*SignalJetConstituentPid)[iJ]))
                        continue;
                    if (fabs((*SignalJetEta)[iJ]) > 1)
                        continue; // if abs(eta) > 1 then continue

                    int matElemIdx;
                    double matElemDr;
                    tie(matElemIdx, matElemDr) = getMatElemIndexDr(dijet, SignalJetMatrixElem, SignalJetMatrixElemPhi, SignalJetMatrixElemEta, (*SignalJetPhi)[iJ], (*SignalJetEta)[iJ]);
                    int matElem = (*SignalJetMatrixElem)[matElemIdx];

                    if (matElemDr < 0.4)
                    {
                        if (matElem == 21)
                        {
                            ZGHists[3]->Fill((*SignalJetZG)[iJ], (*Weight)[0]);
                        }
                        else
                        {
                            ZGHists[2]->Fill((*SignalJetZG)[iJ], (*Weight)[0]);
                        }
                    }
                    ZGHists[f]->Fill((*SignalJetZG)[iJ], (*Weight)[0]);
                }
            }
            File.Close();
        }

        for (int bin = 1; bin <= 8; bin++)
        {
            if (dijet)
            {
                csvFile << "dijet"
                        << "," << 0.1 + bin * 0.05 - 0.025 << ",0.025," << ZGHists[f]->GetBinContent(bin) << "," << ZGHists[f]->GetBinError(bin) << "\n";
            }
            else
            {
                csvFile << "photonjet"
                        << "," << 0.1 + bin * 0.05 - 0.025 << ",0.025," << ZGHists[f]->GetBinContent(bin) << "," << ZGHists[f]->GetBinError(bin) << "\n";
            }
        }
    }

    for (int h = 2; h <= 3; h++)
    {
        string label;
        if (h == 2)
        {
            label = "quark";
        }
        else
        {
            label = "gluon";
        }
        for (int bin = 1; bin <= 8; bin++)
        {
            csvFile << label << "," << 0.1 + bin * 0.05 - 0.025 << ",0.025," << ZGHists[h]->GetBinContent(bin) << "," << ZGHists[h]->GetBinError(bin) << "\n";
        }
    }

    return 0;
}
