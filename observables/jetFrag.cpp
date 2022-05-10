// #include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <dirent.h>
using namespace std;

#include "ObsHelper.h"
#include "CommandLine.h"
#include "PlotHelper4.h"
#include "SetStyle.h"
#include "uti.h"
#include <math.h>
#include <set>

// NOTE: this is export for jupyter notebook

// g++ jetFrag.cpp $(root-config --cflags --libs) -O2 -o "jetFrag.exe"
// ./jetFrag.exe -input "/data/kying/final80/pp80,/data/kying/final80/ppPhotonJet80" -minpt 80 -maxpt 100
// ./jetFrag.exe -input "/data/kying/final80/PbPbWide80_0_10,/data/kying/final80/PbPbWidePhotonJet80_0_10" -minpt 80 -maxpt 100

int main(int argc, char *argv[])
{
    SetThesisStyle();

    CommandLine CL(argc, argv);

    vector<string> InputFileFolders = CL.GetStringVector("input");
    double MinPT = CL.GetDouble("minpt", 0);
    double MaxPT = CL.GetDouble("maxpt", 500);
    // string kappasFile = CL.Get("kappas");
    double trackPtCut = 0;

    int NBins = 10; // let's just do 1-5
    const Int_t NTopics = 2;
    TH1D *ZHists[NTopics * 3];
    TH1D *XiHists[NTopics * 3];
    double yMin = INFINITY;
    double yMax = -INFINITY;

    for (int i = 0; i < NTopics * 3; i++)
    {
        // TODO make sure axes labels are right
        XiHists[i] = new TH1D("Xi", "Jet Fragmentation;\\xi;1/Nd\\xi", 10, 0, 5);
        XiHists[i]->SetStats(0);
        ZHists[i] = new TH1D("Z", "Jet Fragmentation;z;1/Ndz", 10, 0, 1);
        ZHists[i]->SetStats(0);
    }

    bool pp = ispp(InputFileFolders[0]);
    set<Int_t> chargedIds = {11, 13, 15, 17, 24, 34, 37,
                             211, 213, 215, 217, 219,
                             321, 323, 325, 327, 329,
                             421, 423, 425, 431, 433, 435,
                             521, 523, 525, 541, 543, 545};
    int totalJetsGluon = 0;
    int totalJetsQuark = 0;
    int totalJets[2] = {0, 0};
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

                    int NParticles = (*SignalJetConstituentPt)[iJ].size();

                    int matElemIdx;
                    double matElemDr;
                    tie(matElemIdx, matElemDr) = getMatElemIndexDr(dijet, SignalJetMatrixElem, SignalJetMatrixElemPhi, SignalJetMatrixElemEta, (*SignalJetPhi)[iJ], (*SignalJetEta)[iJ]);
                    int matElem = (*SignalJetMatrixElem)[matElemIdx];

                    for (int iP = 0; iP < NParticles; iP++)
                    {
                        // if particle is not charged then continue
                        if (chargedIds.find((*SignalJetConstituentPid)[iJ][iP]) == chargedIds.end())
                            continue;

                        double r = getR((*SignalJetConstituentEta)[iJ][iP], (*SignalJetEta)[iJ], (*SignalJetConstituentPhi)[iJ][iP], (*SignalJetPhi)[iJ]);
                        double z = (*SignalJetConstituentPt)[iJ][iP] / (*SignalJetPt)[iJ] * cos(r);
                        double xi = -log(z);
                        XiHists[f]->Fill(xi, (*Weight)[0]);
                        ZHists[f]->Fill(z, (*Weight)[0]);
                        if (matElemDr < 0.4)
                        {
                            if (matElem == 21)
                            {
                                // gluon
                                XiHists[3]->Fill(xi, (*Weight)[0]);
                                ZHists[3]->Fill(z, (*Weight)[0]);
                            }
                            else
                            {
                                XiHists[2]->Fill(xi, (*Weight)[0]);
                                ZHists[2]->Fill(z, (*Weight)[0]);
                            }
                        }
                    }

                    if (matElemDr < 0.4)
                    {
                        totalJets[f] += (*Weight)[0];
                        if (matElem == 21)
                        {
                            totalJetsGluon += (*Weight)[0];
                        }
                        else
                        {
                            totalJetsQuark += (*Weight)[0];
                        }
                    }
                }
            }
            File.Close();
            // csvFile << "numJets_" + folder << "," << totalJets[f] << "\n";
        }
        XiHists[f]->Scale(1. / totalJets[f]);
        ZHists[f]->Scale(1. / totalJets[f]);

        for (int bin = 1; bin <= 10; bin++)
        {
            if (dijet)
            {
                csvFile << "xi_dijet," << -0.25 + bin * 0.5 << "," << 0.25 << "," << XiHists[f]->GetBinContent(bin) << "," << XiHists[f]->GetBinError(bin) << "\n";
                csvFile << "z_dijet," << -0.25 + bin * 0.5 << "," << 0.25 << "," << ZHists[f]->GetBinContent(bin) << "," << ZHists[f]->GetBinError(bin) << "\n";
            }
            else
            {
                csvFile << "xi_photonjet," << -0.25 + bin * 0.5 << "," << 0.25 << "," << XiHists[f]->GetBinContent(bin) << "," << XiHists[f]->GetBinError(bin) << "\n";
                csvFile << "z_photonjet," << -0.25 + bin * 0.5 << "," << 0.25 << "," << ZHists[f]->GetBinContent(bin) << "," << ZHists[f]->GetBinError(bin) << "\n";
            }
        }
    }

    XiHists[2]->Scale(1. / totalJetsQuark);
    XiHists[3]->Scale(1. / totalJetsGluon);
    ZHists[2]->Scale(1. / totalJetsQuark);
    ZHists[3]->Scale(1. / totalJetsGluon);

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
        for (int bin = 1; bin <= 10; bin++)
        {
            csvFile << "xi_" + label << "," << -0.25 + bin * 0.5 << "," << 0.25 << "," << XiHists[h]->GetBinContent(bin) << "," << XiHists[h]->GetBinError(bin) << "\n";
            csvFile << "z_" + label << "," << -0.25 + bin * 0.5 << "," << 0.25 << "," << ZHists[h]->GetBinContent(bin) << "," << ZHists[h]->GetBinError(bin) << "\n";
        }
    }

    return 0;
}
