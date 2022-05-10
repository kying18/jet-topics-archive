#include <iostream>
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

// Note: this is export for jupyter notebook

// g++ jetShape.cpp $(root-config --cflags --libs) -O2 -o "jetShape.exe"
// ./jetShape.exe -input "/data/kying/final80/pp80,/data/kying/final80/ppPhotonJet80" -minpt 80 -maxpt 100 -trkpt 0
// ./jetShape.exe -input "/data/kying/final80/PbPbWide80_0_10,/data/kying/final80/PbPbWidePhotonJet80_0_10" -minpt 80 -maxpt 100 -trkpt 0

// tracks are charged particles
// use particleseta, particlesphi, signaljeteta, signaljetphi
// dphi needs to be between -pi and pi
// jet shape = 1/dr 1/Njet sum[[pT_track for track in tracks if track in dr bucket] / pT_jet for jet in jets]
// in the example, they used dr = 0.05

// plot the sum pT distribution from tracks
// as a function of delta R with respect to the jet
// and delta R is defined as dR = sqrt(delta eta ^2 + delta phi ^2)
// delta eta = track eta - jet eta
// delta phi = track phi - jet phi
// (of course, need to be careful here about the delta phi calculation, don’t miss the factor of 2 pi here)
// we want the delta phi to be between -pi and pi
// runEmmi should be able to do it
// so jet shape is like
// the constituent momentum spectrum as a function of dR

int main(int argc, char *argv[])
{
    SetThesisStyle();

    CommandLine CL(argc, argv);

    vector<string> InputFileFolders = CL.GetStringVector("input");
    double MinPT = CL.GetDouble("minpt", 0);
    double MaxPT = CL.GetDouble("maxpt", 500);
    double trackPtCut = CL.GetDouble("trkpt", 1.0);
    double dr = CL.GetDouble("dr", 0.05);     // bins are [0, 0.05), [0.05, 0.10), [0.10, 0.15), ... [0.25, 0.30)
    double MaxR = CL.GetDouble("maxr", 0.30); // max value of max on rightmost bin

    int NBins = (int)ceil(MaxR / dr);
    Int_t NTopics = 2;
    Double_t x[NBins], y[NTopics * 3][NBins]; // 01 original, 23 q/g truth, 45 topics
    Double_t ex[NBins], ey[NTopics * 3][NBins];
    string legendLabels[NTopics * 3];
    double yMin = INFINITY;
    double yMax = -INFINITY;

    bool pp = ispp(InputFileFolders[0]);

    ofstream csvFile;
    string name = "data_output_path/" + to_string((int)MinPT) + "_trackpt" + to_string((int)trackPtCut) + (pp ? "_pp" : "_pbpb") + ".csv";
    csvFile.open(name, ios_base::app);

    int count[NTopics + 1][NBins]; // 0 will be for the folder, 1 quark, 2 gluon
    double rho[NTopics + 1][NBins];
    double rhoErrors[NTopics + 1][NBins];
    for (int i = 0; i < NTopics + 1; i++)
    {
        fill(count[i], count[i] + NBins, 0);
        fill(rho[i], rho[i] + NBins, 0);
        fill(rhoErrors[i], rhoErrors[i] + NBins, 0);
    }

    for (int f = 0; f < NTopics; f++)
    {
        string folder = InputFileFolders[f];
        cout << folder << endl;

        bool dijet = isDijet(folder);
        // bool pp = ispp(folder);

        struct dirent *de; // to get files in path
        DIR *dir = opendir(folder.c_str());
        if (dir == NULL) // opendir returns NULL if couldn't open directory
        {
            cout << "Could not open directory: " << folder << endl;
            return 0;
        }

        // INITIALIZING THE BUCKETS FOR JET SHAPE
        // int count [NBins];
        // double rho [NBins];
        fill(count[0], count[0] + NBins, 0);
        fill(rho[0], rho[0] + NBins, 0);

        int fileCount = 0;
        while ((de = readdir(dir)) != NULL)
        {
            if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, ".."))
                continue;
            // if (fileCount >= 1) continue; // REMOVE
            string fileName = de->d_name;
            TFile File((folder + "/" + fileName).c_str());

            // cout << "Filename: " << fileName << endl;

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
            // cout << "Entry count... " << to_string(EntryCount) << endl;
            for (int iE = 0; iE < EntryCount; iE++)
            {
                Tree->GetEntry(iE);

                if (SignalJetPt == nullptr)
                    continue;

                // get max pt jet indices for the event
                int NJet = SignalJetPt->size();
                // cout << "Jets in entry: " << NJet << endl;
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

                    double sumTrackPts[NTopics + 1][NBins];
                    double sumTrackPtsErrors[NTopics + 1][NBins];
                    for (int i = 0; i < NTopics + 1; i++)
                    {
                        fill(sumTrackPts[i], sumTrackPts[i] + NBins, 0);
                        fill(sumTrackPtsErrors[i], sumTrackPtsErrors[i] + NBins, 0);
                    }

                    int NParticles = (*SignalJetConstituentPt)[iJ].size();

                    // int matElem = (*SignalJetMatrixElem)[iJ]; // todo this is wrong!! need to calculate the matrix elem based on getmatelemindexdr
                    int matElemIdx;
                    double matElemDr;
                    tie(matElemIdx, matElemDr) = getMatElemIndexDr(dijet, SignalJetMatrixElem, SignalJetMatrixElemPhi, SignalJetMatrixElemEta, (*SignalJetPhi)[iJ], (*SignalJetEta)[iJ]);
                    int matElem = (*SignalJetMatrixElem)[matElemIdx];

                    for (int iP = 0; iP < NParticles; iP++)
                    {
                        if ((*SignalJetConstituentPt)[iJ][iP] > trackPtCut)
                        {
                            double r = getR((*SignalJetConstituentEta)[iJ][iP], (*SignalJetEta)[iJ], (*SignalJetConstituentPhi)[iJ][iP], (*SignalJetPhi)[iJ]);

                            for (int bin = 0; bin < NBins; bin++)
                            {
                                if ((r >= bin * dr) & (r < (bin + 1) * dr))
                                {
                                    if (matElemDr < 0.4)
                                    {
                                        if (matElem == 21)
                                        {
                                            // gluon
                                            sumTrackPts[2][bin] += (*SignalJetConstituentPt)[iJ][iP];
                                            sumTrackPtsErrors[2][bin] += pow((*SignalJetConstituentPt)[iJ][iP], 2);
                                        }
                                        else
                                        {
                                            sumTrackPts[1][bin] += (*SignalJetConstituentPt)[iJ][iP];
                                            sumTrackPtsErrors[1][bin] += pow((*SignalJetConstituentPt)[iJ][iP], 2);
                                        }
                                    }
                                    sumTrackPts[0][bin] += (*SignalJetConstituentPt)[iJ][iP];
                                    sumTrackPtsErrors[0][bin] += pow((*SignalJetConstituentPt)[iJ][iP], 2);
                                }
                            }
                        }
                    }

                    for (int bin = 0; bin < NBins; bin++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            sumTrackPtsErrors[i][bin] = sqrt(sumTrackPtsErrors[i][bin]) / (*SignalJetPt)[iJ];
                        }
                        if (matElemDr < 0.4)
                        {
                            if (matElem == 21)
                            {
                                // gluon
                                rho[2][bin] += sumTrackPts[2][bin] * (*Weight)[0] / (*SignalJetPt)[iJ];
                                rhoErrors[2][bin] += pow(sumTrackPtsErrors[2][bin], 2);
                                count[2][bin] += (*Weight)[0];
                            }
                            else
                            {
                                rho[1][bin] += sumTrackPts[1][bin] * (*Weight)[0] / (*SignalJetPt)[iJ];
                                rhoErrors[1][bin] += pow(sumTrackPtsErrors[1][bin], 2);
                                count[1][bin] += (*Weight)[0];
                            }
                        }
                        rho[0][bin] += sumTrackPts[0][bin] * (*Weight)[0] / (*SignalJetPt)[iJ];
                        rhoErrors[0][bin] += pow(sumTrackPtsErrors[0][bin], 2);
                        count[0][bin] += (*Weight)[0];
                    }
                }
            }
            File.Close();
        }
        for (int bin = 0; bin < NBins; bin++)
        {
            rho[0][bin] *= 1 / dr * 1 / count[0][bin];
            cout << (bin + 0.5) * dr << " " << rho[0][bin] << endl;
            x[bin] = (bin + 0.5) * dr;
            ex[bin] = 0.5 * dr;
            y[f][bin] = rho[0][bin];
            ey[f][bin] = 1 / dr * 1 / count[0][bin] * sqrt(rhoErrors[0][bin]);

            if (dijet)
            {
                csvFile << "dijet"
                        << "," << x[bin] << "," << ex[bin] << "," << y[f][bin] << "," << ey[f][bin] << "\n";
            }
            else
            {
                csvFile << "photonjet"
                        << "," << x[bin] << "," << ex[bin] << "," << y[f][bin] << "," << ey[f][bin] << "\n";
            }
        }
    }

    for (int i = 0; i < 2; i++)
    {
        for (int bin = 0; bin < NBins; bin++)
        {
            rho[i + 1][bin] *= 1 / dr * 1 / count[i + 1][bin];

            y[i + 2][bin] = rho[i + 1][bin];

            ey[i + 2][bin] = 1 / dr * 1 / count[i + 1][bin] * sqrt(rhoErrors[i + 1][bin]);

            if (i == 0)
                csvFile << "quark"
                        << ",";
            if (i == 1)
                csvFile << "gluon"
                        << ",";
            csvFile << x[bin] << "," << ex[bin] << "," << y[i + 2][bin] << "," << ey[i + 2][bin] << "\n";
        }
    }

    csvFile.close();

    return 0;
}