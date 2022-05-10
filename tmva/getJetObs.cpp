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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodLD.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Factory.h"

// g++ getJetObs.cpp $(root-config --cflags --libs) -O2 -o "getJetObs.exe"
// ./getJetObs.exe -input "/data/kying/final80/ppPhotonJet80,/data/kying/final80/pp80" -minpt 80 -maxpt 100 -output /data/kying/ml-export/pp_pt80100
// ./getJetObs.exe -input "/data/kying/final80/PbPbWidePhotonJet80_0_10,/data/kying/final80/PbPbWide80_0_10" -minpt 80 -maxpt 100 -output /data/kying/ml-export/pbpb_pt80100

bool oppositeHemisphere(double particlePhi, double jetPhi)
{
    if (particlePhi >= 0)
    {
        double dPhi = particlePhi - jetPhi;
        if (dPhi < -M_PI)
            dPhi = dPhi + 2 * M_PI;
        if (dPhi > +M_PI)
            dPhi = dPhi - 2 * M_PI;

        if (dPhi < -M_PI / 2 || dPhi > M_PI / 2)
        {
            return true;
        }
        return false;
    }
    // todo: what if there is no leading Z?
    return false; // false for now.. don't count it
}

double getR(double etaMatrix, double etaJet, double phiMatrix, double phiJet)
{
    // this is the reconstructed mat elem's radial distance from jet axis (defined by etaJet and phiJet)
    double deta = etaMatrix - etaJet;
    double dphi = phiMatrix - phiJet;

    // dphi should be between -pi and pi
    if (dphi < -M_PI)
    {
        dphi = 2 * M_PI - abs(dphi); // if phi is < -pi, then 2pi - value of dphi is new
    }
    else if (dphi > M_PI)
    {
        dphi = -(2 * M_PI - dphi); // if phi > pi, then 2pi - dphi is new value, but need to negate
    }

    // if (!(dphi >= -M_PI && dphi <= M_PI)) {
    //     cout << "dphi: " << dphi << "\n";
    // }
    assert(dphi >= -M_PI && dphi <= M_PI);

    return sqrt(pow(deta, 2) + pow(dphi, 2));
}

bool isDijet(string folderName)
{
    if (folderName.find("PhotonJet") != string::npos)
        return false;
    return true;
}

pair<int, int> getMaxPtJetIndices(int NJet, bool dijet, double LeadingPhotonPhi, vector<double> *SignalJetPhi, vector<double> *SignalJetPt, vector<vector<double>> *SignalJetConstituentPt, vector<vector<double>> *SignalJetConstituentPid)
{
    pair<int, int> MaxJetPtIndices(-1, -1);
    pair<double, double> MaxJetPts(0, 0);

    for (int iJ = 0; iJ < NJet; iJ++)
    {
        // if photon jet and jet is not in opposite hemisphere, then we skip this jet
        if (!dijet && !oppositeHemisphere(LeadingPhotonPhi, (*SignalJetPhi)[iJ]))
        {
            // cout << "leading photon phi: " << LeadingPhotonPhi << ", signal jet phi: " << (*SignalJetPhi)[iJ] << endl;
            continue;
        }
        // if (hasLargePhotonPtFrac((*SignalJetConstituentPt)[iJ], (*SignalJetConstituentPid)[iJ])) continue;

        // cout << "Signal Jet Pt: " << (*SignalJetPt)[iJ] << endl;

        if ((*SignalJetPt)[iJ] > MaxJetPts.first)
        {
            // move current first to second, insert current jet as first
            MaxJetPts.second = MaxJetPts.first;
            MaxJetPtIndices.second = MaxJetPtIndices.first;
            MaxJetPts.first = (*SignalJetPt)[iJ];
            MaxJetPtIndices.first = iJ;
        }
        else if ((*SignalJetPt)[iJ] > MaxJetPts.second)
        {
            MaxJetPts.second = (*SignalJetPt)[iJ];
            MaxJetPtIndices.second = iJ;
        }
    }
    return MaxJetPtIndices;
}

pair<int, double> getMatElemIndexDr(bool dijet, vector<int> *SignalJetMatrixElem, vector<double> *SignalJetMatrixElemPhi, vector<double> *SignalJetMatrixElemEta, double SignalJetPhi, double SignalJetEta)
{
    // returns the int rep index of the matrix elem in the vector and the double rep DR of the jet to the matrix elem
    // all of these should be size 2
    assert(SignalJetMatrixElem->size() == 2 && SignalJetMatrixElemPhi->size() == 2 && SignalJetMatrixElemEta->size() == 2);
    int index;
    double dr;

    if (!dijet)
    {
        // if photon jet, then return the matrix element that is not the photon
        index = ((*SignalJetMatrixElem)[0] != 22 ? 0 : 1);
        if ((*SignalJetMatrixElem)[1] != 22)
        {
            return make_pair(-1, -1);
        }
        dr = getR((*SignalJetMatrixElemEta)[index], SignalJetEta, (*SignalJetMatrixElemPhi)[index], SignalJetPhi);
        return make_pair(index, dr);
    }
    // otherwise if dijet, compute the r value and select the smaller
    // cout << "check 2" << "\n";
    double index_0_r = getR((*SignalJetMatrixElemEta)[0], SignalJetEta, (*SignalJetMatrixElemPhi)[0], SignalJetPhi);
    // cout << "check 3" << "\n";
    double index_1_r = getR((*SignalJetMatrixElemEta)[1], SignalJetEta, (*SignalJetMatrixElemPhi)[1], SignalJetPhi);

    index = (index_0_r <= index_1_r ? 0 : 1);
    dr = (index_0_r <= index_1_r ? index_0_r : index_1_r);
    return make_pair(index, dr);
}

bool hasLargePhotonPtFrac(vector<double> ConstituentPts, vector<double> ConstituentPids)
{
    assert(ConstituentPts.size() == ConstituentPids.size());
    // get the largest photon pt fraction and cut if larger than 0.8
    double MaxPhotonPt = 0;
    double TotalPt = 0;
    double CurrentPt;

    for (int i = 0; i < ConstituentPts.size(); i++)
    {
        TotalPt += ConstituentPts[i];
        if (ConstituentPids[i] == 22)
        {
            CurrentPt = ConstituentPts[i];
            if (CurrentPt > MaxPhotonPt)
                MaxPhotonPt = CurrentPt;
        }
    }
    double PhotonPtFrac = MaxPhotonPt / TotalPt;
    return (PhotonPtFrac > 0.8);
}

double getAngularity(double kappa, double beta, vector<double> &constituentPt, double jetPt, vector<double> &constituentDr, double R = 0.4)
{
    double lambda = 0;
    for (int i = 0; i < constituentPt.size(); i++)
    {
        double z = constituentPt[i] / jetPt;
        double theta = constituentDr[i] / R;
        lambda += (pow(z, kappa) * pow(theta, beta));
    }
    return lambda;
}

vector<double> getObservableVector(vector<double> &constituentPt, double jetPt, vector<double> &constituentDr)
{
    double multiplicity = getAngularity(0, 0, constituentPt, jetPt, constituentDr, 0.4); // including the 100 as a scaling factor
    double pTD = getAngularity(2, 0, constituentPt, jetPt, constituentDr, 0.4);
    double LHA = getAngularity(1, 0.5, constituentPt, jetPt, constituentDr, 0.4);
    double width = getAngularity(1, 1, constituentPt, jetPt, constituentDr, 0.4);
    double mass = getAngularity(1, 2, constituentPt, jetPt, constituentDr, 0.4);

    // getting vector representing multiplicity
    double v[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10, 10-20, 20+
    for (int i = 0; i < constituentPt.size(); i++)
    {
        double trackPt = constituentPt[i];
        if (trackPt >= 20)
        {
            v[11] += 1;
        }
        else if (trackPt >= 10)
        {
            v[10] += 1;
        }
        else if (trackPt >= 9)
        {
            v[9] += 1;
        }
        else if (trackPt >= 8)
        {
            v[8] += 1;
        }
        else if (trackPt >= 7)
        {
            v[7] += 1;
        }
        else if (trackPt >= 6)
        {
            v[6] += 1;
        }
        else if (trackPt >= 5)
        {
            v[5] += 1;
        }
        else if (trackPt >= 4)
        {
            v[4] += 1;
        }
        else if (trackPt >= 3)
        {
            v[3] += 1;
        }
        else if (trackPt >= 2)
        {
            v[2] += 1;
        }
        else if (trackPt >= 1)
        {
            v[1] += 1;
        }
        else if (trackPt >= 0)
        {
            v[0] += 1;
        }
    }

    vector<double> returnVector = {multiplicity, pTD, LHA, width, mass};
    for (int i = 0; i < 12; i++)
    {
        returnVector.push_back(v[i]);
    }
    return returnVector;
}

vector<double> getDRVector(double eta, vector<double> &constituentEta, double phi, vector<double> &constituentPhi)
{
    vector<double> drVector;
    double dr;
    for (int i = 0; i < constituentEta.size(); i++)
    {
        // cout << "check 4" << "\n";
        dr = getR(constituentEta[i], eta, constituentPhi[i], phi);
        drVector.push_back(dr);
    }
    return drVector;
}

void writeVecToCsv(vector<double> &obsVector, ofstream &outFile1, ofstream &outFile2, int matrixElem, bool qg = true)
{
    double quarkList[] = {1, 2, 3, 4, 5, 6};

    if (qg)
    {
        if (matrixElem == 21)
        {
            for (const auto &e : obsVector)
                outFile2 << e << " ";
            outFile2 << "\n";
        }
        else if (find(begin(quarkList), end(quarkList), matrixElem) != end(quarkList))
        {
            for (const auto &e : obsVector)
                outFile1 << e << " ";
            outFile1 << "\n";
        }
    }
    else
    {
        // let's make this so we only print to outFile1
        for (const auto &e : obsVector)
            outFile1 << e << " ";
        outFile1 << "\n";
    }
}

int main(int argc, char *argv[])
{
    SetThesisStyle();
    CommandLine CL(argc, argv);

    vector<string> inputFileFolders = CL.GetStringVector("input");
    string OutputBase = CL.Get("output");

    double MinPT = CL.GetDouble("minpt", 0);
    double MaxPT = CL.GetDouble("maxpt", 500);
    // bool qg = (bool) CL.GetInt("qg", 1);
    // cout << qg << endl;
    // double MaxTrackPT = CL.GetInt("maxtrkpt", 100); // this needs to be int like 200 or something

    // string trackPtCuts = "0:1:4:10";
    string outputLabel = OutputBase;

    ofstream outFileQuark(outputLabel + "_quark.txt", ios_base::app);
    ofstream outFileGluon(outputLabel + "_gluon.txt", ios_base::app);
    ofstream outFileDijet(outputLabel + "_dijet.txt", ios_base::app);
    ofstream outFilePhotonJet(outputLabel + "_photonjet.txt", ios_base::app);
    ofstream outFileDijetQuark(outputLabel + "_dijet_quark_truth.txt", ios_base::app);
    ofstream outFileDijetGluon(outputLabel + "_dijet_gluon_truth.txt", ios_base::app);
    ofstream outFilePhotonJetQuark(outputLabel + "_photonjet_quark_truth.txt", ios_base::app);
    ofstream outFilePhotonJetGluon(outputLabel + "_photonjet_gluon_truth.txt", ios_base::app);

    for (string folder : inputFileFolders)
    {
        bool dijet = isDijet(folder);

        struct dirent *de; // to get files in path
        DIR *dr = opendir(folder.c_str());
        if (dr == NULL) // opendir returns NULL if couldn't open directory
        {
            cout << "Could not open directory: " << folder << endl;
            return 0;
        }

        int fileCount = 0; // REMOVE
        while ((de = readdir(dr)) != NULL)
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
            vector<double> *SignalJetMatrixElemPt = nullptr;
            vector<double> *SignalJetMatrixElemPx = nullptr;
            vector<double> *SignalJetMatrixElemPy = nullptr;
            vector<double> *SignalJetMatrixElemMass = nullptr;
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
            Tree->SetBranchAddress("SignalJetMatrixElemPt", &SignalJetMatrixElemPt);
            Tree->SetBranchAddress("SignalJetMatrixElemPx", &SignalJetMatrixElemPx);
            Tree->SetBranchAddress("SignalJetMatrixElemPy", &SignalJetMatrixElemPy);
            Tree->SetBranchAddress("SignalJetMatrixElemMass", &SignalJetMatrixElemMass);
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

                if (MaxJetPtIndices.first >= 0 && MaxJetPts.first >= MinPT && MaxJetPts.first <= MaxPT)
                {
                    if (hasLargePhotonPtFrac((*SignalJetConstituentPt)[MaxJetPtIndices.first], (*SignalJetConstituentPid)[MaxJetPtIndices.first]))
                        continue;
                    if (fabs((*SignalJetEta)[MaxJetPtIndices.first]) < 1)
                    {
                        int matElemIdx;
                        double matElemDr;
                        tie(matElemIdx, matElemDr) = getMatElemIndexDr(dijet, SignalJetMatrixElem, SignalJetMatrixElemPhi, SignalJetMatrixElemEta, (*SignalJetPhi)[MaxJetPtIndices.first], (*SignalJetEta)[MaxJetPtIndices.first]);
                        if (matElemIdx == -1)
                            continue;
                        int matElem = (*SignalJetMatrixElem)[matElemIdx];

                        // double etaMatrix, double etaJet, double phiMatrix, double phiJet
                        vector<double> dr = getDRVector((*SignalJetEta)[MaxJetPtIndices.first], (*SignalJetConstituentEta)[MaxJetPtIndices.first], (*SignalJetPhi)[MaxJetPtIndices.first], (*SignalJetConstituentPhi)[MaxJetPtIndices.first]);
                        vector<double> obsVector = getObservableVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], (*SignalJetPt)[MaxJetPtIndices.first], dr);

                        // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], outFileQuark, outFileGluon, matElem, true);
                        writeVecToCsv(obsVector, outFileQuark, outFileGluon, matElem, true);
                        if (!dijet)
                        {
                            // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], outFilePhotonJet, outFileDijet, matElem, false);
                            writeVecToCsv(obsVector, outFilePhotonJet, outFileDijet, matElem, false);
                            if (matElemDr <= 0.4)
                            {
                                // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], outFilePhotonJetQuark, outFilePhotonJetGluon, matElem, true);
                                writeVecToCsv(obsVector, outFilePhotonJetQuark, outFilePhotonJetGluon, matElem, true);
                            }
                        }
                        else
                        {
                            // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], outFileDijet, outFilePhotonJet, matElem, false);
                            writeVecToCsv(obsVector, outFileDijet, outFilePhotonJet, matElem, false);
                            if (matElemDr <= 0.4)
                            {
                                // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.first], outFileDijetQuark, outFileDijetGluon, matElem, true);
                                writeVecToCsv(obsVector, outFileDijetQuark, outFileDijetGluon, matElem, true);
                            }
                        }
                    }
                }
                if (dijet && MaxJetPtIndices.second >= 0 && MaxJetPts.second >= MinPT && MaxJetPts.second <= MaxPT)
                {
                    if (hasLargePhotonPtFrac((*SignalJetConstituentPt)[MaxJetPtIndices.second], (*SignalJetConstituentPid)[MaxJetPtIndices.second]))
                        continue;
                    if (fabs((*SignalJetEta)[MaxJetPtIndices.second]) < 1)
                    {
                        int matElemIdx;
                        double matElemDr;
                        tie(matElemIdx, matElemDr) = getMatElemIndexDr(dijet, SignalJetMatrixElem, SignalJetMatrixElemPhi, SignalJetMatrixElemEta, (*SignalJetPhi)[MaxJetPtIndices.second], (*SignalJetEta)[MaxJetPtIndices.second]);
                        if (matElemIdx == -1)
                            continue;
                        int matElem = (*SignalJetMatrixElem)[matElemIdx];

                        // double etaMatrix, double etaJet, double phiMatrix, double phiJet
                        vector<double> dr = getDRVector((*SignalJetEta)[MaxJetPtIndices.second], (*SignalJetConstituentEta)[MaxJetPtIndices.second], (*SignalJetPhi)[MaxJetPtIndices.second], (*SignalJetConstituentPhi)[MaxJetPtIndices.second]);
                        vector<double> obsVector = getObservableVector((*SignalJetConstituentPt)[MaxJetPtIndices.second], (*SignalJetPt)[MaxJetPtIndices.second], dr);

                        writeVecToCsv(obsVector, outFileQuark, outFileGluon, matElem, true);
                        writeVecToCsv(obsVector, outFileDijet, outFilePhotonJet, matElem, false);
                        // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.second], outFileQuark, outFileGluon, matElem, true);
                        // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.second], outFileDijet, outFilePhotonJet, matElem, false);
                        if (matElemDr <= 0.4)
                        {
                            // getLDMultiplicityVector((*SignalJetConstituentPt)[MaxJetPtIndices.second], outFileDijetQuark, outFileDijetGluon, matElem, true);
                            writeVecToCsv(obsVector, outFileDijetQuark, outFileDijetGluon, matElem, true);
                        }
                    }
                }
            }
        }
    }

    return 0;
}
