#include <fstream>
#include <iterator>
#include <cmath>
#include <dirent.h>
#include <math.h>
#include <set>
#include <vector>
#include <cassert>
using namespace std;

bool oppositeHemisphere(double particlePhi, double jetPhi) {
    if(particlePhi >= 0) {
        double dPhi = particlePhi - jetPhi;
        if(dPhi < -M_PI)   dPhi = dPhi + 2 * M_PI;
        if(dPhi > +M_PI)   dPhi = dPhi - 2 * M_PI;

        if(dPhi < - M_PI / 2 || dPhi > M_PI / 2) {
            return true;
        }
        return false;
    }

    return false;
}

bool isDijet(string folderName) {
    if (folderName.find("PhotonJet") != string::npos) return false;
    if (folderName.find("photonjet") != string::npos) return false;
    return true;
}

bool ispp(string folderName) {
    if (folderName.find("pp") != string::npos) return true;
    return false;
}

bool hasLargePhotonPtFrac(vector<double>ConstituentPts, vector<double>ConstituentPids) {
    assert (ConstituentPts.size() == ConstituentPids.size());
    // get the largest photon pt fraction and cut if larger than 0.8
    double MaxPhotonPt = 0;
    double TotalPt = 0;
    double CurrentPt;

    for (int i=0; i < ConstituentPts.size(); i++) {
        TotalPt += ConstituentPts[i];
        if (ConstituentPids[i] == 22) {
            CurrentPt = ConstituentPts[i];
            if (CurrentPt > MaxPhotonPt) MaxPhotonPt = CurrentPt;
        }
    }
    double PhotonPtFrac = MaxPhotonPt / TotalPt;
    return (PhotonPtFrac > 0.8);
}

pair<int, int> getMaxPtJetIndices(int NJet, bool dijet, double LeadingPhotonPhi, vector<double>*SignalJetPhi, vector<double>*SignalJetPt, vector<vector<double>>*SignalJetConstituentPt, vector<vector<double>>*SignalJetConstituentPid) {
    pair<int, int> MaxJetPtIndices (-1, -1);
    pair<double, double> MaxJetPts (0, 0);

    for (int iJ = 0; iJ < NJet; iJ++) {
        // if photon jet and jet is not in opposite hemisphere, then we skip this jet
        if (!dijet && !oppositeHemisphere(LeadingPhotonPhi, (*SignalJetPhi)[iJ])) {
            // cout << "leading photon phi: " << LeadingPhotonPhi << ", signal jet phi: " << (*SignalJetPhi)[iJ] << endl;
            continue; 
        }
        // if (hasLargePhotonPtFrac((*SignalJetConstituentPt)[iJ], (*SignalJetConstituentPid)[iJ])) continue;

        // cout << "Signal Jet Pt: " << (*SignalJetPt)[iJ] << endl;

        if ((*SignalJetPt)[iJ] > MaxJetPts.first) {
            // move current first to second, insert current jet as first
            MaxJetPts.second = MaxJetPts.first;
            MaxJetPtIndices.second = MaxJetPtIndices.first;
            MaxJetPts.first = (*SignalJetPt)[iJ];
            MaxJetPtIndices.first = iJ;
        } else if ((*SignalJetPt)[iJ] > MaxJetPts.second) {
            MaxJetPts.second = (*SignalJetPt)[iJ];
            MaxJetPtIndices.second = iJ;
        }
    }
    return MaxJetPtIndices;
}

double getR(double etaMatrix, double etaJet, double phiMatrix, double phiJet) {
    // this is the reconstructed mat elem's radial distance from jet axis (defined by etaJet and phiJet)
    double deta = etaMatrix - etaJet;
    double dphi = phiMatrix - phiJet;

    // dphi should be between -pi and pi
    if (dphi < -M_PI) {
        dphi = 2*M_PI - abs(dphi);  // if phi is < -pi, then 2pi - value of dphi is new
    } else if (dphi > M_PI) {
        dphi = - (2*M_PI - dphi); // if phi > pi, then 2pi - dphi is new value, but need to negate
    }

    assert(dphi >= -M_PI && dphi <= M_PI);

    return sqrt(pow(deta,2) + pow(dphi, 2));
}

pair<int, double> getMatElemIndexDr(bool dijet, vector<int>*SignalJetMatrixElem, vector<double>*SignalJetMatrixElemPhi, vector<double>*SignalJetMatrixElemEta, double SignalJetPhi, double SignalJetEta) {
    // returns the int rep index of the matrix elem in the vector and the double rep DR of the jet to the matrix elem
    // all of these should be size 2
    assert (SignalJetMatrixElem->size() == 2 && SignalJetMatrixElemPhi->size() == 2 && SignalJetMatrixElemEta->size() == 2);
    int index;
    double dr;

    if (!dijet) {
        // if photon jet, then return the matrix element that is not the photon
        index = ((*SignalJetMatrixElem)[0] != 22? 0 : 1);
        dr = getR((*SignalJetMatrixElemEta)[index], SignalJetEta, (*SignalJetMatrixElemPhi)[index], SignalJetPhi);
        return make_pair(index, dr);
    }
    // otherwise if dijet, compute the r value and select the smaller
    double index_0_r = getR((*SignalJetMatrixElemEta)[0], SignalJetEta, (*SignalJetMatrixElemPhi)[0], SignalJetPhi);
    double index_1_r = getR((*SignalJetMatrixElemEta)[1], SignalJetEta, (*SignalJetMatrixElemPhi)[1], SignalJetPhi);

    index = (index_0_r <= index_1_r ? 0 : 1);
    dr = (index_0_r <= index_1_r ? index_0_r : index_1_r);
    return make_pair(index, dr);

}