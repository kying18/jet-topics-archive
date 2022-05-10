#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int RunLDApplication(string input, string output, string label)
{
   TMVA::Tools::Instance();

   TFile dataFile(input.c_str());
   TTree *dataTree = (TTree *)dataFile.Get("Tree");

   // TString outfileName = "LDResults/pp_pt100120_trkpt0,1,4,10_dijet.root";
   TFile *outputFile = TFile::Open(output.c_str(), "RECREATE");

   TTree *tree = new TTree("tree", "treelibrated tree");

   TMVA::Tools::Instance();
   TMVA::Reader *reader = new TMVA::Reader("V:Color:!Silent");

   Float_t var[4]; // TODO LDA change here

   // TODO LDA change here
   // reader->AddVariable("multiplicity",&var[0]);
   // reader->AddVariable("pTD",&var[1]);
   // reader->AddVariable("LHA",&var[2]);
   // reader->AddVariable("width",&var[3]);
   // reader->AddVariable("mass",&var[4]);

   reader->AddVariable("pt01", &var[0]);
   reader->AddVariable("pt14", &var[1]);
   reader->AddVariable("pt410", &var[2]);
   reader->AddVariable("pt10inf", &var[3]);

   // reader->BookMVA("LDA method", ("dataset/weights/" + label + "TMVAClassification_LD.weights.xml").c_str());
   // reader->BookMVA("BDT method", ("dataset/weights/" + label + "TMVAClassification_BDT.weights.xml").c_str());
   reader->BookMVA("MLP method", ("dataset/weights/" + label + "TMVAClassification_MLP.weights.xml").c_str());

   Double_t userVar[4]; // TODO LDA change here

   // TODO LDA change here
   // dataTree->SetBranchAddress("multiplicity", &userVar[0]);
   // dataTree->SetBranchAddress("pTD", &userVar[1]);
   // dataTree->SetBranchAddress("LHA", &userVar[2]);
   // dataTree->SetBranchAddress("width", &userVar[3]);
   // dataTree->SetBranchAddress("mass", &userVar[4]);

   dataTree->SetBranchAddress("pt01", &userVar[0]);
   dataTree->SetBranchAddress("pt14", &userVar[1]);
   dataTree->SetBranchAddress("pt410", &userVar[2]);
   dataTree->SetBranchAddress("pt10inf", &userVar[3]);

   // Float_t LDA_response;
   // tree->Branch("LDA_response", &LDA_response);
   // Float_t BDT_response;
   // tree->Branch("BDT_response", &BDT_response);
   Float_t MLP_response;
   tree->Branch("MLP_response", &MLP_response);

   for (Long64_t ievt = 0; ievt < dataTree->GetEntries(); ievt++)
   {

      if (ievt % 100000 == 0)
         std::cout << "--- ... Processing event: " << ievt << std::endl;

      dataTree->GetEntry(ievt);

      // TODO LDA change here
      var[0] = userVar[0];
      var[1] = userVar[1];
      var[2] = userVar[2];
      var[3] = userVar[3];
      // var[4]=userVar[4];

      // LDA_response = reader->EvaluateMVA("LDA method");
      // BDT_response = reader->EvaluateMVA("BDT method");
      MLP_response = reader->EvaluateMVA("MLP method");

      tree->Fill();
   }

   tree->Write();

   return 0;
}

int main(int argc, char **argv)
{
   return RunLDApplication(argv[0], argv[1], argv[2]);
}