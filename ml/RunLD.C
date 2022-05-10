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

int RunLD(string filePath, string outfileName, string label, string signal, string background)
{
   TMVA::Tools::Instance();

   TFile Signal((filePath + "_" + signal + ".root").c_str());
   TFile Background((filePath + "_" + background + ".root").c_str());

   TTree *TreeSignal = (TTree *)Signal.Get("Tree");
   TTree *TreeBackground = (TTree *)Background.Get("Tree");

   // TString outfileName = "LDResults/pp_pt80100_trkpt0,1,4,10.root";
   TFile *outputFile = TFile::Open(outfileName.c_str(), "RECREATE");

   TMVA::Factory *factory = new TMVA::Factory((label + "TMVAClassification").c_str(), outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   // TODO LDA bins change here
   // dataloader->AddVariable("multiplicity", 'D');
   // dataloader->AddVariable("pTD", 'D');
   // dataloader->AddVariable("LHA", 'D');
   // dataloader->AddVariable("width", 'D');
   // dataloader->AddVariable("mass", 'D');
   dataloader->AddVariable("pt01", 'D');
   dataloader->AddVariable("pt14", 'D');
   dataloader->AddVariable("pt410", 'D');
   dataloader->AddVariable("pt10inf", 'D');

   dataloader->AddSignalTree(TreeSignal);
   dataloader->AddBackgroundTree(TreeBackground);

   dataloader->PrepareTrainingAndTestTree("", "",
                                          "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

   // factory->BookMethod(dataloader, TMVA::Types::kLD, "LD",
   //                     "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
   // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
   //                     "BoostType=AdaBoost:SeparationType=GiniIndex:Ntrees=200");
   factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP",
                       "!H:!V:NeuronType=tanh:NCycles=100:HiddenLayers=N-1:TrainingMethod=BP");

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   outputFile->Close();

   delete factory;
   delete dataloader;

   if (!gROOT->IsBatch())
      TMVA::TMVAGui(outfileName.c_str());

   return 0;
}

int main(int argc, char **argv)
{
   return RunLD(argv[0], argv[1], argv[2], argv[3], argv[4]);
}