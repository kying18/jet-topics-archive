# Jet Topics Additional Scripts

In this folder, you'll find the scripts that generated the data used for the topic modeling project (https://arxiv.org/abs/2204.00641). This includes the following:


## `runLundPlane` script
Extracts jet information given pu14 files (modified from the original in the `JetToyHI` repo to include matrix elements and certain particle information)


## `observables` folder:
Files for extracting the substructure observables (jet shape, jet fragmentation, jet mass, and jet splitting fraction). This folder is pretty self explanatory, there are 4 cpp files representing the 4 observables. To get the observable values in each, just compile and run the files.

The specific commands are:
```
g++ jetX.cpp $(root-config --cflags --libs) -O2 -o "jetX.exe"

./jetX.exe -input "data_folder_path/pp80,data_folder_path/ppPhotonJet80" -minpt 80 -maxpt 100
```

These will output values to a csv file (look in the files to change the output folder name.. sorry, didn't implement that part as an input option)

## `ml` folder
Files for extracting the features of the jets to input into the machine learning algorithms (LHA, pTD, constituent multiplicity in pT bins, etc.)

The flow of these files can be found in `lda_{x}.sh` or `tmvaJetObs.sh` (yeah sorry, going forward, these bash files can be cleaned up and merged):

1) First get txt files. Each row represents a jet, each column represents a data column. The data column might be the multiplicity of each pT bin (`getLDMultiplicity.cpp`) or it might be an observable (`getJetObs.cpp`).
2) Then, we turn the txt files into root trees thanks to `TextToTree.cpp`
3) Now, run the TMVA algorithm. The first one I ran was LDA (linear discriminant analysis), so the naming stuck. However, `RunLD.C` runs the TMVA algo on your trees and gets a classifier, which saves to the `dataset/weights` folder.
4) In order to apply the algo to the entire dataset, run `RunLDApplication.C`.
5) Finally, we want the csv file representing the histograms so we can input this into the jet topic modeling machinery. That is at `getTopicModelingCsv.cpp`.

NOTE: if you prefer Python for machine learning (which I do), you can also directly use the txt files that get generated and use those to train a model elsewhere. I do this using a jupyter notebook, which I uploaded as `mlnotebook.ipynb` as an example

I started an overleaf file containing preliminary findings about the ML portion of the project. View it here: https://www.overleaf.com/read/ywtqwtnsnxst

If you wish to extend results, please reach out to me at me{at}kylieying.com or contact Professor Yen-Jie Lee to contact me.

Good luck!!