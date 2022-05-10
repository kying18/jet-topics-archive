# input="pp" # pp or PbPbWide
pthat=80
tmva="lda"
folder_prefix="TMVAObsSqrt"

count=0

declare -a suffixes=("_quark" "_gluon" "_dijet" "_photonjet" "_dijet_quark_truth" "_dijet_gluon_truth" "_photonjet_quark_truth" "_photonjet_gluon_truth")

>nohup.out

g++ getJetObs.cpp $(root-config --cflags --libs) -O2 -o "getJetObs.exe"
g++ getTopicModelingCsv.cpp $(root-config --cflags --libs) -O2 -o "getTopicModelingCsv.exe"

mkdir -p /data/kying/lda/${folder_prefix}Txt
mkdir -p /data/kying/lda/${folder_prefix}Root
mkdir -p /data/kying/lda/${folder_prefix}Csvs
mkdir -p /data/kying/lda/${folder_prefix}Results
mkdir -p /data/kying/lda/${folder_prefix}Application
mkdir -p ~/jet-topics/data/${folder_prefix}

declare -a inputs=("pp" "PbPbWide")
for minpt in `seq 80 20 100`
do 
    (
        for input in ${inputs[@]}
        do
            
            maxpt=$(($minpt + 20))
            
            # echo $input $minpt $maxpt
            # sleep 5

            ##########################################################################################################################################
            # already did this 
            ##########################################################################################################################################
            # first thing to do is to run these to get the text files
            for suffix in ${suffixes[@]}
            do
                > /data/kying/lda/${folder_prefix}Txt/${input}_pt${minpt}${maxpt}${suffix}.txt
            done

            if [ "$input" == "pp" ]; then
                input_files="/data/kying/final80/ppPhotonJet80,/data/kying/final80/pp80"
            else
                input_files="/data/kying/final80/PbPbWidePhotonJet80_0_10,/data/kying/final80/PbPbWide80_0_10"
            fi


            ./getJetObs.exe -input $input_files -minpt $minpt -maxpt $maxpt -output /data/kying/lda/${folder_prefix}Txt/${input}_pt${minpt}${maxpt}

            # then turn txt into root tree
            for suffix in ${suffixes[@]}
            do
                > /data/kying/lda/${folder_prefix}Root/${input}_pt${minpt}${maxpt}${suffix}.root
                cat /data/kying/lda/${folder_prefix}Txt/${input}_pt${minpt}${maxpt}${suffix}.txt | ./TextToTree /data/kying/lda/${folder_prefix}Root/${input}_pt${minpt}${maxpt}${suffix}.root 5 "multiplicity:pTD:LHA:width:mass"
            done

            echo "Got root trees"

            # # remove the txt files bc they take up way too much space
            # for suffix in ${suffixes[@]}
            # do
            #     rm -rf ./LDMultTxt/${input}_pt${minpt}${maxpt}${suffix}.txt
            # done
            ##########################################################################################################################################

            # now run the LDA analysis on the quark/gluon trees but only if pp
            if [ "$input" == "pp" ]; then
                >/data/kying/lda/${folder_prefix}Results/${tmva}_${input}_pt${minpt}${maxpt}.root
                root -b -l -q 'RunLD.C("'/data/kying/lda/${folder_prefix}Root/${input}_pt${minpt}${maxpt}'","'/data/kying/lda/${folder_prefix}Root/${tmva}_${input}_pt${minpt}${maxpt}'.root","'${input}_pt${minpt}${maxpt}_'","gluon","quark")'
                echo "Ran LDA..."
            fi

            # root -l 'RunLD.C("LDMultRoot/pp_pt100120_trkpt0:1:4:10","LDResults/pp_pt100120_trkpt0:1:4:10.root")'

            # and then run the LDA application on dijet/photon jet inputs
            for suffix in ${suffixes[@]}
            do
                >/data/kying/lda/${folder_prefix}Application/${tmva}_${input}_pt${minpt}${maxpt}${suffix}.root
                root -b -l -q 'RunLDApplication.C("'/data/kying/lda/${folder_prefix}Root/${input}_pt${minpt}${maxpt}${suffix}.root'","'/data/kying/lda/${folder_prefix}Application/${tmva}_${input}_pt${minpt}${maxpt}${suffix}.root'","'pp_pt${minpt}${maxpt}_'")'
            done
            # root -l 'RunLDApplication.C("/data/kying/lda/LDMultRoot/pp_pt80100_trkpt2:4:10_photonjet.root","/data/kying/lda/LDApplication/bdt_pp_pt80100_trkpt2:4:10_photonjet.root","pp_pt80100_trkpt2:4:10_")'
            # root -l 'RunLDApplication.C("LDMultRoot/pp_pt80100_trkpt0:1:4:10_photonjet.root","LDApplication/pp_pt80100_trkpt0:1:4:10_photonjet.root")
            echo "Ran LDA application"

            ##########################################################################################################################################

            # # remove the root files bc they take up way too much space
            # for suffix in ${suffixes[@]}
            # do
            #     rm -rf LDMultRoot/${input}_pt${minpt}${maxpt}${suffix}.root
            # done

            # and get the final output histogram
            label_prefix="pbpb80_0_10_wide"
            if [ "$input" == "pp" ]; then
                >/data/kying/lda/${folder_prefix}Csvs/${tmva}_pt${minpt}${maxpt}.csv
                label_prefix="pp80"
            fi

            declare hist_suffixes=("" "_quark_truth" "_gluon_truth")
            for suffix in "${hist_suffixes[@]}"
            do
                echo $suffix
                ./getTopicModelingCsv.exe -input /data/kying/lda/${folder_prefix}Application/${tmva}_${input}_pt${minpt}${maxpt}_photonjet${suffix}.root -labels ${label_prefix}_photonjet${suffix} -output /data/kying/lda/${folder_prefix}Csvs/${tmva}_pt${minpt}${maxpt}.csv
                ./getTopicModelingCsv.exe -input /data/kying/lda/${folder_prefix}Application/${tmva}_${input}_pt${minpt}${maxpt}_dijet${suffix}.root -labels ${label_prefix}${suffix} -output /data/kying/lda/${folder_prefix}Csvs/${tmva}_pt${minpt}${maxpt}.csv
            done

            # # remove the root files bc they take up way too much space
            # for suffix in "${suffixes[@]}"
            # do
            #     rm -rf LDApplication/${input}_pt${minpt}${maxpt}${suffix}.root
            # done


        done
        cp /data/kying/lda/${folder_prefix}Csvs/${tmva}_pt${minpt}${maxpt}.csv ~/jet-topics/data/${folder_prefix}
    ) &
    let count+=1
    [[ $((count%5)) -eq 0 ]] && wait
done