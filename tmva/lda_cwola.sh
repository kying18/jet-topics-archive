# input="pp" # pp or PbPbWide
pthat=80
# dijet_input = "/data/kying/final80/pp80"
# photonjet_input = "/data/kying/final80/ppPhotonJet80"
root_input_trkpt_cuts="0:1:4:10"
trkpt_cuts="0:1:4:10" # TODO LDA bins change here
# minpt=100
# maxpt=120
tmva="mlp"

count=0

declare -a suffixes=("_dijet" "_photonjet" "_dijet_quark_truth" "_dijet_gluon_truth" "_photonjet_quark_truth" "_photonjet_gluon_truth")

>nohup.out

g++ getLDMultiplicity.cpp $(root-config --cflags --libs) -O2 -o "getLDMultiplicity.exe"
g++ getTopicModelingCsv.cpp $(root-config --cflags --libs) -O2 -o "getTopicModelingCsv.exe"

declare -a inputs=("pp" "PbPbWide")
for minpt in `seq 80 20 100`
do 
    (
        for input in ${inputs[@]}
        do
            
            maxpt=$(($minpt + 20))
            
            # echo $input $minpt $maxpt
            # sleep 5

            ########################################################################
            # # we should already have these files!
            # # first thing to do is to run these to get the text files
            ########################################################################
            # for suffix in ${suffixes[@]}
            # do
            #     > /data/kying/lda/LDMultTxt/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}_cwola.txt
            # done

            # if [ "$input" == "pp" ]; then
            #     input_files="/data/kying/final80/ppPhotonJet80,/data/kying/final80/pp80"
            # else
            #     input_files="/data/kying/final80/PbPbWidePhotonJet80_0_10,/data/kying/final80/PbPbWide80_0_10"
            # fi

            # ./getLDMultiplicity.exe -input $input_files -minpt $minpt -maxpt $maxpt -output /data/kying/lda/LDMultTxt/${input}_pt${minpt}${maxpt}
            ########################################################################


            ########################################################################
            # # should already have these
            # # then turn txt into root tree
            ########################################################################
            # for suffix in ${suffixes[@]}
            # do
            #     > /data/kying/lda/LDMultRoot/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.root
            #     cat /data/kying/lda/LDMultTxt/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.txt | ./TextToTree /data/kying/lda/LDMultRoot/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.root 4 "pt01:pt14:pt410:pt10inf"
            # done

            # echo "Got root trees"
            ########################################################################

            # # remove the txt files bc they take up way too much space
            # for suffix in ${suffixes[@]}
            # do
            #     rm -rf ./LDMultTxt/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.txt
            # done

            # now run the LDA analysis on the quark/gluon trees
            # if [ "$input" == "pp" ]; then
            >/data/kying/lda/LDResults/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola.root
            root -b -l -q 'RunLD.C("'/data/kying/lda/LDMultRoot/${input}_pt${minpt}${maxpt}_trkpt${root_input_trkpt_cuts}'","'/data/kying/lda/LDResults/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}'_cwola.root","'${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola_'","dijet","photonjet")'
            echo "Ran TMVA..."
            # fi

            # root -l 'RunLD.C("LDMultRoot/pp_pt100120_trkpt0:1:4:10","LDResults/pp_pt100120_trkpt0:1:4:10.root")'

            # and then run the LDA application on dijet/photon jet inputs
            for suffix in ${suffixes[@]}
            do
                >/data/kying/lda/LDApplication/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}_cwola.root
                root -b -l -q 'RunLDApplication.C("'/data/kying/lda/LDMultRoot/${input}_pt${minpt}${maxpt}_trkpt${root_input_trkpt_cuts}${suffix}.root'","'/data/kying/lda/LDApplication/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}_cwola.root'","'${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola_'")'
            done
            # root -l 'RunLDApplication.C("LDMultRoot/pp_pt80100_trkpt0:1:4:10_photonjet.root","LDApplication/pp_pt80100_trkpt0:1:4:10_photonjet.root")
            echo "Ran TMVA application"

            # # remove the root files bc they take up way too much space
            # for suffix in ${suffixes[@]}
            # do
            #     rm -rf LDMultRoot/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.root
            # done

            # and get the final output histogram
            label_prefix="pbpb80_0_10_wide"
            if [ "$input" == "pp" ]; then
                >/data/kying/lda/LDMultCsvs/${tmva}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola.csv
                label_prefix="pp80"
            fi

            declare hist_suffixes=("" "_quark_truth" "_gluon_truth")
            for suffix in "${hist_suffixes[@]}"
            do
                echo $suffix
                ./getTopicModelingCsv.exe -input /data/kying/lda/LDApplication/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_photonjet${suffix}_cwola.root -labels ${label_prefix}_photonjet${suffix} -output /data/kying/lda/LDMultCsvs/${tmva}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola.csv
                ./getTopicModelingCsv.exe -input /data/kying/lda/LDApplication/${tmva}_${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_dijet${suffix}_cwola.root -labels ${label_prefix}${suffix} -output /data/kying/lda/LDMultCsvs/${tmva}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola.csv
            done

            # # remove the root files bc they take up way too much space
            # for suffix in "${suffixes[@]}"
            # do
            #     rm -rf LDApplication/${input}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}${suffix}.root
            # done
        done
        cp /data/kying/lda/LDMultCsvs/${tmva}_pt${minpt}${maxpt}_trkpt${trkpt_cuts}_cwola.csv ~/jet-topics/data
    ) &
    let count+=1
    [[ $((count%5)) -eq 0 ]] && wait
done