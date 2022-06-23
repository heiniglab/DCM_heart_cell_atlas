#!/bin/bash

#-------------
#   Solo quantification
#   This script was run from the cluster VM (with conda envir solo)
#-------------

# runs on 5k cells ~35mins (Tesla-T4)

conda activate solo
path_origin=/home/elindbe/fast/AG_Huebner/huebner3/ANALYSES/20190926_gp_BadOyenhausen/scanpy/Global_touse/individual_h5ads_20200430
path_model_json_file=/home/elindbe/fast/AG_Huebner/huebner3/ANALYSES/20190926_gp_BadOyenhausen/scanpy/scripts/Global/solo_Step2.2/model_json_file

cd ${path_origin}

while read sample; do
    echo -e "\n###\nProcessing sample ${sample}\n###\n"
    
    new_name=`echo ${sample} | sed -e 's/\(.*\)\//\1/g'`
    
    mkdir $path_origin/${new_name}/solo_out    
    cd $path_origin/${new_name}
    solo -o $path_origin/${new_name}/solo_out $path_model_json_file $path_origin/${new_name}/${new_name}.h5ad
    echo -e "\n---\nFinished with ${sample}\n---\n"

done </home/elindbe/fast/AG_Huebner/huebner3/ANALYSES/20190926_gp_BadOyenhausen/scanpy/scripts/Global/solo_Step2.2/all_samples.txt
