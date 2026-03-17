#!/bin/bash

P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
cd ${P_DIR}

old_participants=("2E_001" "2E_002" "2E_003" "2E_004" "2E_005" "2E_006" "2E_008" "2E_009" "2E_010" "2E_011" "2E_012" "2E_013" "2E_015" "2E_016" "2E_017" "2E_018" "2E_019" "2E_020" "2E_021" "2E_022" "2E_024" "2E_025" "2E_027" "2E_028" "2E_029" "2E_030")
all_participants=("2E_022" "2E_029" "2E_023" "2E_031" "2E_032" "2E_033" "2E_034" "2E_035" "2E_038" "2E_039" "2E_040" "2E_041" "2E_042" "2E_043" "2E_044" "2E_045" "2E_047" "2E_048" "2E_049" "2E_050" "2E_051" "2E_052" "2E_053" "2E_054" "2E_055" "2E_056" "2E_057" "2E_066" "2E_070" "2E_075" "2E_084" "2E_085" "2E_090" "2E_095" "2E_096" "2E_097" "2E_098" "2E_099" "2E_100" "2E_102" "2E_104" "2E_105" "2E_106" "2E_108" "2E_109" "2E_112" "2E_115" "2E_124" "2E_126" "2E_118" "2E_131" "2E_133" "2E_134")

fmri_participants=("2E_001" "2E_002" "2E_003" "2E_004" "2E_005" "2E_006" "2E_008" "2E_009" "2E_010" "2E_011" "2E_012" "2E_013" "2E_016" "2E_017" "2E_018" "2E_019" "2E_020" "2E_022" "2E_024" "2E_025" "2E_027" "2E_028" "2E_029" "2E_030" "2E_023" "2E_031" "2E_032" "2E_033" "2E_034" "2E_035" "2E_038" "2E_039" "2E_040" "2E_041" "2E_042" "2E_043" "2E_044" "2E_045" "2E_048" "2E_050" "2E_051" "2E_052" "2E_053" "2E_054" "2E_055" "2E_056" "2E_057" "2E_066" "2E_070" "2E_075" "2E_084" "2E_085" "2E_090" "2E_095" "2E_096" "2E_097" "2E_098" "2E_099" "2E_100" "2E_102" "2E_104" "2E_105" "2E_106" "2E_108" "2E_109" "2E_112" "2E_115" "2E_124" "2E_126" "2E_118" "2E_131" "2E_133" "2E_134")

ALLLL_participants=("2E_001" "2E_002" "2E_003" "2E_004" "2E_005" "2E_006" "2E_008" "2E_009" "2E_010" "2E_011" "2E_012" "2E_013" "2E_015" "2E_016" "2E_017" "2E_018" "2E_019" "2E_020" "2E_021" "2E_022" "2E_024" "2E_025" "2E_027" "2E_028" "2E_029" "2E_030" "2E_023" "2E_031" "2E_032" "2E_033" "2E_034" "2E_035" "2E_038" "2E_039" "2E_040" "2E_041" "2E_042" "2E_043" "2E_044" "2E_045" "2E_047" "2E_048" "2E_049" "2E_050" "2E_051" "2E_052" "2E_053" "2E_054" "2E_055" "2E_056" "2E_057" "2E_066" "2E_070" "2E_075" "2E_084" "2E_085" "2E_090" "2E_095" "2E_096" "2E_097" "2E_098" "2E_099" "2E_100" "2E_102" "2E_104" "2E_105" "2E_106" "2E_108" "2E_109" "2E_112" "2E_115" "2E_124" "2E_126" "2E_118" "2E_131" "2E_133" "2E_134")
#######################
# T1 processing
#######################
participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    file=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_anat-T1.nii.gz
    echo ${file}
    
    SRC=T1-preprocess-02_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q JM-GPU -pe smp 56 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/01_anat-T1_V2.sh ${participant_id} ${P_DIR} ${file}

done


#######################
# T2 processing
#######################
participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    file=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_anat-T2.nii.gz
    echo ${file}
    
    SRC=T2-preprocess-02_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q JM,CCOM,UI -pe smp 20 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/03_anat-T2_V2.sh ${participant_id} ${P_DIR} ${file}
done


#######################
# Registration to MNI
#######################
participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    
    SRC=ants-reg_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q JM,UI -pe smp 14 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/04_anat-MNI-registration_V2.sh ${participant_id}
done

### define these after looking at their registration output images
# choose a good method per participant and put their ID here
SyN_participants=("2E_001" "2E_002" "2E_003" "2E_004" "2E_005" "2E_006" "2E_008" "2E_009" "2E_010" "2E_011" "2E_013" "2E_015" "2E_016" "2E_017" "2E_018" "2E_019" "2E_020" "2E_021" "2E_022" "2E_024" "2E_029" "2E_030" "2E_023" "2E_031" "2E_032" "2E_033" "2E_034" "2E_035" "2E_038" "2E_039" "2E_040" "2E_041" "2E_042" "2E_043" "2E_044" "2E_045" "2E_047" "2E_048" "2E_049" "2E_050" "2E_051" "2E_052" "2E_053" "2E_054" "2E_055" "2E_056" "2E_057" "2E_066" "2E_070" "2E_075" "2E_084" "2E_085" "2E_090" "2E_095" "2E_096" "2E_097" "2E_098" "2E_099" "2E_100" "2E_102" "2E_104" "2E_105" "2E_106" "2E_108" "2E_109" "2E_022" "2E_029" "2E_124" "2E_126" "2E_118" "2E_131" "2E_133" "2E_134")
long_participants=("2E_012" "2E_025" "2E_027" "2E_112")

if [[ " ${SyN_participants[@]} " =~ " $participant_id " ]]; then
  echo "SyN"
fi


###################
# fMRI processing #
###################

#######################
# REST processing
#######################
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
cd ${P_DIR}

participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    file=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_fMRI-REST1.nii.gz
    
    task_NAME="REST1"
    rm -rf ${P_DIR}/data/derivatives/func/${participant_id}/run-2/${task_NAME}
    
    SRC=${task_NAME}-V2-preprocess_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q UI,JM -pe smp 12 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/07_func-preprocessing_V2.sh ${participant_id} ${P_DIR} ${file} ${task_NAME}
done

participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    file=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_fMRI-REST2.nii.gz
    
    task_NAME="REST2"
    rm -rf ${P_DIR}/data/derivatives/func/${participant_id}/run-2/${task_NAME}
    
    SRC=${task_NAME}-V2-preprocess_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q UI,JM,CCOM -pe smp 12 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/07_func-preprocessing_V2.sh ${participant_id} ${P_DIR} ${file} ${task_NAME}
done

#######################
# PS-VC processing
#######################
participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    file=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_fMRI-PSVC.nii.gz
    
    task_NAME="PS-VC"
    rm -rf ${P_DIR}/data/derivatives/func/${participant_id}/run-2/${task_NAME}
    
    SRC=${task_NAME}-V2-preprocess_${participant_id}
    rm logs/${SRC}.log
    qsub -cwd -q JM,UI,CCOM -pe smp 80 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/07_func-preprocessing_V2.sh ${participant_id} ${P_DIR} ${file} ${task_NAME}
done


## AFNI subject-level analysis
good_PSVC_R1=("2E_022" "2E_029" "2E_023" "2E_031" "2E_032" "2E_033" "2E_034" "2E_035" "2E_038" "2E_039" "2E_040" "2E_041" "2E_042" "2E_043" "2E_044" "2E_045" "2E_048" "2E_050" "2E_051" "2E_052" "2E_053" "2E_054" "2E_055" "2E_056" "2E_057" "2E_066" "2E_070" "2E_075" "2E_084" "2E_085" "2E_090" "2E_095" "2E_096" "2E_097" "2E_098" "2E_099" "2E_100" "2E_102" "2E_104" "2E_105" "2E_106" "2E_108" "2E_109" "2E_112" "2E_115" "2E_124" "2E_126" "2E_118" "2E_131" "2E_133" "2E_134")

for participant_id in "${good_PSVC_R1[@]}"; do
    echo ${participant_id}
    
    SRC=AFNI-PSVC_${participant_id}
    rm logs/${SRC}.log
    
    qsub -cwd -q JM,UI,CCOM -pe smp 20 \
        -N ${SRC} \
        -o logs/${SRC}.log \
        -j y -ckpt user \
        src/processing-pipelines/08_01_func-psvc-AFNI-subject-level-analysis.sh ${participant_id}
done


