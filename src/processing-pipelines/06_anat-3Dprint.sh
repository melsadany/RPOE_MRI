###############
# 3D printing #
###############
# making a mesh
# runs local on Topaz
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
cd ${P_DIR}

participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    
    sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/synthseg
    cd ${sub_DIR}
    
    
    flirt \
        -in ../E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz \
        -ref cortex-mask.nii.gz \
        -out resampled_brain.nii.gz \
        -applyxfm -usesqform
    fslmaths resampled_brain.nii.gz -mas cortex-mask.nii.gz brain.nii.gz
    
    /Dedicated/jmichaelson-wdata/msmuhammad/workbench/nii2mesh/src/nii2mesh brain.nii.gz -r 0.5 final.stl
done

