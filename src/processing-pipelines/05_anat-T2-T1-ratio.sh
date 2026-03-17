#################
# T1w/T2w ratio #
#################
# runs local on Topaz
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
cd ${P_DIR}

# Function to limit the number of parallel jobs
limit_jobs() {
    local max_jobs=$1
    while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
        sleep 0.1
    done
}
MAX_JOBS=20

participants=("2E_131" "2E_133" "2E_134")
for participant_id in "${participants[@]}"; do
    echo ${participant_id}
    
    sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2
    cd ${sub_DIR}
    
    T2=T2/D_sub-${participant_id}_anat-T2_ANTs-bfcorr.nii.gz
    T1=E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz
    T1_BM=F_sub-${participant_id}_anat-T1_hd-bet-brain_mask.nii.gz
    
    rm -rf anat-metrics2
    mkdir -p anat-metrics2
    cd anat-metrics2
    
    # extract brain
    fslmaths ../$T1 -mul ../$T1_BM t1w-brain.nii.gz
    fslmaths ../$T2 -mul ../$T1_BM t2w-brain.nii.gz
    
    
    ## clip the values making sure there's no negative values
    3dcalc -a t1w-brain.nii.gz -expr 'step(a)*a' -prefix t1w-brain_clipped.nii.gz
    3dcalc -a t2w-brain.nii.gz -expr 'step(a)*a' -prefix t2w-brain_clipped.nii.gz
    
    ## normalize
    min_t1=$(3dmaskave -quiet -min -mask ../$T1_BM t1w-brain_clipped.nii.gz)
    max_t1=$(3dmaskave -quiet -max -mask ../$T1_BM t1w-brain_clipped.nii.gz)
    min_t2=$(3dmaskave -quiet -min -mask ../$T1_BM t2w-brain_clipped.nii.gz)
    max_t2=$(3dmaskave -quiet -max -mask ../$T1_BM t2w-brain_clipped.nii.gz)
    
    3dcalc -a t1w-brain_clipped.nii.gz -b ../$T1_BM \
        -expr "b * (a - ${min_t1}) / (${max_t1} - ${min_t1})" \
        -prefix t1w-brain_clipped_norm.nii.gz
    3dcalc -a t2w-brain_clipped.nii.gz -b ../$T1_BM \
        -expr "b * (a - ${min_t2}) / (${max_t2} - ${min_t2})" \
        -prefix t2w-brain_clipped_norm.nii.gz
    
    
    ## divide and save
    3dcalc -a t1w-brain_clipped_norm.nii.gz -b t2w-brain_clipped_norm.nii.gz -c ../$T1_BM -expr 'c * (a/b)' -prefix t1w-t2w_ratio.nii.gz
    
    
    # resample the labeled brain mask
    labeled_BM=../synthseg2/E_sub-${participant_id}_anat-T1_ANTs-bfcorr_synthseg.nii.gz
    labeled_resampled=E_sub-${participant_id}_anat-T1_ANTs-bfcorr_synthseg-resampled.nii.gz
    flirt -in $labeled_BM \
        -ref ../$T1 \
        -out $labeled_resampled \
        -applyxfm -usesqform -interp nearestneighbour
    
    # extract summary stats per ROI and save
    output_csv=sub-${participant_id}_t1w-t2w-summary-stats_2.csv
    3dROIstats -mask ${labeled_resampled} -nzmedian -nzmean -nzminmax t1w-t2w_ratio.nii.gz > $output_csv
done