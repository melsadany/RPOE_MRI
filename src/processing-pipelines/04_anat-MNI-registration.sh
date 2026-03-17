# this script is for registering the T1w to MNI space

module load stack/2022.2
module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH

participant_id=$1

: <<"trial"
participant_id='2E_126'
trial


### ANTs registration for all participants
MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz
MNI_B=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain.nii.gz
MNI_BM=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain_mask.nii.gz
MNI_L=/Dedicated/jmichaelson-wdata/msmuhammad/refs/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_reor.nii.gz

##
echo $participant_id
##


P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
sub_DIR=${P_DIR}/data/derivatives/anat-registration/${participant_id}
mkdir -p ${sub_DIR}
cd ${sub_DIR}
T1_IMAGE=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz
T1_BM=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/F_sub-${participant_id}_anat-T1_hd-bet-brain_mask.nii.gz


### clean the folder
mkdir -p ${sub_DIR}/archive
mv ${sub_DIR}/*.nii.gz ${sub_DIR}/archive/.
mv ${sub_DIR}/*.txt ${sub_DIR}/archive/.
mv ${sub_DIR}/*.png ${sub_DIR}/archive/.
mv ${sub_DIR}/*.mat ${sub_DIR}/archive/.
mv ${sub_DIR}/*.h5 ${sub_DIR}/archive/.

    
    
### ANTs SyN registration
antsRegistrationSyN.sh \
     -d 3 \
     -f ${MNI_W} \
     -m ${T1_IMAGE} \
     -o ${sub_DIR}/sub-${participant_id}-T1-to-MNI-SyN_ \
     -x ${MNI_BM} \
     -n 12
cp ${sub_DIR}/sub-${participant_id}-T1-to-MNI-SyN_Warped.nii.gz ${sub_DIR}/sub-${participant_id}-T1-to-MNI-SyN.nii.gz


# make a PNG of the registration result and save it
slicer ${MNI_W} sub-${participant_id}-T1-to-MNI-SyN.nii.gz -a ${participant_id}_reg-SyN.png


### ANTs long registration
#scriptDir=/Shared/pinc/sharedopt/apps
#module load stack/2022.2
antsRegistration \
    --dimensionality 3 \
    --output ${sub_DIR}/sub-${participant_id}-T1-to-MNI-long_ \
    --initial-moving-transform [${MNI_W},${T1_IMAGE},1] \
    --transform Rigid[0.1] \
    --metric Mattes[${MNI_W},${T1_IMAGE},1,32,Regular,0.25] \
    --convergence [2000x2000x2000x2000x2000,1e-6,10] \
    --smoothing-sigmas 4x3x2x1x0vox \
    --shrink-factors 8x8x4x2x1 \
    --transform Affine[0.1] \
    --metric Mattes[${MNI_W},${T1_IMAGE},1,32,Regular,0.25] \
    --convergence [1000x1000x1000x1000x1000,1e-6,10] \
    --smoothing-sigmas 4x3x2x1x0vox \
    --shrink-factors 8x8x4x2x1 \
    --transform Syn[0.1,3,0] \
    --metric CC[${MNI_W},${T1_IMAGE},1,4,Regular,0.25] \
    --convergence [1000x1000x1000x1000x1000,1e-6,10] \
    --smoothing-sigmas 4x3x2x1x0vox \
    --shrink-factors 8x8x4x2x1 \
    --use-histogram-matching 1 \
    --verbose 1 \
    --random-seed 13983981 \
    --winsorize-image-intensities [0.005,0.995] \
    --write-composite-transform 1



# apply transformation matrix
antsApplyTransforms \
    -d 3 \
    -n BSpline[3] \
    -i ${T1_IMAGE} \
    -o ${sub_DIR}/sub-${participant_id}-T1-to-MNI-long.nii.gz \
    -t ${sub_DIR}/sub-${participant_id}-T1-to-MNI-long_0GenericAffine.mat \
    -r ${MNI_W}


slicer ${MNI_W} sub-${participant_id}-T1-to-MNI-long.nii.gz -a ${participant_id}_reg-long.png


######################################
# corpus callosum and other DWI ROIs #
######################################
# 

# Step 1: Apply the inverse transformation to bring labeled regions back to the native space

# I needed to resample the labeled nii to match my reference of MNI
# only donr once
#FMRIB_labels=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/atlases/JHU/JHU-ICBM-labels
#3dresample \
#    -master $MNI_W \
#    -input ${FMRIB_labels}-1mm.nii.gz \
#    -prefix ${FMRIB_labels}-0.5mm_ME.nii.gz \
#    -rmode NN

FMRIB_labels=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/atlases/JHU/JHU-ICBM-labels-0.5mm_ME.nii.gz

# apply the SyN inverse
warp=${sub_DIR}/sub-${participant_id}-T1-to-MNI-SyN_1Warp.nii.gz
affine=${sub_DIR}/sub-${participant_id}-T1-to-MNI-SyN_0GenericAffine.mat
antsApplyTransforms \
    -d 3 \
    -i ${FMRIB_labels} \
    -r ${T1_IMAGE} \
    -o ${sub_DIR}/sub-${participant_id}-FMRIB-labels_SyN.nii.gz \
    -t [${affine},1] \
    -t ${warp}
  
# apply the long inverse
antsApplyTransforms \
    -d 3 \
    -i ${FMRIB_labels} \
    -r ${T1_IMAGE} \
    -o ${sub_DIR}/sub-${participant_id}-FMRIB-labels_long.nii.gz \
    -t [${sub_DIR}/sub-${participant_id}-T1-to-MNI-long_0GenericAffine.mat,1] \
    -v

methods=("SyN" "long")
for method in "${methods[@]}"; do
    labeled=${sub_DIR}/sub-${participant_id}-FMRIB-labels_${method}.nii.gz
    
    OUTPUT_CSV=FMRIB-vols-${method}.csv
    LabelGeometryMeasures 3 \
        ${labeled} > ${OUTPUT_CSV}
        
    echo "Volume extraction complete for ${method}"

done


#########
echo "------------DONE------------"


