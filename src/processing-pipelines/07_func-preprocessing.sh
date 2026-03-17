# this script is for preprocessing the fMRI data

module load stack/2022.2
#module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/ants/2.5.1/bin:$PATH


export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=80
export OMP_NUM_THREADS=80

###### input data 
participant_id=$1
P_DIR=$2
task_IMAGE=$3
task_NAME=$4

###### trial data
<< "trial"
participant_id=2E_134
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
task_NAME="REST1"
#task_NAME="PS-VC"
task_IMAGE=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_fMRI-${task_NAME}.nii.gz
trial
######
# participant-specific data
sub_DIR=${P_DIR}/data/derivatives/func/${participant_id}/run-3/${task_NAME}
mkdir -p ${sub_DIR}
cd ${sub_DIR}
T1_IMAGE=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz
T2_IMAGE=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/T2/F_sub-${participant_id}_anat-T2_T1-rigid-reg.nii.gz
T1_MASK=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/F_sub-${participant_id}_anat-T1_hd-bet-brain_mask.nii.gz
T1_BRAIN=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/F_sub-${participant_id}_anat-T1_hd-bet-brain.nii.gz

###### standard data
# get standard reference brain 
MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_2mm.nii.gz
MNI_BRAIN=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_2mm_brain.nii.gz
MNI_BRAIN_M=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_2mm_brain_mask.nii.gz


###### reference volume for processes 
# get the middle volume in the functional series 
#halfPoint=$(fslhd ${task_IMAGE} | grep "^dim4" | awk '{print int($2/2)}')

# Extract the TR from the image json
extract_tr=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/extract_tr.sh
TR=$(${extract_tr} ${task_IMAGE%.nii.gz}.json)

out_basename=sub-${participant_id}_${task_NAME}


#########
# STEP1 #
#########
# slice-timing correction

# It is a debate of when to do slicetiming correction (whether before or after motion correction)
# I relied on the answer here: https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/fmri/preprocessing/slice_timing/ 
# and here https://neurostars.org/t/why-motion-correction-first-fmriprep/2844

# either use the extracted slice-timing vector or the standard one in the data directory
# some participants didn't have the PSVC slice timing vector written in their JSON    

# unified
slicetimer -i ${task_IMAGE} \
    -o A_${out_basename}_slice-timing-corrected.nii.gz \
    -r ${TR} \
    --tcustom=${P_DIR}/data/slice-timing_${task_NAME}.txt

if [ -f A_${out_basename}_slice-timing-corrected.nii.gz ]; then
    echo "------Done with slice timing correction"
fi

#########
# STEP2 #
#########
# RPI reorientation 
3dresample \
    -orient rpi \
    -overwrite \
    -input A_${out_basename}_slice-timing-corrected.nii.gz \
    -prefix B_${out_basename}_reorient-RPI.nii.gz

if [ -f B_${out_basename}_reorient-RPI.nii.gz ]; then
    echo "------Done with RPI reorientation"
fi

# deoblique is mainly to align the image more squarely in the native space
# deoblique doesn't overwrite existing file, so delete if it exists
if [ -f C_${out_basename}_deoblique.nii.gz ]; then
    rm C_${out_basename}_deoblique.nii.gz
fi
3dWarp \
    -deoblique \
    -prefix C_${out_basename}_deoblique.nii.gz \
    B_${out_basename}_reorient-RPI.nii.gz

if [ -f C_${out_basename}_deoblique.nii.gz ]; then
    echo "------Done with the deoblique step"
fi



#########
# STEP3 #
#########
# denoise the image

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=20
export OMP_NUM_THREADS=20
3dDespike \
    -prefix D_${out_basename}_denoise.nii.gz \
    C_${out_basename}_deoblique.nii.gz

if [ -f D_${out_basename}_denoise.nii.gz ]; then
    echo "------Done with the denoising step"
fi



#########
# STEP4 #
#########
# motion correction

## extract a reference volume
ref_vol=E_${out_basename}_ref-vol.nii.gz
fslroi D_${out_basename}_denoise.nii.gz ${ref_vol} 0 1

TS=D_${out_basename}_denoise.nii.gz
TS_MC=F_${out_basename}_moco.nii.gz
TS_null=F_${out_basename}_null.nii.gz

antsMotionCorr \
  -d 3 -o [${TS_null},${TS_MC},${ref_vol}] \
  -m MI[${ref_vol},${TS},1,32,Regular,0.2] \
  -t Rigid[0.25] \
  -i 50x25x10 \
  -s 2x1x0 \
  -f 3x2x1 \
  -u 1 -e 1

# get the first volume of the moco
TS_first_vol=F_${out_basename}_moco-first-vol.nii.gz
fslroi ${TS_MC} ${TS_first_vol} 0 1
  

if [[ -f ${TS_MC} ]]; then
    echo "------Done with the motion correction step"
fi



#########
# STEP5 #
#########
# get brain mask
TS_mask=G_${out_basename}_first-vol-masked
bet ${TS_first_vol} ${TS_mask} -m -v
slicer ${TS_first_vol} ${TS_mask}_mask.nii.gz -a G_bet-mask-overlay

# mask the entire 4D image
TS_MC_masked=G_${out_basename}_moco-masked.nii.gz
fslmaths ${TS_MC} -mas ${TS_mask}_mask.nii.gz ${TS_MC_masked}


if [[ -f ${TS_MC_masked} ]]; then
    echo "------Done with the brain mask step"
fi



#########
# STEP6 #
#########
# coreg to native space
antsRegistration \
  --dimensionality 3 \
  --output [H_${out_basename}_func2T1_SyNDC_, H_${out_basename}_func2T1_SyNDC_Warped.nii.gz] \
  --interpolation Linear \
  --use-histogram-matching 1 \
  --initial-moving-transform [${T1_BRAIN}, ${TS_mask}.nii.gz, 1] \
  --transform Rigid[0.1] \
  --metric MI[${T1_BRAIN}, ${TS_mask}.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [1000x500x250x100,1e-6,10] \
  --shrink-factors 8x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox \
  --transform SyN[0.1,3,0] \
  --metric MI[${T1_BRAIN}, ${TS_mask}.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [100x70x50x20,1e-6,10] \
  --shrink-factors 6x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox \
  --restrict-deformation 0x1x0
slicer ${T1_IMAGE} H_${out_basename}_func2T1_SyNDC_Warped.nii.gz -a H_TS-T1-SyNDC-overlay.png



if [[ -f H_${out_basename}_func2T1_SyNDC_Warped.nii.gz ]]; then
    echo "------Done with the registration to native step"
fi




#########
# STEP7 #
#########
# reg to ref atlas


## func to T1 to MNI
antsRegistration \
  --dimensionality 3 \
  --output [I_${out_basename}_func2T1_SyNDC2MNI-SyN_, I_${out_basename}_func2T1_SyNDC2MNI-SyN_Warped.nii.gz] \
  --interpolation Linear \
  --use-histogram-matching 1 \
  --initial-moving-transform [${MNI_BRAIN}, H_${out_basename}_func2T1_SyNDC_Warped.nii.gz, 1] \
  --transform Rigid[0.1] \
  --metric MI[${MNI_BRAIN}, H_${out_basename}_func2T1_SyNDC_Warped.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [1000x500x250x100,1e-6,10] \
  --shrink-factors 8x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox \
  --transform SyN[0.1,3,0] \
  --metric MI[${MNI_BRAIN},H_${out_basename}_func2T1_SyNDC_Warped.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [100x70x50x20,1e-6,10] \
  --shrink-factors 6x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox
slicer ${MNI_BRAIN} I_${out_basename}_func2T1_SyNDC2MNI-SyN_Warped.nii.gz -S 2 1000 I_${out_basename}_func2T1_SyNDC2MNI-SyN-overlay.png


antsApplyTransforms \
    -d 3 \
    -e 3 \
    -n Linear \
    -i ${TS_MC_masked} \
    -o J_${out_basename}_MNI-reg_final.nii.gz \
    -t I_${out_basename}_func2T1_SyNDC2MNI-SyN_1Warp.nii.gz \
    -t I_${out_basename}_func2T1_SyNDC2MNI-SyN_0GenericAffine.mat \
    -t H_${out_basename}_func2T1_SyNDC_1Warp.nii.gz \
    -t H_${out_basename}_func2T1_SyNDC_0GenericAffine.mat \
    -r ${MNI_BRAIN} \
    -v



if [[ -f J_${out_basename}_MNI-reg_final.nii.gz ]]; then
    echo "------Done with the registration to MNI 2mm step"
fi





#########
# STEP8 #
#########
# extract ReHo 
3dReHo -prefix K_${out_basename}_MNI-reg-ReHo.nii.gz \
       -inset J_${out_basename}_MNI-reg_final.nii.gz \
       -nneigh 27
SCHF_atlas=/Dedicated/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz
3dROIstats -mask ${SCHF_atlas} K_${out_basename}_MNI-reg-ReHo.nii.gz > K_${out_basename}_MNI-reg-ReHo_ROI-stats.txt




#########
# STEP9 #
#########
# extract fALFF 
if [ $task_NAME != "PS-VC" ]; then
    3dRSFC -input J_${out_basename}_MNI-reg_final.nii.gz \
        -prefix L_${out_basename}_MNI-reg-RSFC \
        -band 0.01 0.08 \
        -nodetrend -mask ${MNI_BRAIN_M}
    3dTcat -prefix L_${out_basename}_MNI-reg-fALFF.nii.gz \
        "L_${out_basename}_MNI-reg-RSFC_fALFF+orig.BRIK[0]"
    3dTcat -prefix L_${out_basename}_MNI-reg-ALFF.nii.gz \
        "L_${out_basename}_MNI-reg-RSFC_ALFF+orig.BRIK[0]"
    3dTcat -prefix L_${out_basename}_MNI-reg-mALFF.nii.gz \
        "L_${out_basename}_MNI-reg-RSFC_mALFF+orig.BRIK[0]"
    3dTcat -prefix L_${out_basename}_MNI-reg-RSFA.nii.gz \
        "L_${out_basename}_MNI-reg-RSFC_RSFA+orig.BRIK[0]"        
    3dTcat -prefix L_${out_basename}_MNI-reg-mRSFA.nii.gz \
        "L_${out_basename}_MNI-reg-RSFC_mRSFA+orig.BRIK[0]"
    rm L_${out_basename}_MNI-reg-RSFC_*.HEAD
fi



#########
echo "------------DONE------------"

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
### CLEAN

# participant-specific data
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
task_IMG=${sub_DIR}/G_sub-${participant_id}_${task_NAME}_moco-masked.nii.gz
reg_task_IMG=${sub_DIR}/J_sub-${participant_id}_${task_NAME}_MNI-reg_final.nii.gz
T1_IMG=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz


O_sub_DIR=${sub_DIR}
sub_DIR=${P_DIR}/data/derivatives/func/USE_THIS/${participant_id}/${task_NAME}
rm -rf ${sub_DIR}
mkdir -p ${sub_DIR}
cd ${sub_DIR}


out_basename=sub-${participant_id}_${task_NAME}
########################################################################################################
# copy the registered task nii here
cp ${task_IMG} A_${out_basename}.nii.gz
cp ${reg_task_IMG} B_${out_basename}_MNI-reg.nii.gz
cp ${T1_IMG} C_sub-${participant_id}_T1.nii.gz
########################################################################################################


