#!/bin/bash

module load stack/2022.2
module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH

# extract input data
participant_id=$1
P_DIR=$2
export inputScan=$3

: <<'trial'
participant_id=2E_134
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
inputScan=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_anat-T2.nii.gz
trial

sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2/T2
mkdir -p ${sub_DIR}

# getting needed tools
conda activate ENA

# MNI reference paths
MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz
MNI_B=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain.nii.gz
MNI_BM=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain_mask.nii.gz

cd ${sub_DIR}
outputNameBase=sub-${participant_id}_anat-T2


# order of processing
# 1. reoriting to RPI orientation
# 2. Denoise
# 3. rigid alignment to T1
# 4. intensity non-uniformaty correction using N4 and T1 brain mask


if [ ! -d $sub_DIR ]; then
  mkdir $sub_DIR
elif [ -f ${sub_DIR}/D_${outputNameBase}_anat-T2_ANTs-bfcorr ]; then
  echo "Preprocessing has already been run"
  exit 0
fi


################################################################################
################################################################################
################################################################################
# Step 1: reorient to RPI
3dresample \
    -orient rpi \
    -overwrite \
    -input $inputScan \
    -prefix A_${outputNameBase}_reorient_RPI.nii.gz
if [ -f A_${outputNameBase}_reorient_RPI.nii.gz ]; then
    echo "------Done with RPI reorientation"
fi


################################################################################
################################################################################
################################################################################
# Step 2: denoising
# not including mask in denoising since it hasn't been generated yet
DenoiseImage \
    -d 3 \
    -n Rician \
    -s 1 -p 1 -r 2 -v 1 \
    -i A_${outputNameBase}_reorient_RPI.nii.gz \
    -o [B_${outputNameBase}_reorient_RPI_denoise.nii.gz,B_${outputNameBase}_reorient_RPI_noise.nii.gz]


: <<'gibbsfix'
# some participants had gibbs ringing, so I tried doing some smoothing hoping it can fix that
fslmaths B_${outputNameBase}_reorient_RPI_denoise.nii.gz \
    -s 0.5 \
    B_${outputNameBase}_reorient_RPI_denoise.nii.gz
gibbsfix

if [ -f B_${outputNameBase}_reorient_RPI_denoise.nii.gz ]; then
    echo "------Done with denoising step"
fi


################################################################################
################################################################################
################################################################################
# Step 3: rigid alignment to T1 scan
# run ants rigid registration on data
T1_scan=${sub_DIR}/../E_sub-${participant_id}_anat-T1_ANTs-bfcorr.nii.gz
antsRegistration \
  --dimensionality 3 \
  --output [C_${outputNameBase}_T22T1-SyNDC_, C_${outputNameBase}_T22T1-SyNDC_Warped.nii.gz] \
  --interpolation Linear \
  --use-histogram-matching 1 \
  --initial-moving-transform [${T1_scan}, B_${outputNameBase}_reorient_RPI_denoise.nii.gz, 1] \
  --transform Rigid[0.1] \
  --metric MI[${T1_scan}, B_${outputNameBase}_reorient_RPI_denoise.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [1000x500x250x100,1e-6,10] \
  --shrink-factors 8x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox \
  --transform SyN[0.1,3,0] \
  --metric MI[${T1_scan}, B_${outputNameBase}_reorient_RPI_denoise.nii.gz, 1, 32, Regular, 0.25] \
  --convergence [100x70x50x20,1e-6,10] \
  --shrink-factors 6x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox \
  --restrict-deformation 0x1x0

slicer ${T1_scan} C_${outputNameBase}_T22T1-SyNDC_Warped.nii.gz -a C_T2-T1-reg-overlay.png
if [ -f C_${outputNameBase}_T22T1-SyNDC_Warped.nii.gz ]; then
    echo "------Done with registration to T1"
fi

################################################################################
################################################################################
################################################################################
# Step 4: bias correction #1
# ANTs intensity correction
T1_BM=${sub_DIR}/../F_sub-${participant_id}_anat-T1_hd-bet-brain_mask.nii.gz
N4BiasFieldCorrection \
    -d 3 \
    -i C_${outputNameBase}_T22T1-SyNDC_Warped.nii.gz \
    -x ${T1_BM} \
    -r 1 -s 4 \
    -c [50x50x50x50,0.0] \
    -b [200,3] \
    -t [0.15,0.01,200] \
    -o [D_${outputNameBase}_ANTs-bfcorr.nii.gz,D_${outputNameBase}_ANTs-bf.nii.gz]
if [ -f D_${outputNameBase}_ANTs-bfcorr.nii.gz ]; then
    echo "------Done with N4Bias correction"
fi


#########
echo "------------DONE------------"


