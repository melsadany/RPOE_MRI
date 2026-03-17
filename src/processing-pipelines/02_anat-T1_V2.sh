#!/bin/bash
# the qsub command 
# participant_id='2E_031';sub_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/anat/${participant_id};mkdir -p ${sub_DIR};inputScan=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/raw/clean/${participant_id}/sub-${participant_id}_anat-T1.nii.gz;qsub -cwd -q JM,UI,CCOM -pe smp 14 -N T1-preprocess_${participant_id} -o logs/T1-preprocess_${participant_id}.log -j y -ckpt user src/processing-pipelines/01_anat-T1.sh ${participant_id} ${sub_DIR} ${file}

module load stack/2022.2
module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:$PATH


# extract input data
participant_id=$1
P_DIR=$2
export inputScan=$3

: <<'trial'
participant_id='2E_134'
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
inputScan=${P_DIR}/data/raw/clean/${participant_id}/sub-${participant_id}_anat-T1.nii.gz
trial

sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2
mkdir -p ${sub_DIR}

# getting needed tools
conda activate ENA

# MNI reference paths
MNI_W=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm.nii.gz
MNI_B=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain.nii.gz
MNI_BM=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_0.5mm_brain_mask.nii.gz

# order of processing
# 1. reoriting to RPI orientation
# 2. Denoise
# 3. brain extraction #1
# 4. intensity non-uniformaty correction using N4 and mask #1
# 5. brain extraction again #2
# 6. N4 bias correction again #2
# 7. brain extraction again #3
# 8. brain extraction using ANTs #4
# 9. fast segmentation
# 10: SynthSeg brain segmentation 
# 11: STL file creation

cd ${sub_DIR}
outputNameBase=sub-${participant_id}_anat-T1

if [ ! -d $sub_DIR ]; then
  mkdir $sub_DIR
elif [ -f ${sub_DIR}/F_${outputNameBase}_ANTs-brain.nii.gz ]; then
  echo "Preprocessing has already been run"
  exit 0
fi


# Step 1: reorient to RPI
3dresample \
    -orient rpi \
    -overwrite \
    -input $inputScan \
    -prefix A_${outputNameBase}_reorient_RPI.nii.gz
if [ -f A_${outputNameBase}_reorient_RPI.nii.gz ]; then
    echo "------Done with RPI reorientation"
fi



# Step 2: denoising
# not including mask in denoising since it hasn't been generated yet
# Denoising done by ANTs
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



# Step 3: brain mask #1
#pre_mask=C_${outputNameBase}_MNI-rigid-reg.nii.gz
pre_mask=B_${outputNameBase}_reorient_RPI_denoise.nii.gz
#3dAutomask \
#    -prefix D_${outputNameBase}_mask-AUTO.nii.gz \
#    -clfrac 0.2 -q \
#    ${pre_mask}
#CopyImageHeaderInformation ${pre_mask} \
#    D_${outputNameBase}_mask-AUTO.nii.gz \
#    D_${outputNameBase}_mask-AUTO.nii.gz 1 1 1
hd-bet -i ${pre_mask} -o D_${outputNameBase}_hd-bet-brain.nii.gz
if [ -f D_${outputNameBase}_hd-bet-brain.nii.gz ]; then
    echo "------Done with initial HD-BET mask #1"
fi



# Step 4: bias correction #1
# ANTs intensity correction
N4BiasFieldCorrection \
    -d 3 \
    -i ${pre_mask} \
    -x D_${outputNameBase}_hd-bet-brain.nii.gz \
    -r 1 -s 4 \
    -c [50x50x50x50,0.0] \
    -b [200,3] \
    -t [0.15,0.01,200] \
    -o [E_${outputNameBase}_ANTs-bfcorr.nii.gz,E_${outputNameBase}_ANTs-bf.nii.gz]
if [ -f E_${outputNameBase}_ANTs-bfcorr.nii.gz ]; then
    echo "------Done with N4Bias correction #1"
fi



# Step 5: brain mask #2
hd-bet -i E_${outputNameBase}_ANTs-bfcorr.nii.gz -o D_${outputNameBase}_hd-bet-brain.nii.gz
if [ -f D_${outputNameBase}_hd-bet-brain.nii.gz ]; then
    echo "------Done with HD-BET mask #2"
fi



# Step 6: bias correction #2
# ANTs intensity correction
N4BiasFieldCorrection \
    -d 3 \
    -i E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    -x D_${outputNameBase}_hd-bet-brain.nii.gz \
    -r 1 -s 4 \
    -c [50x50x50x50,0.0] \
    -b [200,3] \
    -t [0.15,0.01,200] \
    -o [E_${outputNameBase}_ANTs-bfcorr.nii.gz,E_${outputNameBase}_ANTs-bf.nii.gz]
if [ -f E_${outputNameBase}_ANTs-bfcorr.nii.gz ]; then
    echo "------Done with N4Bias correction #2"
fi


# AFNI intensity correction; takes a long time but good output
: <<"AFNIBC"
3dUnifize \
    -prefix E_${outputNameBase}_AFNI-bfcorr.nii.gz \
    -input E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    -GM \
    -Urad 30
if [ -f E_${outputNameBase}_AFNI-bfcorr.nii.gz ]; then
    echo "------Done with 3dUnifize Bias correction"
fi
AFNIBC


# Step 7: brain mask #3
hd-bet -i E_${outputNameBase}_ANTs-bfcorr.nii.gz -o F_${outputNameBase}_hd-bet-brain.nii.gz
if [ -f F_${outputNameBase}_hd-bet-brain.nii.gz ]; then
    echo "------Done with HD-BET #3"
fi




# Step 8: ANTs brain mask #4
LABEL=brain
TEMPLATE="OASIS" 
DIR_TEMPLATE=/Dedicated/jmichaelson-wdata/msmuhammad/refs/${TEMPLATE}

antsBrainExtraction.sh \
    -d 3 \
    -a E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    -e ${DIR_TEMPLATE}/T_template0.nii.gz \
    -m ${DIR_TEMPLATE}/T_template0_BrainCerebellumProbabilityMask.nii.gz \
    -f ${DIR_TEMPLATE}/T_template0_BrainCerebellumRegistrationMask.nii.gz \
    -o F_${outputNameBase}_ANTs-mask-


CopyImageHeaderInformation E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    F_${outputNameBase}_ANTs-mask-BrainExtractionMask.nii.gz \
    F_${outputNameBase}_ANTs-mask-BrainExtractionMask.nii.gz 1 1 1
mv F_${outputNameBase}_ANTs-mask-BrainExtractionMask.nii.gz \
    F_${outputNameBase}_ANTs-mask-.nii.gz
mv F_${outputNameBase}_ANTs-mask-BrainExtractionBrain.nii.gz \
    F_${outputNameBase}_ANTs-brain.nii.gz
if [ -f F_${outputNameBase}_ANTs-brain.nii.gz ]; then
    echo "------Done with ANTs brain extraction"
fi



# Step 9: FAST tissue segmentation
# for HD-BET brain
fast -v F_${outputNameBase}_hd-bet-brain.nii.gz

# extract volumetrics of these tissues
# do this locally. the function on argon isn't the same as the one locally.
label_stats_output=F_${outputNameBase}_fast-volumetric-stats.txt
LabelGeometryMeasures 3 F_${outputNameBase}_hd-bet-brain_seg.nii.gz > ${label_stats_output}


#########
echo "------------DONE------------"





###############################################################################################
###############################################################################################
###############################################################################################
########################################### PART 2 ############################################
###############################################################################################
###############################################################################################
###############################################################################################

module load stack/2022.2
module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH

# extract input data
participant_id=$1
P_DIR=$2

: <<'trial'
participant_id='2E_124'
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
trial

sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2
mkdir -p ${sub_DIR}

# getting needed tools
conda activate ENA

# order of processing
# 10: SynthSeg brain segmentation 
# 11: STL file creation

cd ${sub_DIR}
outputNameBase=sub-${participant_id}_anat-T1


# Step 10: SynthSeg brain segmentation
# it needs Freesurfer with a new version. at least 7.4.1

export FREESURFER_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1
export FSFAST_HOME=/Users/msmuhammad/packages/FreeSurfer/7.4.1/fsfast
export SUBJECTS_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/subjects
export MNI_DIR=/Users/msmuhammad/packages/FreeSurfer/7.4.1/mni
export PATH=${FREESURFER_HOME}/bin:$PATH

mkdir -p synthseg2
mri_synthseg \
    --i E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    --parc \
    --o synthseg2/ \
    --vol synthseg2/vols.csv \
    --threads 4

# Step 11: STL mesh
# extract cortex regions first
mri_binarize --i synthseg/E_${outputNameBase}_ANTs-bfcorr_synthseg.nii.gz \
    --match 2 3 31 41 42 63 72 77 51 52 13 12 43 50 4 11 26 58 49 10 17 18 53 54 44 5 80 14 30 62 \
    --o synthseg/cortex-mask.nii.gz



#########
echo "------------DONE------------"

###############################################################################################
###############################################################################################
###############################################################################################
########################################### PART 3 ############################################
###############################################################################################
###############################################################################################
###############################################################################################

module load stack/2022.2
module load ants
module load FSL 
export PATH=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07:/Users/msmuhammad/packages/FreeSurfer/7.4.1:$PATH

# extract input data
participant_id=$1
P_DIR=$2

: <<'trial'
participant_id='2E_118'
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
trial

sub_DIR=${P_DIR}/data/derivatives/anat/${participant_id}/run-2
mkdir -p ${sub_DIR}

# getting needed tools
conda activate ENA

# order of processing
# 12: DL+DiReCT cortical thickness and processing

cd ${sub_DIR}
outputNameBase=sub-${participant_id}_anat-T1


# Step 12: DL+DiReCT brain segmentation and cortical thickness

mkdir -p DL_DiReCT
dl+direct \
    --subject ${participant_id} \
    --bet E_${outputNameBase}_ANTs-bfcorr.nii.gz \
    DL_DiReCT

slicer DL_DiReCT/T1w_norm.nii.gz DL_DiReCT/T1w_norm_seg.nii.gz \
    -a DL_DiReCT/T1w_norm_seg.png



#########
echo "------------DONE------------"



