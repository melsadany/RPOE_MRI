source ~/.bashrc
conda deactivate

<< "argon"
conda activate theone
AFNI=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07
export PATH=$AFNI:$PATH
argon

SCHF_atlas=/Dedicated/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz

# prepare a folder of symlinks for each task
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
MAP_DIR=${P_DIR}/data/derivatives/func/AFNI-group-level-maps
task="word_association__word"
W_DIR=${MAP_DIR}/${task}/resid_to_lang
mkdir -p ${W_DIR}
cd ${W_DIR}

task_AFNI_res=${MAP_DIR}/${task}/${task}_MEMA_resZ+orig

### STEP1 ###
# Extract each subject's residual map using their actual name
for i in {0..50}; do
    subj=$(3dinfo -label "${task_AFNI_res}[$i]" | sed 's/-Res:Z//g')
    3dcalc -a "${task_AFNI_res}[$i]" -expr 'a' -prefix "res_${subj}.nii.gz"
done

### STEP2 ###
# Extract mean residuals from ROIs
3dROIstats -mask ${SCHF_atlas} res_*.nii.gz > roi_residuals.txt

