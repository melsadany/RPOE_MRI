source ~/.bashrc
conda deactivate
###
# AFNI relies on $PATH variable, so make sure you start in a new terminal
###
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
MAP_DIR=${P_DIR}/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs


task="fALFF"
W_DIR=${MAP_DIR}/${task}
cd ${W_DIR}

###################
mkdir -p $W_DIR/masked $W_DIR/final

## mask the input
for subj in $(cat subj-list.txt); do
  echo $subj
  3dcalc -a "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/${subj}/run-3/REST1/L_sub-${subj}_REST1_MNI-reg-fALFF.nii.gz" \
    -b final_mask.nii.gz \
    -expr 'a*b' \
    -prefix masked/${subj}_masked.nii.gz
done

## Add tiny constant to avoid exact zeros
for subj in $(cat subj-list.txt); do
  echo $subj
  3dcalc -a masked/${subj}_masked.nii.gz \
    -expr 'a + 0.000001' \
    -prefix final/${subj}.nii.gz
done

## Run analysis
3dMVM -prefix ${task}_MVM_shared -jobs 40 \
  -mask final_mask.nii.gz \
  -bsVars "MRI_age+sex+shared_PC1" \
  -qVars "MRI_age,shared_PC1" \
  -qVarCenters Z -SS_type 3 \
  -num_glt 1 \
  -gltLabel 1 shared_PC1_t -gltCode 1 'shared_PC1 : 1' \
  -dataTable @AFNI-table.txt


###################
## get files
INPUT_FILES=$(tail -n +2 AFNI-table.txt | awk '{print $14}' | tr '\n' ' ')

## run the multi-regression model
3dMVM -prefix ${task}_multi -jobs 40 \
   -mask final_mask.nii.gz \
   -bsVars "MRI_age+sex+shared_PC1+lang_specific+iq_specific" \
   -qVars "MRI_age,shared_PC1,lang_specific,iq_specific" \
   -qVarCenters Z -SS_type 3 \
   -num_glt 3 \
   -gltLabel 1 shared_PC1_t -gltCode 1 'shared_PC1 : 1' \
   -gltLabel 2 lang_specific_t -gltCode 2 'lang_specific : 1' \
   -gltLabel 3 iq_specific_t -gltCode 3 'iq_specific : 1' \
   -dataTable @AFNI-table.txt
   
   
## extract nifti maps
3dcalc -a ${task}_multi+orig'[3]' -expr 'a' -prefix stat_multi_shared_PC1.nii.gz
3dcalc -a ${task}_multi+orig'[3]' -prefix pval_multi_shared_PC1.nii.gz \
    -expr '2*min(fitt_t2p(a,45), 1-fitt_t2p(a,45))' 
3dFDR -input pval_multi_shared_PC1.nii.gz -prefix FDR_multi_shared_PC1.nii.gz


3dMVM -prefix ${task}_shared -jobs 40 \
   -mask tight_mask.nii.gz \
   -bsVars "MRI_age+sex+shared_PC1" \
   -qVars "MRI_age,shared_PC1" \
   -qVarCenters Z -SS_type 3 \
   -num_glt 1 \
   -gltLabel 1 shared_PC1_t -gltCode 1 'shared_PC1 : 1' \
   -dataTable @AFNI-table.txt


sub_list="AFNI-table.txt"; set_param=""
while read -r subj v1 v2 v3 v4 v5 file; do
    set_param+=" ${subj} ${file}"
done < <(tail -n +2 "$sub_list")

3dttest++ -prefix ${task}_shared_tt -mask final_mask.nii.gz \
  -setA ${task} ${set_param} \
  -covariates ../age-sex-covar.txt \
  -regCov "shared_PC1" \
  -resid ${task}_shared_tt_resid \
  -mask final_mask.nii.gz \
  -Clustsim 
  
  
### viz
export SUBJECTS_DIR=/wdata/msmuhammad/refs/mni-freesurfer
SUBJECT_ID="2mm"
TEMP_DIR=$SUBJECTS_DIR/$SUBJECT_ID
@SUMA_Make_Spec_FS -sid ${SUBJECT_ID} -fspath ${TEMP_DIR}/surf

P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
W_DIR=${P_DIR}/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF

mkdir -p viz
cd ${W_DIR}/viz
vars=("shared_PC1" "lang_specific" "iq_specific")
for var in ${vars[@]}; do
    STAT_MAP=${W_DIR}/tstat_$var.nii.gz

    # Left hemisphere
    3dVol2Surf \
        -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
        -surf_A lh.smoothwm \
        -surf_B lh.pial \
        -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
        -grid_parent ${STAT_MAP} \
        -map_func ave \
        -out_1D ${var}.lh.1D.dset

    # Right hemisphere
    3dVol2Surf \
        -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
        -surf_A rh.smoothwm \
        -surf_B rh.pial \
        -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii \
        -grid_parent ${STAT_MAP} \
        -map_func ave \
        -out_1D ${var}.rh.1D.dset
done

suma -spec $TEMP_DIR/SUMA/std.141.fsaverage_both.spec \
  -sv $TEMP_DIR/SUMA/fsaverage_SurfVol.nii

plotting.plot_glass_brain("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/tstat_shared_PC1.nii.gz",colorbar=True, display_mode='ortho')


## wb
W_DIR=/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF
map=${W_DIR}/tstat_shared_PC1.nii.gz
# Convert NIFTI to CIFTI for surface mapping
wb_command -volume-to-surface-mapping \
    ${map} \
    /path/to/freesurfer/subjects/fsaverage/surf/lh.midthickness.32k_fs_LR.surf.gii \
    your_stat_lh.func.gii \
    -trilinear

wb_command -volume-to-surface-mapping \
    your_stat.nii.gz \
    /path/to/freesurfer/subjects/fsaverage/surf/rh.midthickness.32k_fs_LR.surf.gii \
    your_stat_rh.func.gii \
    -trilinear

# Create scene file for multiple views
cat > visualize.scene << 'EOF'
<?xml version="1.0" ?>
<Scene>
  <MetaData name="Description" value="Your Statistic Map"/>
  <SurfaceDisplaySettings>
    <Surface identifier="lh.midthickness.32k_fs_LR">
      <Overlay identifier="your_stat_lh" display="ON"/>
    </Surface>
    <Surface identifier="rh.midthickness.32k_fs_LR">
      <Overlay identifier="your_stat_rh" display="ON"/>
    </Surface>
  </SurfaceDisplaySettings>
</Scene>
EOF

# Launch Workbench with all views
wb_view \
    -surface /path/to/freesurfer/subjects/fsaverage/surf/lh.midthickness.32k_fs_LR.surf.gii \
    -metric your_stat_lh.func.gii \
    -surface /path/to/freesurfer/subjects/fsaverage/surf/rh.midthickness.32k_fs_LR.surf.gii \
    -metric your_stat_rh.func.gii \
    -scene visualize.scene \
    -show-scene
