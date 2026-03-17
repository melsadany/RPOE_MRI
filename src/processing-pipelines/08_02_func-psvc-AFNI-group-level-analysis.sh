source ~/.bashrc
conda deactivate

<< "argon"
conda activate theone
AFNI=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07
export PATH=$AFNI:$PATH
argon


###### trial data

SCHF_atlas=/Dedicated/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/HCP/fslr32k/cifti/Schaefer2018_100Parcels_17Networks_order.dscalar.nii
MNI_BRAIN_M=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_2mm_brain_mask.nii.gz

# prepare a folder of symlinks for each task
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
MAP_DIR=${P_DIR}/data/derivatives/func/AFNI-group-level-maps
tasks=($(find ${MAP_DIR} -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))


<< "trial"
task="word_association__word"
trial
######

# group-level analysis
for task in "${tasks[@]}"; do
    echo ${task}
    cd ${MAP_DIR}/${task}

    ### 3dMEMA
    rm ${task}_MEMA*
    sub_list="subjects-list.txt"; set_param=""
    while read -r subj coef tstat; do
        set_param+=" ${subj} ${coef} ${tstat}"
    done < <(tail -n +2 "$sub_list")
    
    # use the local installation on topaz. the older versions of AFNI have an issue parsing parameters 
    #   and handling single/double quotes
    3dMEMA -prefix ${task}_MEMA \
        -jobs 20 \
        -set ${task} ${set_param} \
        -max_zeros 4 \
        -covariates ../age-sex-covar.txt \
        -HKtest \
        -residual_Z \
        -verb 1
    
    ## extract the beta coefficient map for the task
    3dcalc -a ${task}_MEMA+orig'[0]' -expr 'a' -prefix ${task}_MEMA_group-Coef.nii.gz
    ## FDR correction for the beta coefficient maps
    df1=$(3dAttribute BRICK_STATAUX ${task}_MEMA+orig.HEAD'[1]' | awk '{print $4}')
    3dcalc -a ${task}_MEMA+orig'[1]' -expr "fitt_t2p(a,${df1})" -prefix ${task}_MEMA_group-pval.nii.gz
    rm *FDR.nii.gz
    3dFDR -input ${task}_MEMA+orig'[1]' -qval -prefix ${task}_MEMA_group-FDR.nii.gz
    # FDR mask
    3dcalc -a ${task}_MEMA_group-FDR.nii.gz -expr 'ispositive(0.005-a)' -prefix ${task}_MEMA_group-FDR-sigmask.nii.gz
    ## clusters
    3dClusterize -inset ${task}_MEMA+orig.BRIK \
        -ithr 1 -idat 0 -bisided p=0.001 -1Dformat -NN 1 \
        -clust_nvox 8 > ${task}_MEMA_group-clusters.txt
    
done
################################################################################
################################################################################
################################################################################
## word association and metrics as covars
task="word_association__word"
echo ${task}
cd ${MAP_DIR}/${task}/plus-language-metrics

### 3dMEMA
rm ${task}_MEMA*
sub_list="subjects-list.txt"; set_param=""
while read -r subj coef tstat; do
    set_param+=" ${subj} ${coef} ${tstat}"
done < <(tail -n +2 "$sub_list")

3dMEMA -prefix ${task}_MEMA \
    -jobs 20 \
    -set ${task} ${set_param} \
    -max_zeros 4 \
    -covariates covars.txt \
    -HKtest \
    -residual_Z \
    -verb 1
    
## extract the beta coefficient map for the task, and language metrics
3dcalc -a ${task}_MEMA+orig'[0]' -expr 'a' -prefix ${task}_MEMA_group-Coef.nii.gz
rm *FDR.nii.gz
3dFDR -input ${task}_MEMA+orig'[1]' -qval -prefix ${task}_MEMA_group-FDR.nii.gz

3dcalc -a ${task}_MEMA+orig'[2]' -expr 'a' -prefix ${task}_MEMA_word-count-Coef.nii.gz
3dFDR -input ${task}_MEMA+orig'[3]' -qval -prefix ${task}_MEMA_word-count-FDR.nii.gz
3dcalc -a ${task}_MEMA+orig'[4]' -expr 'a' -prefix ${task}_MEMA_thinking-time-mean-Coef.nii.gz
3dFDR -input ${task}_MEMA+orig'[5]' -qval -prefix ${task}_MEMA_thinking-time-mean-FDR.nii.gz
3dcalc -a ${task}_MEMA+orig'[6]' -expr 'a' -prefix ${task}_MEMA_onset-Coef.nii.gz
3dFDR -input ${task}_MEMA+orig'[7]' -qval -prefix ${task}_MEMA_onset-FDR.nii.gz

################################################################################
################################################################################
################################################################################
################################################################################
mets=("fALFF" "ReHo")
for met in "${mets[@]}"; do
    echo ${met}
    cd ${MAP_DIR}/${met}
    
    
    ### 3dttest++
    ## ASD_dx
    setA_args=$(cat ASD-files.txt | awk '{print $1, $2}')
    setB_args=$(cat NASD-files.txt | awk '{print $1, $2}')
    3dttest++ \
        -prefix ${met}_ttest_asd-vs-ctrl-w-cov \
        -setA ASD $setA_args \
        -setB NASD $setB_args \
        -covariates age-sex.txt \
        -center covariates
    ## extract the statistic map for the ASD-NASD
    3dcalc -a ${met}_ttest_asd-vs-ctrl-w-cov+orig'[1]' \
       -expr 'a' -prefix ${met}_ttest_asd-vs-ctrl-w-cov-ASD-NASD-stat.nii.gz
    ## pval extraction for the statistic map
    df1=$(3dAttribute BRICK_STATAUX ${met}_ttest_asd-vs-ctrl-w-cov+orig.HEAD'[1]' | awk '{print $4}')
    3dcalc -a ${met}_ttest_asd-vs-ctrl-w-cov+orig'[1]' -expr "fitt_t2p(a,${df1})" \
        -prefix ${met}_ttest_asd-vs-ctrl-w-cov-ASD-NASD-pval.nii.gz
    3dFDR -input ${met}_ttest_asd-vs-ctrl-w-cov+orig'[1]' \
        -qval \
        -prefix ${met}_ttest_asd-vs-ctrl-w-cov-ASD-NASD-FDR.nii.gz
    
    ## pred_ASD
    setA_args=$(cat pred-asd-files.txt | awk '{print $1, $2}')
    3dttest++ \
        -prefix ${met}_ttest_pred-asd \
        -setA ASD $setA_args \
        -covariates age-sex-pred-asd.txt \
        -center covariates
    ## extract the statistic map for the ASD-NASD
    3dcalc -a ${met}_ttest_pred-asd+orig'[3]' \
       -expr 'a' -prefix ${met}_ttest_pred-asd-ASD-NASD-stat.nii.gz
    ## pval extraction for the statistic map
    df1=$(3dAttribute BRICK_STATAUX ${met}_ttest_pred-asd+orig.HEAD'[3]' | awk '{print $4}')
    3dcalc -a ${met}_ttest_pred-asd+orig'[3]' -expr "fitt_t2p(a,${df1})" \
        -prefix ${met}_ttest_pred-asd-ASD-NASD-pval.nii.gz
    3dFDR -input ${met}_ttest_pred-asd+orig'[3]' \
        -qval \
        -prefix ${met}_ttest_pred-asd-ASD-NASD-FDR.nii.gz
        
    ## pred ASD and CP
    setA_args=$(cat pred-all-files.txt | awk '{print $1, $2}')
    3dttest++ \
        -prefix ${met}_ttest_pred-all \
        -setA pred $setA_args \
        -covariates age-sex-pred-all.txt \
        -center covariates
    ## extract the statistic map for the ASD-NASD
    3dcalc -a ${met}_ttest_pred-all+orig'[3]' \
       -expr 'a' -prefix ${met}_ttest_pred-all-ASD-stat.nii.gz
    3dcalc -a ${met}_ttest_pred-all+orig'[5]' \
       -expr 'a' -prefix ${met}_ttest_pred-all-CP-stat.nii.gz
    ## pval extraction for the statistic map
    df1=$(3dAttribute BRICK_STATAUX ${met}_ttest_pred-all+orig.HEAD'[3]' | awk '{print $4}')
    3dcalc -a ${met}_ttest_pred-all+orig'[3]' -expr "fitt_t2p(a,${df1})" \
        -prefix ${met}_ttest_pred-all-ASD-pval.nii.gz
    3dFDR -input ${met}_ttest_pred-all+orig'[3]' \
        -qval \
        -prefix ${met}_ttest_pred-all-ASD-FDR.nii.gz
    
    
    


    ### 3dMVM
    rm ${met}_MVM*
    3dMVM -prefix ${met}_MVM -bsVars "ASD_dx+MRI_age+sexM" -qVars "MRI_age" -num_glt 1 -gltLabel 1 ASD_vs_NASD -gltCode 1 'ASD_dx : 1*1 -1*0' -dataTable @mvm-table.txt -jobs 30
    ## extract the statistic map for the ASD-NASD
    3dcalc -a ${met}_MVM+orig'[5]' \
       -expr 'a' -prefix ${met}_MVM-ASD-NASD-stat.nii.gz
       
    ## pval extraction for the statistic map
    df1=$(3dAttribute BRICK_STATAUX ${met}_MVM+orig.HEAD'[5]' | awk '{print $4}')
    3dcalc -a ${met}_MVM+orig'[5]' -expr "fitt_t2p(a,${df1})" \
        -prefix ${met}_MVM-ASD-NASD-pval.nii.gz
    3dFDR -input ${met}_MVM+orig'[5]' \
        -qval \
        -prefix ${met}_MVM-ASD-NASD-FDR.nii.gz
    # FDR mask
    3dcalc -a ${met}_MVM-ASD-NASD-FDR.nii.gz -expr 'ispositive(0.005-a)' \
       -prefix ${met}_MVM-ASD-NASD-FDR-sigmask.nii.gz

done

