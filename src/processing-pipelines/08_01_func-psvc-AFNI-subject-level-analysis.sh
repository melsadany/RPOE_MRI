module load FSL
AFNIDIR=/Dedicated/jmichaelson-wdata/msmuhammad/workbench/afni/22.3.07
export PATH=${AFNIDIR}:${ANTSDIR}:${PATH}
conda activate ENA

participant_id=$1

###### trial data
<< "trial"
participant_id=2E_126
trial
######

echo ${participant_id}
P_DIR=/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri
sub_DIR=${P_DIR}/data/derivatives/func/USE_THIS/${participant_id}/PS-VC
MNI_BRAIN_M=/Dedicated/jmichaelson-wdata/msmuhammad/refs/fsl-data/standard/MNI152_T1_2mm_brain_mask.nii.gz
task_MNI=${sub_DIR}/B_sub-${participant_id}_PS-VC_MNI-reg.nii.gz

rm -rf ${sub_DIR}/sla
mkdir -p ${sub_DIR}/sla
cd ${sub_DIR}/sla
    
STIM_FOLDER=${P_DIR}/data/derivatives/PSVC-task-TS/del5/afni
################################################################################
################################################################################
################################################################################
## drop the first 5 volumes
task_MNI_del=${sub_DIR}/B_sub-${participant_id}_PS-VC_MNI-reg_del5.nii.gz
3dTcat -prefix ${task_MNI_del} "${task_MNI}[5..$]"



<< "normalize"
## get the mean activity of PSVC, then calculate the Percent Signal Change
##    you basically divide the signal by the average
3dTstat -mean -prefix ${sub_DIR}/B_sub-${participant_id}_PS-VC_MNI-reg_mean.nii.gz \
    ${task_MNI}
task_scaled=${sub_DIR}/E_sub-${participant_id}_PS-VC_MNI-reg_scaled-PSC.nii.gz
3dcalc -a ${sub_DIR}/B_sub-${participant_id}_PS-VC_MNI-reg.nii.gz \
    -b ${sub_DIR}/B_sub-${participant_id}_PS-VC_MNI-reg_mean.nii.gz \
    -expr "((a - b) / b) * 100" \
    -prefix ${task_scaled}
normalize
    


###################################
## get motion corrections params ##
###################################
raw_MC_param=${P_DIR}/data/derivatives/func/${participant_id}/run-3/PS-VC/F_sub-${participant_id}_PS-VC_null.nii.gzMOCOparams.csv
tail -n +7 ${raw_MC_param} | cut -d',' -f3-8 > MOCO-params-trimmed.1D


    
##################
## major design ##
##################
mkdir -p ${sub_DIR}/sla/major
cd ${sub_DIR}/sla/major
3dDeconvolve -input ${task_MNI_del} \
    -polort A \
    -num_stimts 7 \
    -stim_times_AM1 1 ${STIM_FOLDER}/instructions.txt "dmUBLOCK(1)" -stim_label 1 instructions \
    -stim_times 2 ${STIM_FOLDER}/baseline.txt "BLOCK(3,1)" -stim_label 2 baseline \
    -stim_times 3 ${STIM_FOLDER}/PS_samediff.txt "BLOCK(3,1)" -stim_label 3 PS_samediff \
    -stim_times 4 ${STIM_FOLDER}/RAN.txt "BLOCK(16,1)" -stim_label 4 RAN \
    -stim_times 5 ${STIM_FOLDER}/semantic_coherence.txt "BLOCK(8,1)" -stim_label 5 semantic_coherence \
    -stim_times 6 ${STIM_FOLDER}/word_association__word.txt "BLOCK(5,1)" -stim_label 6 word_association_word \
    -stim_times 7 ${STIM_FOLDER}/word_association__number.txt "BLOCK(5,1)" -stim_label 7 word_association_number \
    -ortvec ${sub_DIR}/sla/MOCO-params-trimmed.1D motion \
    -nobucket -jobs 20 -mask ${MNI_BRAIN_M} -x1D major
    
# run 3dREMLfit file
3dREMLfit -input ${task_MNI_del} \
    -matrix major.xmat.1D \
    -mask ${MNI_BRAIN_M} \
    -Rbuck major \
    -tout

echo "done with major model design"
################################################################################
#################    
# minor design ##
#################
mkdir -p ${sub_DIR}/sla/minor
cd ${sub_DIR}/sla/minor

3dDeconvolve -input ${task_MNI_del} \
    -polort A \
    -num_stimts 9 \
    -stim_times 1 ${STIM_FOLDER}/baseline.txt "BLOCK(3,1)" -stim_label 1 baseline \
    -stim_times 2 ${STIM_FOLDER}/PS_samediff__face.txt "BLOCK(3,1)" -stim_label 2 PS_samediff__face \
    -stim_times 3 ${STIM_FOLDER}/PS_samediff__symbol.txt "BLOCK(3,1)" -stim_label 3 PS_samediff__symbol \
    -stim_times 4 ${STIM_FOLDER}/RAN.txt "BLOCK(16,1)" -stim_label 4 RAN \
    -stim_times 5 ${STIM_FOLDER}/semantic_coherence__coherent.txt "BLOCK(8,1)" -stim_label 5 semantic_coherence__coherent \
    -stim_times 6 ${STIM_FOLDER}/semantic_coherence__incoherent.txt "BLOCK(8,1)" -stim_label 6 semantic_coherence__incoherent \
    -stim_times 7 ${STIM_FOLDER}/word_association__word.txt "BLOCK(5,1)" -stim_label 7 word_association_word \
    -stim_times 8 ${STIM_FOLDER}/word_association__number.txt "BLOCK(5,1)" -stim_label 8 word_association_number \
    -stim_times_AM1 9 ${STIM_FOLDER}/instructions.txt "dmUBLOCK(1)" -stim_label 9 instructions \
    -ortvec ${sub_DIR}/sla/MOCO-params-trimmed.1D motion \
    -nobucket -jobs 20 -mask ${MNI_BRAIN_M} -x1D minor
    
# run 3dREMLfit file
3dREMLfit -input ${task_MNI_del} \
    -matrix minor.xmat.1D \
    -mask ${MNI_BRAIN_M} \
    -Rbuck minor \
    -tout

echo "done with minor model design"
################################################################################
##########################    
# minor-detailed design ##
##########################
mkdir -p ${sub_DIR}/sla/minor_det
cd ${sub_DIR}/sla/minor_det

3dDeconvolve -input ${task_MNI_del} \
    -polort A \
    -num_stimts 50 \
    -stim_times 1 ${STIM_FOLDER}/baseline.txt "BLOCK(3,1)" -stim_label 1 baseline \
    -stim_times 2 ${STIM_FOLDER}/PS_samediff__same.txt "BLOCK(3,1)" -stim_label 2 PS_samediff__same \
    -stim_times 3 ${STIM_FOLDER}/PS_samediff__diff.txt "BLOCK(3,1)" -stim_label 3 PS_samediff__diff \
    -stim_times 4 ${STIM_FOLDER}/RAN.txt "BLOCK(16,1)" -stim_label 4 RAN \
    -stim_times 5 ${STIM_FOLDER}/semantic_coherence__PNW.txt "BLOCK(8,1)" -stim_label 5 semantic_coherence__PNW \
    -stim_times 6 ${STIM_FOLDER}/semantic_coherence__algorithm.txt "BLOCK(8,1)" -stim_label 6 semantic_coherence__algorithm \
    -stim_times 7 ${STIM_FOLDER}/semantic_coherence__approve.txt "BLOCK(8,1)" -stim_label 7 semantic_coherence__approve \
    -stim_times 8 ${STIM_FOLDER}/semantic_coherence__arrogant.txt "BLOCK(8,1)" -stim_label 8 semantic_coherence__arrogant \
    -stim_times 9 ${STIM_FOLDER}/semantic_coherence__chlorine.txt "BLOCK(8,1)" -stim_label 9 semantic_coherence__chlorine \
    -stim_times 10 ${STIM_FOLDER}/semantic_coherence__christ.txt "BLOCK(8,1)" -stim_label 10 semantic_coherence__christ \
    -stim_times 11 ${STIM_FOLDER}/semantic_coherence__concert.txt "BLOCK(8,1)" -stim_label 11 semantic_coherence__concert \
    -stim_times 12 ${STIM_FOLDER}/semantic_coherence__concrete.txt "BLOCK(8,1)" -stim_label 12 semantic_coherence__concrete \
    -stim_times 13 ${STIM_FOLDER}/semantic_coherence__cylinder.txt "BLOCK(8,1)" -stim_label 13 semantic_coherence__cylinder \
    -stim_times 14 ${STIM_FOLDER}/semantic_coherence__felony.txt "BLOCK(8,1)" -stim_label 14 semantic_coherence__felony \
    -stim_times 15 ${STIM_FOLDER}/semantic_coherence__lol.txt "BLOCK(8,1)" -stim_label 15 semantic_coherence__lol \
    -stim_times 16 ${STIM_FOLDER}/semantic_coherence__lush.txt "BLOCK(8,1)" -stim_label 16 semantic_coherence__lush \
    -stim_times 17 ${STIM_FOLDER}/semantic_coherence__mixed.txt "BLOCK(8,1)" -stim_label 17 semantic_coherence__mixed \
    -stim_times 18 ${STIM_FOLDER}/semantic_coherence__sailing.txt "BLOCK(8,1)" -stim_label 18 semantic_coherence__sailing \
    -stim_times 19 ${STIM_FOLDER}/semantic_coherence__shopping.txt "BLOCK(8,1)" -stim_label 19 semantic_coherence__shopping \
    -stim_times 20 ${STIM_FOLDER}/semantic_coherence__shriek.txt "BLOCK(8,1)" -stim_label 20 semantic_coherence__shriek \
    -stim_times 21 ${STIM_FOLDER}/semantic_coherence__species.txt "BLOCK(8,1)" -stim_label 21 semantic_coherence__species \
    -stim_times 22 ${STIM_FOLDER}/semantic_coherence__spinach.txt "BLOCK(8,1)" -stim_label 22 semantic_coherence__spinach \
    -stim_times 23 ${STIM_FOLDER}/semantic_coherence__susan.txt "BLOCK(8,1)" -stim_label 23 semantic_coherence__susan \
    -stim_times 24 ${STIM_FOLDER}/semantic_coherence__sweater.txt "BLOCK(8,1)" -stim_label 24 semantic_coherence__sweater \
    -stim_times 25 ${STIM_FOLDER}/semantic_coherence__symptoms.txt "BLOCK(8,1)" -stim_label 25 semantic_coherence__symptoms \
    -stim_times 26 ${STIM_FOLDER}/word_association__answer.txt "BLOCK(5,1)" -stim_label 26 word_association__answer \
    -stim_times 27 ${STIM_FOLDER}/word_association__boy.txt "BLOCK(5,1)" -stim_label 27 word_association__boy \
    -stim_times 28 ${STIM_FOLDER}/word_association__brother.txt "BLOCK(5,1)" -stim_label 28 word_association__brother \
    -stim_times 29 ${STIM_FOLDER}/word_association__cross.txt "BLOCK(5,1)" -stim_label 29 word_association__cross \
    -stim_times 30 ${STIM_FOLDER}/word_association__cut.txt "BLOCK(5,1)" -stim_label 30 word_association__cut \
    -stim_times 31 ${STIM_FOLDER}/word_association__drop.txt "BLOCK(5,1)" -stim_label 31 word_association__drop \
    -stim_times 32 ${STIM_FOLDER}/word_association__eat.txt "BLOCK(5,1)" -stim_label 32 word_association__eat \
    -stim_times 33 ${STIM_FOLDER}/word_association__fast.txt "BLOCK(5,1)" -stim_label 33 word_association__fast \
    -stim_times 34 ${STIM_FOLDER}/word_association__gold.txt "BLOCK(5,1)" -stim_label 34 word_association__gold \
    -stim_times 35 ${STIM_FOLDER}/word_association__hand.txt "BLOCK(5,1)" -stim_label 35 word_association__hand \
    -stim_times 36 ${STIM_FOLDER}/word_association__heavy.txt "BLOCK(5,1)" -stim_label 36 word_association__heavy \
    -stim_times 37 ${STIM_FOLDER}/word_association__home.txt "BLOCK(5,1)" -stim_label 37 word_association__home \
    -stim_times 38 ${STIM_FOLDER}/word_association__horse.txt "BLOCK(5,1)" -stim_label 38 word_association__horse \
    -stim_times 39 ${STIM_FOLDER}/word_association__month.txt "BLOCK(5,1)" -stim_label 39 word_association__month \
    -stim_times 40 ${STIM_FOLDER}/word_association__plant.txt "BLOCK(5,1)" -stim_label 40 word_association__plant \
    -stim_times 41 ${STIM_FOLDER}/word_association__ride.txt "BLOCK(5,1)" -stim_label 41 word_association__ride \
    -stim_times 42 ${STIM_FOLDER}/word_association__spring.txt "BLOCK(5,1)" -stim_label 42 word_association__spring \
    -stim_times 43 ${STIM_FOLDER}/word_association__table.txt "BLOCK(5,1)" -stim_label 43 word_association__table \
    -stim_times 44 ${STIM_FOLDER}/word_association__turn.txt "BLOCK(5,1)" -stim_label 44 word_association__turn \
    -stim_times 45 ${STIM_FOLDER}/word_association__wall.txt "BLOCK(5,1)" -stim_label 45 word_association__wall \
    -stim_times 46 ${STIM_FOLDER}/word_association__water.txt "BLOCK(5,1)" -stim_label 46 word_association__water \
    -stim_times 47 ${STIM_FOLDER}/word_association__white.txt "BLOCK(5,1)" -stim_label 47 word_association__white \
    -stim_times 48 ${STIM_FOLDER}/word_association__wish.txt "BLOCK(5,1)" -stim_label 48 word_association__wish \
    -stim_times 49 ${STIM_FOLDER}/word_association__work.txt "BLOCK(5,1)" -stim_label 49 word_association__work \
    -stim_times_AM1 50 ${STIM_FOLDER}/instructions.txt "dmUBLOCK(1)" -stim_label 50 instructions \
    -ortvec ${sub_DIR}/sla/MOCO-params-trimmed.1D motion \
    -nobucket -jobs 20 -mask ${MNI_BRAIN_M} -x1D minor-detailed
    
# run 3dREMLfit file
3dREMLfit -input ${task_MNI_del} \
    -matrix minor-detailed.xmat.1D \
    -mask ${MNI_BRAIN_M} \
    -Rbuck minor-detailed \
    -tout

echo "done with minor-detailed model design"
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## extract nifti maps for Tstat and coefficients from each model fit

################
## major task ##
################
task_des=${sub_DIR}/sla/major
cd ${task_des}
stim_labels=("instructions" "baseline" "PS_samediff" "RAN" "semantic_coherence" "word_association__word" "word_association__number")
for i in "${!stim_labels[@]}"; do
    stim_label="${stim_labels[$i]}"
    index=$((i*2 + 1)); index2=$((index + 1))
    echo $stim_label; echo $index; echo $index2
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "major+orig[${index}]"
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "major+orig[${index2}]"
done

################
## minor task ##
################
task_des=${sub_DIR}/sla/minor
cd ${task_des}
stim_label="PS_samediff__face"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor+orig[3]"
stim_label="PS_samediff__face"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor+orig[4]"
stim_label="PS_samediff__symbol"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor+orig[5]"
stim_label="PS_samediff__symbol"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor+orig[6]"
stim_label="semantic_coherence__coherent"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor+orig[9]"
stim_label="semantic_coherence__coherent"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor+orig[10]"
stim_label="semantic_coherence__incoherent"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor+orig[11]"
stim_label="semantic_coherence__incoherent"; 3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor+orig[12]"

#########################
## minor-detailed task ##
#########################
task_des=${sub_DIR}/sla/minor_det
cd ${task_des}
stim_labels=("same" "diff")
for i in "${!stim_labels[@]}"; do
    stim_label="PS_samediff__${stim_labels[$i]}"
    index=$((i*2 + 3)); index2=$((index + 1))
    echo $stim_label; echo $index; echo $index2
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor-detailed+orig[${index}]"
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor-detailed+orig[${index2}]"
done

stim_labels=("PNW" "algorithm" "approve" "arrogant" "chlorine" "christ" "concert" "concrete" "cylinder" "felony" "lol" "lush" "mixed" "sailing" "shopping" "shriek" "species" "spinach" "susan" "sweater" "symptoms")
for i in "${!stim_labels[@]}"; do
    stim_label="semantic_coherence__${stim_labels[$i]}"
    index=$((i*2 + 9)); index2=$((index + 1))
    echo $stim_label; echo $index; echo $index2
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor-detailed+orig[${index}]"
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor-detailed+orig[${index2}]"
done

stim_labels=("answer" "boy" "brother" "cross" "cut" "drop" "eat" "fast" "gold" "hand" "heavy" "home" "horse" "month" "plant" "ride" "spring" "table" "turn" "wall" "water" "white" "wish" "work")
for i in "${!stim_labels[@]}"; do
    stim_label="word_association__${stim_labels[$i]}"
    index=$((i*2 + 51)); index2=$((index + 1))
    echo $stim_label; echo $index; echo $index2
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Coef.nii.gz" -mean "minor-detailed+orig[${index}]"
    3dTstat -prefix "${task_des}/sub-${participant_id}_${stim_label}_Tstat.nii.gz" -mean "minor-detailed+orig[${index2}]"
done



echo "DONE with coefficients and Tstat maps extraction per task"

echo "------------------------"
