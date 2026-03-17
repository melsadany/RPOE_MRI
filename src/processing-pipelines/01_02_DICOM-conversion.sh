#!/bin/bash
# script will convert DICOMs to nifti
# you'll need to make sure that the participant dicoms are already decompressed in the sourcedata folder

: <<'decompress'
id="2E_118"
cd /Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/sourcedata/downloaded
# unzip RPACS_E23614.zip
# folder name as ses
mkdir -p ../sub-${id}
mv ${ses} ../sub-${id}/.
ENA
decompress


# 
# Define the base directory
base_dir="/Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/sourcedata"

# Define the output directory
output_dir="/Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/rawdata"

# Ensure output directory exists
mkdir -p "$output_dir"

# Function to convert DICOM to NIfTI using dcm2niix
convert_dicom_to_nifti() {
    local dicom_dir=$1
    local output_subdir=$2
    local prefix=$3

    # Create the output subdirectory if it doesn't exist
    mkdir -p "$output_subdir"

    # Command to convert DICOM to NIfTI using dcm2niix
    dcm2niix -z y -f "${prefix}_%s" -o "$output_subdir" "$dicom_dir"
}


# Map of scan names to their respective output subdirectories
declare -A scan_map=(
    ["MPRAGE_T1"]="anat"
    ["Sag_CUBE_T2"]="anat"
    ["DTI_32_DIR"]="dwi"
    ["DTI_Rev_PE"]="dwi"
    ["ORIG__DTI_32_DIR"]="dwi"
    ["fMRI_REST_Run_1"]="func"
    ["fMRI_REST_Run_2"]="func"
    ["PS_VC"]="func"
    ["VIDEO"]="func"
    ["3_Plane_Loc"]="other"
    ["3_Plane_Loc_FGRE"]="other"
    ["__A___X"]="other"
    ["GE_HOS_FOV28"]="other"
    ["Mean_Epi__242_"]="other"
    ["Mean_Epi__121_"]="other"
)



# Traverse the directory structure
# define particpants list
participants_ls=("2E_118")

for te_id in "${participants_ls[@]}"; do
    participant_dir=${base_dir}/sub-${te_id}
    # Extract te_td from the participant_id
    #te_id=$(echo "$participant_id" | sed -n 's/^sub-\(.*\)_ses-.*/\1/p')
    #ses=$(echo "$participant_id" | sed -n 's/.*_ses-\(.*\)/\1/p')

    for session in $(ls ${participant_dir}); do
        session_dir=${participant_dir}/${session}
        scans_dir=${session_dir}/scans
        for scan_name in $(ls --color=never "${scans_dir}"); do
            scan_dir="${scans_dir}/${scan_name}/resources/DICOM/files"
            # Determine the correct output subdirectory based on scan name
            scan_key=$(echo ${scan_name} | cut -d'-' -f2-)
            scan_key=$(echo ${scan_key} | sed 's/[^a-zA-Z0-9_]//g')
            echo $scan_name
            
            output_subdir=${output_dir}/sub-${te_id}/ses-${session}/${scan_map[${scan_key}]}
            prefix=sub-${te_id}_ses-${ses}_${scan_key}

            # Convert DICOM to NIfTI
            convert_dicom_to_nifti "$scan_dir" "$output_subdir" "$prefix"
            
        done
    done
done

echo "DICOM to NIfTI conversion completed."
