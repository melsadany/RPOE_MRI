import sys
import ast
import cbig_network_correspondence as cnc
import pandas as pd
################################################################################
################################################################################
## get input data
atlas_names_list = ["MG360J12","AS400Y17","WS90_14","TY17","XS268_8","AL20","TL12","EG17"]
csv_path = sys.argv[1]
# csv_path = '/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/NCT/NCT-files.csv'
df = pd.read_csv(csv_path)

################################################################################
################################################################################
################################################################################
# Loop through all rows and process
for i, row in df.iterrows():
    print(f"Processing: {row['task']}")
    file_path = row["file_path"]
    config = row["config"]
    output_path = row["output_path"]
    
    ref_params = cnc.compute_overlap_with_atlases.DataParams(config, file_path)
    cnc.compute_overlap_with_atlases.network_correspondence(ref_params, atlas_names_list, output_path)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# #### single-use
# file_path = sys.argv[1]
# config = sys.argv[2]
# atlas_names_list = ["MG360J12","AS400Y17","WS90_14","TY17","XS268_8","AL20","TL12","EG17"]
# output_path = sys.argv[3]
# 
# ## trial
# file_path = '/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/semantic_coherence__felony/semantic_coherence__felony_MEMA_group-FDR.nii.gz'
# config = '/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/NCT/base-config'
# 
# # construct DataParams object based on the data file path and config
# ref_params = cnc.compute_overlap_with_atlases.DataParams(config, file_path)
# # compute the overlap with atlases and save the results
# cnc.compute_overlap_with_atlases.network_correspondence(ref_params, atlas_names_list, output_path)
################################################################################
################################################################################
