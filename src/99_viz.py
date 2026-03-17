import os
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn import plotting
from nilearn.image import math_img
from scipy.ndimage import label
import matplotlib.pyplot as plt


project_dir = "/wdata/msmuhammad/projects/RPOE/mri"
tasks = ["PS_samediff__diff", "PS_samediff__face","PS_samediff__same","PS_samediff__symbol","PS_samediff","RAN","semantic_coherence__algorithm","semantic_coherence__approve","semantic_coherence__arrogant","semantic_coherence__chlorine","semantic_coherence__christ","semantic_coherence__coherent","semantic_coherence__concert","semantic_coherence__concrete","semantic_coherence__cylinder","semantic_coherence__felony","semantic_coherence__incoherent","semantic_coherence__lol","semantic_coherence__lush","semantic_coherence__mixed","semantic_coherence__PNW","semantic_coherence__sailing","semantic_coherence__shopping","semantic_coherence__shriek","semantic_coherence__species","semantic_coherence__spinach","semantic_coherence__susan","semantic_coherence__sweater","semantic_coherence__symptoms","semantic_coherence","word_association__answer","word_association__boy","word_association__brother","word_association__cross","word_association__cut","word_association__drop","word_association__eat","word_association__fast","word_association__gold","word_association__hand","word_association__heavy","word_association__home","word_association__horse","word_association__month","word_association__number","word_association__plant","word_association__ride","word_association__spring","word_association__table","word_association__turn","word_association__wall","word_association__water","word_association__white","word_association__wish","word_association__word","word_association__work"]

afni_group_res_files = pd.DataFrame({
    "F_file": [os.path.join(project_dir, "data/derivatives/func/AFNI-group-level-maps", task, f"{task}-intercept_F.nii.gz") for task in tasks[1:]],
    "P_file": [os.path.join(project_dir, "data/derivatives/func/AFNI-group-level-maps", task, f"{task}-intercept_pval.nii.gz") for task in tasks[1:]],
    "task": tasks[1:]
})
afni_group_res_files["F_exists"] = afni_group_res_files["F_file"].apply(os.path.exists)
afni_group_res_files["P_exists"] = afni_group_res_files["P_file"].apply(os.path.exists)

# Create output directory for plots
output_dir = os.path.join(project_dir, "figs/AFNI-group-level")
os.makedirs(output_dir, exist_ok=True)

# Loop through the DataFrame to plot existing files
for _, row in afni_group_res_files.iterrows():
    if row["F_exists"]:  # Only plot if the file exists
        output_path = os.path.join(output_dir, f"{row['task']}-intercept_F.png")
        # Plot and save
        display = plotting.plot_glass_brain(row["F_file"], threshold=20, display_mode='lzry', cmap="cold_hot", colorbar=True, alpha=0.7)
        display.savefig(output_path)
        display.close()  # Close figure to free memory
        print(f"Saved: {output_path}")




# pvalue
# Parameters
pval_threshold = 0.005  
log_pval_threshold = -np.log10(pval_threshold)
voxel_size = 2  # Voxel size in mm (assuming isotropic)
min_cluster_size = int((4 / voxel_size) ** 3)  # Minimum cluster size (4mm³ ~ at least 2 neighboring voxels)

# Loop through p-value maps
for _, row in afni_group_res_files.iterrows():
    if row["P_exists"]:  # Only process existing files
        pval_path = row["P_file"]
        output_path = os.path.join(output_dir, f"{row['task']}-intercept_pval.png")
        # Load the p-value NIfTI file
        pval_img = nib.load(pval_path)
        # Convert p-values to -log10(p-value)
        log_pval_img = math_img("-np.log10(np.clip(img, 1e-10, 1))", img=pval_img)
        # Threshold the -log10(p) image (keep significant voxels)
        thresholded_img = math_img(f"img > {log_pval_threshold}", img=log_pval_img)
        thresholded_data = thresholded_img.get_fdata()
        # Cluster labeling
        labeled_array, num_features = label(thresholded_data)
        # Remove small clusters
        cluster_sizes = np.bincount(labeled_array.flatten())[1:]  # Ignore background (index 0)
        valid_clusters = np.where(cluster_sizes >= min_cluster_size)[0] + 1  # Labels start from 1
        filtered_data = np.isin(labeled_array, valid_clusters).astype(int) * thresholded_data
        # Convert back to Nifti
        filtered_img = nib.Nifti1Image(filtered_data, affine=pval_img.affine, header=pval_img.header)
        # Plot the filtered significant clusters
        display = plotting.plot_glass_brain(filtered_img, threshold=log_pval_threshold, display_mode='lzry', cmap="magma", colorbar=True, alpha=0.7)
        display.savefig(output_path)
        display.close()
        print(f"Saved clustered p-value map: {output_path}")


for _, row in afni_group_res_files.iterrows():
    if row["P_exists"]:  # Only process existing files
        pval_path = row["P_file"]
        output_path = os.path.join(output_dir, f"{row['task']}-intercept_pval.png")
        # Load the p-value NIfTI file
        pval_img = nib.load(pval_path)
        # Compute -log10(p-value), avoiding log(0) issues
        log_pval_img = math_img("-np.log10(np.clip(img, 1e-10, 1))", img=pval_img)
        # Plot and save
        display = plotting.plot_glass_brain(log_pval_img, threshold=2, display_mode='lzry', cmap="magma", colorbar=True, alpha=0.7)
        display.savefig(output_path)
        display.close()  # Free memory
        print(f"Saved: {output_path}")




plotting.plot_glass_brain(stat_img, title="Task Activation", threshold=10, colorbar=True, display_mode='ortho')
plotting.show()

#plotting.plot_glass_brain(stat_img, threshold=3, display_mode='lzry', cmap="cold_hot", colorbar=True, alpha=0.7)










plotting.plot_glass_brain("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/2E_134/PS-VC/sla/major/sub-2E_134_word_association__word_Tstat.nii.gz",threshold=2,colorbar=True, display_mode='ortho')


project_dir = "/wdata/msmuhammad/projects/RPOE/mri"
task = "PS_samediff__face"
stat_img = os.path.join(project_dir, "data/derivatives/func/AFNI-group-level-maps", task, f"{task}_MVM-ASD-NASD-stat.nii.gz")
plotting.plot_stat_map(stat_img,threshold=2,display_mode="ortho",cmap="magma",
    cut_coords=[-38, -22, -16],title="plotting threshold=1")
plt.show()


plotting.plot_glass_brain(stat_img,threshold=2, display_mode='lzry', cmap="magma", colorbar=True, alpha=0.7)
