import glob
import os
import nibabel as nib
import numpy as np
from nilearn import plotting
import matplotlib.pyplot as plt

# ------------------------------
# 1. Gather subject file paths
# ------------------------------
base_dir = '/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/USE_THIS/'
pattern = os.path.join(base_dir, '*', 'PS-VC', 'sla', 'major', 
                       'sub-*_word_association__word_Tstat.nii.gz')
file_list = sorted(glob.glob(pattern))
print(f'Found {len(file_list)} subject files.')


# ------------------------------
# 2. Load group t-map
# ------------------------------
group_img = nib.load("/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/word_association__word/word_association__word_MEMA_group-Coef.nii.gz")
group_data = group_img.get_fdata()


# ------------------------------
# 3. Load all subject maps into a 4D array
# ------------------------------
subject_data_list = []
for f in file_list:
    img = nib.load(f)
    data = img.get_fdata()
    subject_data_list.append(data)
# Stack along the last axis -> shape (x, y, z, n_subjects)
subject_data = np.stack(subject_data_list, axis=-1)
n_subjects = subject_data.shape[-1]
print(f'Subject data shape: {subject_data.shape}')


# ------------------------------
# 4. Compute COUNT of participants exceeding group t
# ------------------------------
# Boolean array of exceedances, sum over subjects
exceed_count = (subject_data > group_data[..., np.newaxis]).sum(axis=-1)
# proportion = (subject_data > group_data[..., np.newaxis]).mean(axis=-1)


# ------------------------------
# 5. Save count map as NIfTI
# ------------------------------
count_img = nib.Nifti1Image(exceed_count, affine=group_img.affine, header=group_img.header)
nib.save(count_img, 'participant_count_exceeding_group_t.nii.gz')
# proportion_img = nib.Nifti1Image(proportion, affine=group_img.affine, header=group_img.header)
# nib.save(proportion_img, '/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/AFNI-group-level-maps/word_association__word/proportion_exceeding_group_t.nii.gz')


# ------------------------------
# 6. Visualize and save to PDF
# ------------------------------
plotting.plot_stat_map(
    count_img,
    threshold=0,
    cmap='hot',
    title=f'Number of participants exceeding group t (max={n_subjects})',
    output_file='/wdata/msmuhammad/projects/RPOE/mri/figs/count-participants-exc-group-t.pdf'
)
plt.show()
# plotting.plot_stat_map(proportion_img, threshold=0, cmap='hot',
#                        title='Proportion of participants > group t')
# plt.show()
plotting.plot_stat_map(
    group_img,
    cmap='hot',
    title=f'Number of participants exceeding group t (max={n_subjects})',
    output_file='/wdata/msmuhammad/projects/RPOE/mri/figs/group-t.pdf'
)
