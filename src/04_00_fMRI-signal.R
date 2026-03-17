################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/make-nifti-array.R"))
################################################################################
# script is for extracting the BOLD signal for the fMRI
################################################################################
################################################################################
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# get a list of file names with data extracted that you can read
# participants with good reg
good.r1 <- c(paste0("2E_00", c(1:6,8,9)),
             paste0("2E_0", c(10:13,16:20,22,24,25,27:30,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
             paste0("2E_", c(100,102,104:106,108,109,112,115,124,126,118,131,133,134)))
good.r2 <- c()

good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))
################################################################################
################################################################################
################################################################################
################################################################################
### read fMRI files that are registered to MNI, reformat their data and save

## REST
registerDoMC(cores = 6)
foreach(i = c(1:length(good.r1))) %dopar% {
  id <- good.r1[i]
  base.dir <- paste0("data/derivatives/func/USE_THIS/", id, "/REST1/")
  out.file <- paste0(base.dir, "B_sub-", id, "_REST1_MNI-reg-voxel-ts-R.rds")
  if (!file.exists(out.file)) {
    df <- make_img_array(file = paste0(base.dir, "B_sub-", id, "_REST1_MNI-reg.nii.gz"),
                         dim = 4, labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
    # save the voxel intensities
    write_rds(df, out.file, compress = "gz")
    gc()
    return(NULL)
  } else {
    print(paste0("Participant: ", id, " is already processed"))
  }
}

## PSVC
foreach(i = seq(1,length(good.psvc))) %dopar% {
  id <- good.psvc[i]
  base.dir <- paste0("data/derivatives/func/USE_THIS/", id, "/PS-VC/")
  out.file <- paste0(base.dir, "B_sub-", id, "_PS-VC_MNI-reg-voxel-ts-R.rds")
  if (!file.exists(out.file)) {
    df <- make_img_array(file = paste0(base.dir, "B_sub-", id, "_PS-VC_MNI-reg.nii.gz"),
                         dim = 4, labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
    # save the voxel intensities
    write_rds(df, out.file, compress = "gz")
    gc()
    return(NULL)
  } else {
    print(paste0("Participant: ", id, " is already processed"))
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## read the AFNI calculated fALFF
foreach(i = c(1:length(good.r1))) %dopar% {
  id <- good.r1[i]
  base.dir <- paste0("data/derivatives/func/USE_THIS/", id)
  out.file <- paste0(base.dir, "/REST-voxel-fALFF.rds")
  
  if (!file.exists(out.file)) {
    df <- make_img_array(file = paste0(project.dir, "/data/derivatives/func/",
                                       id,"/run-3/REST1/L_sub-",id,"_REST1_MNI-reg-fALFF.nii.gz"),
                         dim = 3, labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds") %>%
      dplyr::select(starts_with("mni"), anat_structure, cerebr_structure, roi_name, fALFF = intensity)
    colnames(df)[7] <- paste0("sub__", id)
    # save the values intensities
    write_rds(df, out.file, compress = "gz")
    gc()
    return(NULL)
  } else {
    print(paste0("Participant: ", id, " is already processed"))
  }
}
################################################################################
################################################################################
## read the ReHo
foreach(i = c(1:length(good.r1))) %dopar% {
  id <- good.r1[i]
  base.dir <- paste0("data/derivatives/func/USE_THIS/", id)
  out.file <- paste0(base.dir, "/REST-voxel-ReHo.rds")
  
  if (!file.exists(out.file)) {
    df <- make_img_array(file = paste0(project.dir, "/data/derivatives/func/",
                                       id,"/run-3/REST1/K_sub-",id,"_REST1_MNI-reg-ReHo.nii.gz"),
                         dim = 3, labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds") %>%
      dplyr::select(starts_with("mni"), anat_structure, cerebr_structure, roi_name, intensity)
    colnames(df)[7] <- paste0("sub__", id)
    # save the values intensities
    write_rds(df, out.file, compress = "gz")
    gc()
    return(NULL)
  } else {
    print(paste0("Participant: ", id, " is already processed"))
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
##### combine fALFF and save
all.rest <- c(good.r1, good.r2)
registerDoMC(cores = 6)
all.rest.falff <- foreach(i = 1:length(all.rest), .combine = full_join) %dopar% {
  id <- all.rest[i]; file.n <- paste0("data/derivatives/func/USE_THIS/", id, "/REST-voxel-fALFF.rds")
  df <- read_rds(file.n); return(df)
}
write_rds(all.rest.falff, "data/derivatives/R-func/REST-fALFF-voxel-wise.rds", compress = "gz")
gc()
##### combine ReHo and save
registerDoMC(cores = 6)
all.rest.reho <- foreach(i = 1:length(all.rest), .combine = full_join) %dopar% {
  id <- all.rest[i]; file.n <- paste0("data/derivatives/func/USE_THIS/", id, "/REST-voxel-ReHo.rds")
  df <- read_rds(file.n); return(df)
}
write_rds(all.rest.reho, "data/derivatives/R-func/REST-ReHo-voxel-wise.rds", compress = "gz")
gc()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################