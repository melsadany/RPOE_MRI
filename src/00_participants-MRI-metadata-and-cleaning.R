################################################################################
#                    building metadata and cleaning MRI data                   #
################################################################################
# this is to identify what kind of data is collected per participant
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(oro.nifti, lib.loc = sub("tximpute", "ENA", lib.location))
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri"
setwd(project.dir)
################################################################################
################################################################################
# make metadata of what's collected per ID
# meta for participant ids
meta <- data.frame(sub = list.dirs("/Dedicated/jmichaelson-sdata/MRI/RPOE/RPOE_MR/rawdata/", 
                                   recursive = F, full.names = F)) %>%
  mutate(te_id = sub("sub-", "", sub),
         raw_dir = paste0("/Dedicated/jmichaelson-sdata/MRI/RPOE/RPOE_MR/rawdata/", sub))
# filter the meta to keep participants of interest that you're adding to the list
meta <- meta[c(78:80),]
# get the file path for the LAST collected T1w image
t1 <- data.frame(file = list.files(meta$raw_dir, pattern = "T1.*\\.nii", recursive = T)) %>%
  filter(!grepl("T1a", file)) %>% # drop the weird second nifti files. probably failed dicom conversion from Zeru
  mutate(te_id = sub(".*sub-", "", sub("_MPRAGE.*", "", file)),
         ses = as.numeric(str_sub(sub("ses-", "", sub("/anat.*", "", file)), start = 1, end = 8)),
         te_id = sub("_ses.*", "", te_id)) %>%
  filter(grepl("anat", file)) %>% 
  group_by(te_id) %>% slice_max(order_by = ses, n = 1, with_ties = F) %>% ungroup() %>%
  select(te_id, t1_file = file) %>%
  mutate(te_id = ifelse(te_id == "5134", "2E_066", te_id))
# get the file path for the LAST collected T2w image
t2 <- data.frame(file = list.files(meta$raw_dir, pattern = "T2.*\\.nii", recursive = T)) %>%
  filter(!grepl("T2a", file)) %>% # drop the weird second nifti files. probably failed dicom conversion from Zeru
  mutate(te_id = sub(".*sub-", "", sub("_Sag.*", "", file)),
         ses = as.numeric(str_sub(sub("ses-", "", sub("/anat.*", "", file)), start = 1, end = 8)),
         te_id = sub("_ses.*", "", te_id)) %>%
  filter(grepl("anat", file)) %>% 
  group_by(te_id) %>% slice_max(order_by = ses, n = 1, with_ties = F) %>% ungroup() %>%
  select(te_id, t2_file = file) %>%
  mutate(te_id = ifelse(te_id == "5134", "2E_066", te_id))
# get the file path for the LAST collected PS-VC image
# old files from Zeru have the pattern PS-VC
# new files that I converted have the pattern PS_VC 
psvc <- data.frame(file = list.files(meta$raw_dir, pattern = "PS[-|_]VC", recursive = T)) %>%
  mutate(te_id = sub(".*sub-", "", sub("_PS-VC.*", "", file)),
         ses = as.numeric(str_sub(sub("ses-", "", sub("/func.*", "", sub("/other.*", "", file))), start = 1, end = 8)),
         te_id = sub("_ses.*", "", te_id),
         full_path = paste0("/Dedicated/jmichaelson-sdata/MRI/RPOE/RPOE_MR/rawdata/sub-", te_id, "/", file),
         exsits = ifelse(file.exists(full_path), T, F),
         dim4 = sapply(full_path, function(f) {
           ifelse(grepl("nii", f), 
                  nifti_header(f)@dim_[5], NA)})) 
psvc2 <- psvc %>% filter(dim4 > 1200) %>%
  distinct(te_id, ses, .keep_all = T) %>% 
  select(te_id, psvc_file = file) %>%
  mutate(te_id = ifelse(te_id == "5134", "2E_066", te_id))
# get the file path for the LAST collected resting state run1 image
rest1 <- data.frame(file = list.files(meta$raw_dir, pattern = "REST_Run_1.*\\.nii", recursive = T)) %>%
  mutate(te_id = sub(".*sub-", "", sub("_fMRI.*", "", file)),
         ses = as.numeric(str_sub(sub("ses-", "", sub("/func.*", "", sub("/other.*", "", file))), start = 1, end = 8)),
         te_id = sub("_ses.*", "", te_id),
         full_path = paste0("/Dedicated/jmichaelson-sdata/MRI/RPOE/RPOE_MR/rawdata/sub-", te_id, "/", file),
         exsits = ifelse(file.exists(full_path), T, F),
         dim4 = sapply(full_path, function(f) {
           ifelse(grepl("nii", f), 
                  nifti_header(f)@dim_[5], NA)})) 
rest12 <- rest1 %>% filter(dim4 > 100) %>%
  group_by(te_id) %>% slice_max(order_by = ses, n = 1) %>% ungroup() %>%
  distinct(te_id, ses, .keep_all = T) %>% 
  select(te_id, rest_run1_file = file) %>%
  mutate(te_id = ifelse(te_id == "5134", "2E_066", te_id))
# get the file path for the LAST collected resting state run2 image
rest2 <- data.frame(file = list.files(meta$raw_dir, pattern = "REST_Run_2.*\\.nii", recursive = T)) %>%
  mutate(te_id = sub(".*sub-", "", sub("_fMRI.*", "", file)),
         ses = as.numeric(str_sub(sub("ses-", "", sub("/func.*", "", sub("/other.*", "", file))), start = 1, end = 8)),
         te_id = sub("_ses.*", "", te_id),
         full_path = paste0("/Dedicated/jmichaelson-sdata/MRI/RPOE/RPOE_MR/rawdata/sub-", te_id, "/", file),
         exsits = ifelse(file.exists(full_path), T, F),
         dim4 = sapply(full_path, function(f) {
           ifelse(grepl("nii", f), 
                  nifti_header(f)@dim_[5], NA)})) 
rest22 <- rest2 %>% filter(dim4 > 100) %>%
  group_by(te_id) %>% slice_max(order_by = ses, n = 1) %>% ungroup() %>%
  distinct(te_id, ses, .keep_all = T) %>% 
  select(te_id, rest_run2_file = file) %>%
  mutate(te_id = ifelse(te_id == "5134", "2E_066", te_id))

meta.all.n <- full_join(meta, 
                        t1 %>% drop_na(t1_file)) %>% 
  full_join(t2) %>% 
  full_join(psvc2) %>%
  full_join(rest12) %>%
  full_join(rest22) %>%
  filter(!grepl("test", te_id),
         !grepl("-", te_id)) %>%
  distinct()
# read old meta
meta.all.o <- read_tsv("data/participants-metadata.tsv")

meta.all.n2 <- full_join(meta.all.o, 
                         meta.all.n) # only choose new records here 
write_tsv(meta.all.n2, "data/participants-metadata.tsv")
################################################################################
################################################################################
# copy these images to my folder, with consistent naming
to <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/raw/clean"
# keep new participants to make the symlink only
meta.all <- meta.all.n
for (i in 1:nrow(meta.all)) {
  id <- meta.all$te_id[i]
  # delete old files
  system(paste0("rm -rf ", to, "/", id))
  # make folder
  system(paste0("mkdir -p ", to, "/", id))
  # symlink T1
  ifelse(is.na(meta.all$t1_file[i]), print(paste0(id, ": NO T1")),
         system(paste0("ln -s ", meta.all$raw_dir[i], "/", meta.all$t1_file[i],
                       " ", to, "/", id, "/sub-", id, "_anat-T1.nii.gz", ";",
                       "ln -s ", meta.all$raw_dir[i], "/", 
                       sub("nii\\.gz", "json", meta.all$t1_file[i]),
                       " ", to, "/", id, "/sub-", id, "_anat-T1.json")))
  # symlink T2
  ifelse(is.na(meta.all$t2_file[i]), print(paste0(id, ": NO T2")),
         system(paste0("ln -s ", meta.all$raw_dir[i], "/", meta.all$t2_file[i],
                       " ", to, "/", id, "/sub-", id, "_anat-T2.nii.gz", ";",
                       "ln -s ", meta.all$raw_dir[i], "/", 
                       sub("nii\\.gz", "json", meta.all$t2_file[i]),
                       " ", to, "/", id, "/sub-", id, "_anat-T2.json")))
  # symlink psvc
  ifelse(is.na(meta.all$psvc_file[i]), print(paste0(id, ": NO PSVC")),
         system(paste0("ln -s ", meta.all$raw_dir[i], "/", meta.all$psvc_file[i],
                       " ", to, "/", id, "/sub-", id, "_fMRI-PSVC.nii.gz", ";",
                       "ln -s ", meta.all$raw_dir[i], "/", 
                       sub("nii\\.gz", "json", meta.all$psvc_file[i]),
                       " ", to, "/", id, "/sub-", id, "_fMRI-PSVC.json")))
  # symlink resting run 1
  ifelse(is.na(meta.all$rest_run1_file[i]), print(paste0(id, ": NO REST1")),
         system(paste0("ln -s ", meta.all$raw_dir[i], "/", meta.all$rest_run1_file[i],
                       " ", to, "/", id, "/sub-", id, "_fMRI-REST1.nii.gz", ";",
                       "ln -s ", meta.all$raw_dir[i], "/", 
                       sub("nii\\.gz", "json", meta.all$rest_run1_file[i]),
                       " ", to, "/", id, "/sub-", id, "_fMRI-REST1.json")))
  # symlink resting run 2
  ifelse(is.na(meta.all$rest_run2_file[i]), print(paste0(id, ": NO REST2")),
         system(paste0("ln -s ", meta.all$raw_dir[i], "/", meta.all$rest_run2_file[i],
                       " ", to, "/", id, "/sub-", id, "_fMRI-REST2.nii.gz", ";",
                       "ln -s ", meta.all$raw_dir[i], "/", 
                       sub("nii\\.gz", "json", meta.all$rest_run2_file[i]),
                       " ", to, "/", id, "/sub-", id, "_fMRI-REST2.json")))
}
# write a list of T1 available participants only
write_tsv(meta.all.n2 %>% 
            drop_na(t1_file) %>% 
            filter(te_id %in% paste0("2E_0", c(23, 31:105))) %>% # only RPOE 
            mutate(t1_file = paste0(to, "/", te_id, "/sub-", te_id, "_anat-T1.nii.gz")) %>% 
            select(te_id, t1_file),
          "data/participants-metadata_T1-only.tsv", col_names = F)
# write a list of T2 available participants only
write_tsv(meta.all.n2 %>% 
            drop_na(t2_file) %>% 
            filter(te_id %in% paste0("2E_0", c(23, 31:105))) %>% # only RPOE
            mutate(t2_file = paste0(to, "/", te_id, "/sub-", te_id, "_anat-T2.nii.gz")) %>% 
            select(te_id, t2_file),
          "data/participants-metadata_T2-only-RPOE.tsv", col_names = F)
# write a list of PSVC available participants only
write_tsv(meta.all.n2 %>% 
            drop_na(psvc_file) %>% 
            mutate(psvc_file = paste0(to, "/", te_id, "/sub-", te_id, "_fMRI-PSVC.nii.gz")) %>% 
            select(te_id, psvc_file),
          "data/participants-metadata_psvc-only.tsv", col_names = F)
# write a list of REST1 available participants only
write_tsv(meta.all.n2 %>% 
            drop_na(rest_run1_file) %>% 
            filter(te_id %in% paste0("2E_0", c(23, 31:105))) %>% # only RPOE
            mutate(rest_run1_file = paste0(to, "/", te_id, "/sub-", te_id, "_fMRI-REST1.nii.gz")) %>% 
            select(te_id, rest_run1_file),
          "data/participants-metadata_rest1-only-RPOE.tsv", col_names = F)
# write a list of REST2 available participants only
write_tsv(meta.all.n2 %>% 
            drop_na(rest_run2_file) %>% 
            filter(te_id %in% paste0("2E_0", c(23, 31:105))) %>% # only RPOE
            mutate(rest_run2_file = paste0(to, "/", te_id, "/sub-", te_id, "_fMRI-REST2.nii.gz")) %>% 
            select(te_id, rest_run2_file),
          "data/participants-metadata_rest2-only-RPOE.tsv", col_names = F)

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
