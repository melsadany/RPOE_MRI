################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/make-nifti-array.R"))
library(oro.nifti)
library(fmri)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))
################################################################################
################################################################################
demo <- read_csv("../shared_data/data/demo-full.csv")
m1.m2.s <- read_rds("../shared_data/data/m1m2-sex-corrected.rds")
################################################################################
################################################################################
# identify the files from 3dDeconvolve
tasks <- sub("\\.txt", "", list.files("data/derivatives/PSVC-task-TS/del5/afni/"))[-c(1:3)]

major.t <- tasks[c(5,6,30,55,45)]
minor.t.1 <- tasks[c(3,1,2,4,12,17)]
minor.t.2 <- tasks[c(31:44,46:54,56)]
minor.t.3 <- tasks[c(7:11,15:16,18:29)]
tasks.of.int <- c(major.t, minor.t.1)

afni.res.files <- data.frame(te_id = rep(good.psvc,each=length(tasks.of.int)),task = rep(tasks.of.int,length(good.psvc)), 
                             group = rep(c(rep("major",length(major.t)),rep("minor_det", 2),rep("minor", length(minor.t.1)-2)),length(good.psvc))) %>%
  mutate(file = paste0(project.dir,"/data/derivatives/func/USE_THIS/",te_id,"/PS-VC/sla/",group,"/sub-",te_id,"_",task, "_Tstat.nii.gz"),
         exists = file.exists(file))

table(afni.res.files$exists)
################################################################################
################################################################################
################################################################################
# read the task stats per participant, combine, and save
registerDoMC(cores = 20)

afni.res <- foreach(ii = 1:nrow(afni.res.files), .combine = full_join) %dopar% {
  id <- afni.res.files$te_id[ii]
  task <- afni.res.files$task[ii]
  df <- make_img_array(file = afni.res.files$file[ii], dim = 3, 
                       labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
  df2 <- df %>% 
    select(starts_with("mni_"), roi_name,cerebr_structure, intensity)
  colnames(df2)[6] <- paste0(task, "___", id)
  return(df2)
}
gc()
psave(afni.res, file = "data/derivatives/R-func/psvc-subject-level-task-Tstats-from-afni.rda")
pload("data/derivatives/R-func/psvc-subject-level-task-Tstats-from-afni.rda.pxz")
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
## prep for group-level analysis 

# SAVE AGE,SEX COVARIATES
demo %>% mutate(MRI_age = round(as.numeric(MRI_age)), sexM = case_when(sex == "Male" ~ 1,sex == "Female" ~ 0)) %>% 
  select(te_id, MRI_age, sexM) %>% filter(te_id %in% unique(afni.res.files$te_id)) %>%
  write_delim("data/derivatives/func/AFNI-group-level-maps/age-sex-covar.txt", delim = " ")

###############################
### prepare for AFNI 3dMEMA ###
###############################
registerDoMC(5)
foreach(tt = 1:length(tasks.of.int)) %dopar% {
  task.n <- tasks.of.int[tt]
  afni.t.dir <- paste0(project.dir, "/data/derivatives/func/AFNI-group-level-maps/", task.n)
  t.files <- afni.res.files %>% filter(task == task.n)
  
  # mkdir
  system(paste0("mkdir -p ", afni.t.dir))
  # write the files' list for the 3dMEMA commands
  afni.res.files %>%
    filter(task == task.n) %>%
    mutate(Subj = te_id, Coef_file = file,
           Tstat_file = sub("_Coef", "_Tstat", file)) %>%
    select(Subj, Coef_file, Tstat_file) %>%
    write_tsv(paste0(afni.t.dir, "/subjects-list.txt"))
}


## prep for word association 3dMEMA

# read language metrics
wa.res <- read_rds("data/derivatives/MRI-log/word-association-metrics.rds") %>% filter(!te_id %in% paste0("2E_0",c(22,40,39)))
task.n <- tasks.of.int[4]
afni.t.dir <- paste0(project.dir, "/data/derivatives/func/AFNI-group-level-maps/", task.n,"/plus-language-metrics")
system(paste0("mkdir -p ",afni.t.dir))
inner_join(wa.res,demo %>% mutate(sexM = as.numeric(sex=="Male")) %>% select(te_id,MRI_age,sexM)) %>%rename(Subj=te_id)%>%
  write_delim(paste0(afni.t.dir, "/covars.txt"), delim = " ")
inner_join(wa.res,afni.res.files %>% filter(task == task.n)) %>%
  mutate(Subj = te_id, Coef_file = file,Tstat_file = sub("_Coef", "_Tstat", file)) %>%
  select(Subj, Coef_file, Tstat_file) %>%
  write_tsv(paste0(afni.t.dir, "/subjects-list.txt"))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
