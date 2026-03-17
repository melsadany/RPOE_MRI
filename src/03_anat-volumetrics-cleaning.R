################################################################################
#                          combining volumetrics data                          #
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
################################################################################
# vols from SynthSeg
t.files <- data.frame(file = list.files("data/derivatives/anat/", 
                                        pattern = "vols", 
                                        recursive = T)) %>%
  filter(grepl("run-2", file)) %>%
  mutate(te_id = sub("\\/run-2.*", "", file),
         full_path = paste0(project.dir, "/data/derivatives/anat/", file))
t.files.2 <- t.files %>% filter(grepl("synthseg2",file))
registerDoMC(4)
t.vols.20 <- foreach(i=1:nrow(t.files.2), .combine = rbind) %dopar% {
  df <- read_csv(t.files.2$full_path[i]) %>% mutate(te_id= t.files.2$te_id[i])
  return(df)
}
t.vols.21 <- t.vols.20 %>% select(-subject) %>% relocate(te_id)
rm(t.vols.10);rm(t.vols.20)


write_rds(t.vols.21, "data/derivatives/brain-vols-synthreg-raw.rds")
# system(paste0("ln -s ", project.dir, "/data/derivatives/brain-vols-synthreg-raw.rds ", project.dir, "/../shared_data/data/."))

t.vols.21 <- t.vols.21 %>% select(-`total intracranial`) %>%
  mutate_at(.vars = vars(colnames(t.vols.21)[-c(1:2)]), .funs = function(x) as.numeric(x)/t.vols.21$`total intracranial`)
demo <- read_csv("../shared_data/data/demo-full.csv")
tt2 <- t.vols.21 %>% inner_join(demo) %>%
  drop_na(MRI_age, sex)
t.vols.21$te_id[!t.vols.21$te_id %in% tt2$te_id] # participants dropped from MRI because of missing age or sex

t.vols3 <- tt2 %>%
  mutate_at(.vars = colnames(t.vols.21)[-1],
            .funs = function(x) {
              df = tt2 %>%
                mutate(y=x, MRI_age = ifelse(is.na(MRI_age), nih_age, MRI_age)) %>%
                select(y, MRI_age, sex)
              return(z_from_lm(y = df$y, x = df[,-1]))
            }) %>%
  select(colnames(t.vols.21))

View(t.vols3)
t.vols3%>%pivot_longer(cols = -c(te_id)) %>% ggplot(aes(x = value, color=name))+geom_density(show.legend = F)
write_rds(t.vols3, "data/derivatives/brain-vols-synthreg-age-sex-tot-vol-corrected.rds")
# system(paste0("ln -s ", project.dir, "/data/derivatives/brain-vols-synthreg-age-sex-tot-vol-corrected.rds ",
#               project.dir, "/../shared_data/data/."))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# get white matter, grey matter, CSF vols
fast.files <- data.frame(te_id = list.dirs("data/derivatives/anat", 
                                           full.names = F, recursive = F)) %>%
  mutate(file = paste0(project.dir, "/data/derivatives/anat/",
                       te_id, "/run-2/F_sub-", te_id, 
                       "_anat-T1_fast-volumetric-stats.txt"),
         exists = file.exists(file))
table(fast.files$exists)
registerDoMC(4)
fast.m <- foreach(i=1:nrow(fast.files), .combine = rbind) %dopar% {
  if(!fast.files$exists[i]) {return(NULL)}
  df <- read_table(fast.files$file[i]) %>%
    select(Label, vol_in_vox = 2) %>%
    mutate(te_id= fast.files$te_id[i])
  # head(df)
  return(df)
}
fast.m2 <- fast.m %>%
  relocate(te_id) %>%
  mutate(tissue = case_when(Label == 1 ~ "CSF", Label == 2 ~ "GM", Label == 3 ~ "WM")) %>%
  select(te_id, tissue, vol_in_vox) %>%
  pivot_wider(names_from = tissue, values_from = vol_in_vox, id_cols = "te_id") %>%
  mutate(GM_WM_ratio = GM/WM)
write_rds(fast.m2, "data/derivatives/brain-tissue-vols-raw.rds")
# system(paste0("ln -s ", project.dir, "/data/derivatives/brain-tissue-vols-raw.rds ",
#               project.dir, "/../shared_data/data/."))
################################################################################
################################################################################
# get the t1w/t2w ratio metrics
t1t2.files <- data.frame(te_id = list.dirs("data/derivatives/anat", full.names = F, recursive = F)) %>%
  filter(grepl("2E",te_id)) %>%
  mutate(file = paste0(project.dir, "/data/derivatives/anat/", te_id, "/run-2/anat-metrics2/sub-", te_id, "_t1w-t2w-summary-stats_2.csv"),
         exists = file.exists(file)) 
table(t1t2.files$exists)

synthseg.roi <- read_tsv(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/SynthSeg/data/labels-table_ME.tsv"))
registerDoMC(4)
t1t2.m <- foreach(i=1:nrow(t1t2.files), .combine = rbind) %dopar% {
  if(!t1t2.files$exists[i]) {return(NULL)}
  df <- read_tsv(t1t2.files$file[i]) %>%
    pivot_longer(cols = -c(1:2)) %>%
    mutate(metric = sub("_.*","",name), labels = sub(".*_","",name),
           te_id = t1t2.files$te_id[i]) %>% select(-File,-`Sub-brick`) %>%
    pivot_wider(names_from = metric, values_from = value, id_cols = c(te_id, labels))
  return(df)
}
t1t2.m2 <- t1t2.m %>% mutate(labels = as.numeric(labels)) %>%
  left_join(synthseg.roi) %>%
  filter(!grepl("cerebral cortex", structures))

write_rds(t1t2.m2, "data/derivatives/t1t2-ratios-raw.rds")
# system(paste0("ln -s ", project.dir, "/data/derivatives/t1t2-ratios-raw.rds ",
#               project.dir, "/../shared_data/data/."))
################################################################################
################################################################################
################################################################################
# vols from DL_DiReCT
dldirect.files <- data.frame(te_id = unique(t.files$te_id)) %>%
  mutate(vols_path = paste0(project.dir, "/data/derivatives/anat/",
                            te_id, "/run-2/DL_DiReCT/result-vol.csv"),
         thick_path = paste0(project.dir, "/data/derivatives/anat/",
                             te_id, "/run-2/DL_DiReCT/result-thick.csv"),
         ex = file.exists(vols_path)) %>% filter(ex)

registerDoMC(4)
dldirect.vols <- foreach(i=1:nrow(dldirect.files), .combine = rbind) %dopar% {
  df <- read_csv(dldirect.files$vols_path[i]) %>% 
    rename(te_id= SUBJECT)
  return(df)
}
write_rds(dldirect.vols, "data/derivatives/brain-vols-dldirect-raw.rds")

# thickness from DL_DiReCT
registerDoMC(4)
dldirect.thick <- foreach(i=1:nrow(dldirect.files), .combine = rbind) %dopar% {
  df <- read_csv(dldirect.files$thick_path[i]) %>% 
    rename(te_id= SUBJECT)
  return(df)
}
dldirect.thick %>%
  pivot_longer(cols = -1) %>%
  ggplot(aes(x=value, color = name)) +
  geom_density(show.legend = F) + 
  bw.theme
write_rds(dldirect.thick, "data/derivatives/brain-thickness-dldirect-raw.rds")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## right vs. left
vols <- read_rds("data/derivatives/brain-vols-synthreg-raw.rds")
t1t2 <- read_rds("data/derivatives/t1t2-ratios-raw.rds")

colnames(vols)
vols %>% select(te_id, contains(c("left", "right", "lh-", "rh-"))) %>%
  pivot_longer(cols = contains(c("left","lh")), values_to = "left", names_to = "left_roi")%>%
  pivot_longer(cols = contains(c("right","rh")), values_to = "right", names_to = "right_roi") %>%
  filter(sub("right |ctx-rh-","",right_roi)==sub("left |ctx-lh-","",left_roi)) %>%
  mutate(roi = sub("right |ctx-rh-","",right_roi)) %>%
  ggplot(aes(right,left)) +
  geom_point(shape=1) + geom_smooth(method = "lm") + ggpubr::stat_cor() +
  facet_wrap(~roi, scales = "free") +
  bw.theme
t1t2 %>% pivot_wider(names_from = structures, values_from = NZMed, id_cols = te_id) %>%
  select(te_id, contains(c("left", "right", "lh-", "rh-"))) %>%
  pivot_longer(cols = contains(c("left","lh")), values_to = "left", names_to = "left_roi")%>%
  pivot_longer(cols = contains(c("right","rh")), values_to = "right", names_to = "right_roi") %>%
  filter(sub("right |ctx-rh-","",right_roi)==sub("left |ctx-lh-","",left_roi)) %>%
  mutate(roi = sub("right |ctx-rh-","",right_roi), ratio=right/left) %>%
  ggplot(aes(right,left)) +
  geom_point(shape=1) + geom_smooth(method = "lm") + ggpubr::stat_cor() +
  facet_wrap(~roi, scales = "free") +
  bw.theme
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
