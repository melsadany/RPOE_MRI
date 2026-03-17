################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/workbench/customized-functions/make-nifti-array.R"))
################################################################################
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
################################################################################
## load data
## voxel fALFF
falff.voxel <- read_rds("data/derivatives/R-func/REST-fALFF-voxel-wise.rds") %>%
  rename_at(.vars=vars(starts_with("sub")),.funs=function(x)sub("sub__","",x))
## task activity maps
pload("data/derivatives/R-func/psvc-subject-level-task-Tstats-from-afni.rda.pxz")
p.of.int <- intersect(colnames(falff.voxel)[-c(1:6)],sub(".*___","",colnames(afni.res)))
wat.afni.res <- afni.res %>% select(1:5,starts_with("word_association__word"))%>%
  rename_all(.funs = function(x) sub("word_association__word___","",x)) %>%
  select(1:5,p.of.int)
rm(afni.res);gc()

## group-level task activity
wat.group <- make_img_array("data/derivatives/func/AFNI-group-level-maps/word_association__word/word_association__word_MEMA_group-Coef.nii.gz",
                            dim = 3,labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
wat.group.fdr <- make_img_array("data/derivatives/func/AFNI-group-level-maps/word_association__word/word_association__word_MEMA_group-FDR.nii.gz",
                                dim = 3,labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
wat.group.cl <- inner_join(wat.group %>% select(starts_with("mni"),coef=intensity),
                           wat.group.fdr %>% select(starts_with("mni"),FDR=intensity))
table(wat.group.cl$FDR<0.05);rm(wat.group,wat.group.fdr)

################################################################################
################################################################################
################################################################################
################################################################################
### make a reference brain map
source("../../../workbench/customized-functions/make-nifti-array.R")
lang.network <- make_img_array("../../../data/neurosynth/language_network_association-test_z_FDR_0.01.nii.gz",
                               dim = 3,labels = "../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")
lipkin.lang.net <- make_img_array("../../../refs/lipkin-language/atlases (alternative)/langloc_n806_p<0.01_atlas.nii") %>% select(-starts_with("mni")) %>%
  inner_join(read_rds("../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds"))
lipkin.mdn.net <- make_img_array("../../../refs/lipkin-MDN/MDloc_n691_top10%_atlas.nii") %>% select(-starts_with("mni")) %>%
  inner_join(read_rds("../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds"))


ref <- read_rds("../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds") %>%
  select(starts_with("mni"), cerebr_structure, roi_name)%>% drop_na(cerebr_structure) %>%
  filter(!grepl("Cerebell|Vent|lobu",cerebr_structure,fixed = F))%>%
  left_join(lang.network%>%select(starts_with("mni"), lang_net_int=intensity)) %>%
  left_join(lipkin.lang.net%>%select(starts_with("mni"), lipkin_lang_net_int=intensity)) %>%
  left_join(lipkin.mdn.net%>%select(starts_with("mni"), lipkin_MDN_net_int=intensity)) %>%
  left_join(wat.group.cl %>% rename(task_group_coef=coef,task_group_FDR=FDR)) %>%
  mutate(roi_name=sub("17Networks_","",roi_name),
         hemisphere = case_when(grepl("LH_",roi_name)~"LH",grepl("RH_",roi_name)~"RH",T~""),
         region = sub(".H_","",roi_name),network=sub("_.*","",region),hem_network=ifelse(!is.na(network),paste0(hemisphere,"_",network),NA),
         lang_net=!is.na(lang_net_int),
         visual_area_V1 = grepl("Pericalcarine",cerebr_structure),
         hand_area = grepl("Precentral",cerebr_structure),
         contA = grepl("ContA",roi_name),contB = grepl("ContB",roi_name),contC = grepl("ContC",roi_name),
         FPN = grepl("Cont[A|B|C]",roi_name),
         task_group_sig = task_group_FDR<0.05)
################################################################################
################################################################################
################################################################################
################################################################################
registerDoMC(5)
wat.and.falff.lang.mets <- foreach(ii = 1:length(p.of.int), .combine = rbind)%dopar%{
  id.n <- p.of.int[ii]
  df <- inner_join(wat.afni.res%>%select(starts_with("mni"),task=id.n),
                   falff.voxel%>%select(starts_with("mni"),falff=id.n)) %>%
    drop_na() %>% inner_join(ref) %>% mutate(falff_raw = falff) %>%
    mutate_at(.vars=vars(c(falff,task)), .funs=function(x)scale(x)[,1])
  null.s <- sample(c(1:(2*nrow(df))), size = nrow(df))
  df2 <- df %>% 
    mutate_at(.vars = vars(c(lang_net,visual_area_V1, hand_area, contA, contB, contC, FPN)),
              .funs = function(x)sample(x))
  
  ### euclidean between task and fALFF
  # all brain
  task_to_fALFF_euclidean_ALL = dist(t(df %>% select(falff, task) %>% drop_na())) %>% as.numeric()
  # neurosynth language network
  task_to_fALFF_euclidean_language_net = dist(t(df %>% filter(lang_net) %>% select(falff, task) %>% drop_na())) %>% as.numeric()
  task_to_fALFF_euclidean_language_net_compliment = dist(t(df %>% filter(!lang_net) %>% select(falff, task) %>% drop_na())) %>% as.numeric()
  # by network
  h_nets <- unique(df$hem_network[!is.na(df$hem_network)])
  h_nets_euc <- foreach(ii=1:length(h_nets),.combine = rbind) %dopar% {
    h_net <- h_nets[ii]
    df4 <- df %>% filter(hem_network == h_net)
    data.frame(te_id = id.n, hem_network = h_net, value = dist(t(df4 %>% select(falff, task) %>% drop_na())) %>% as.numeric())
  } %>% pivot_wider(names_from = hem_network, values_from = value) %>%
    rename_at(.vars = vars(-te_id), .funs=function(x) paste0("task_to_fALFF_euclidean__",x))
  # multiply task/fALFF by MDN and lang net from Lipkin
  df5 <- df %>% select(falff_raw, lipkin_lang_net_int) %>% drop_na
  fALFF_lipkin_lang <- df5$falff_raw %*% df5$lipkin_lang_net_int; rm(df5)
  df6 <- df %>% select(falff_raw, lipkin_MDN_net_int) %>% drop_na
  fALFF_lipkin_MDN <- df6$falff_raw %*% df6$lipkin_MDN_net_int;rm(df6)
  df7 <- df %>% select(task, lipkin_lang_net_int) %>% drop_na
  task_lipkin_lang <- (abs(df7$task)/max(abs(df7$task))) %*% df7$lipkin_lang_net_int
  task_lipkin_lang.2 <- (abs(df7$task)) %*% df7$lipkin_lang_net_int
  task_lipkin_lang.3 <- (df7$task/max(df7$task)) %*% df7$lipkin_lang_net_int
  task_lipkin_lang.4 <- (df7$task) %*% df7$lipkin_lang_net_int;rm(df7)
  df8 <- df %>% select(task, lipkin_MDN_net_int) %>% drop_na
  task_lipkin_MDN <- (abs(df8$task)/max(abs(df8$task))) %*% df8$lipkin_MDN_net_int
  task_lipkin_MDN.2 <- (abs(df8$task)) %*% df8$lipkin_MDN_net_int
  task_lipkin_MDN.3 <- (df8$task/max(df8$task)) %*% df8$lipkin_MDN_net_int
  task_lipkin_MDN.4 <- (df8$task) %*% df8$lipkin_MDN_net_int;rm(df8)
  # subject to group euclidean
  task_subj_to_group_euclidean_ALL = dist(t(df %>% select(task, task_group_coef) %>% drop_na())) %>% as.numeric()
  task_subj_to_group_euclidean_task_sig = dist(t(df %>% filter(task_group_sig) %>% select(task, task_group_coef) %>% drop_na())) %>% as.numeric()
  task_subj_to_group_euclidean_task_sig_compliment = dist(t(df %>% filter(!task_group_sig) %>% select(task, task_group_coef) %>% drop_na())) %>% as.numeric()
  
  ## combine
  data.frame(te_id = id.n,
             task_to_fALFF_euclidean_ALL = task_to_fALFF_euclidean_ALL,
             task_to_fALFF_euclidean_language_net = task_to_fALFF_euclidean_language_net,
             task_to_fALFF_euclidean_language_net_compliment = task_to_fALFF_euclidean_language_net_compliment,
             raw_fALFF_lipkin_lang = fALFF_lipkin_lang,raw_fALFF_lipkin_MDN = fALFF_lipkin_MDN,
             task_subj_to_group_euclidean_ALL = task_subj_to_group_euclidean_ALL,
             task_subj_to_group_euclidean_task_sig = task_subj_to_group_euclidean_task_sig,
             task_subj_to_group_euclidean_task_sig_compliment = task_subj_to_group_euclidean_task_sig_compliment,
             abs_max_task_lipkin_lang = task_lipkin_lang,abs_max_task_lipkin_MDN = task_lipkin_MDN,
             abs_task_lipkin_lang = task_lipkin_lang.2,abs_task_lipkin_MDN = task_lipkin_MDN.2,
             max_task_lipkin_lang = task_lipkin_lang.3,max_task_lipkin_MDN = task_lipkin_MDN.3,
             task_lipkin_lang = task_lipkin_lang.4,task_lipkin_MDN = task_lipkin_MDN.4) %>%
    inner_join(h_nets_euc) 
}
################################################################################
################################################################################
################################################################################
################################################################################
## check age sex correlations
df9 <- wat.and.falff.lang.mets %>% 
  left_join(read_csv("../shared_data/data/demo-full.csv") %>% 
              mutate(sexM=as.numeric(sex=="Male")) %>% select(te_id, MRI_age, sexM))
corr.table(df9[,-1]) %>% filter(!V1 %in% colnames(wat.and.falff.lang.mets),V2 %in% colnames(wat.and.falff.lang.mets)) %>%
  ggplot(aes(V1,V2,fill=r,label=case_when(FDR<0.05~"*",pval<0.05~small.circle,T~"")))+
  geom_tile()+redblu.col.gradient.2(label="r")+my.guides+
  geom_text()+bw.theme+null_labs
ggsave2("figs/task-falff-euc-or-dotprod-to-others-and-themselves__age-sex.png",6,12)
df10 <- df9 %>% mutate_at(.vars = vars(colnames(wat.and.falff.lang.mets)[-1]),.funs=function(x){
  df00 <- cbind(df9 %>% select(MRI_age, sexM)) %>% mutate(y=x)
  z_from_lm(df00[,1],df00[,-11])
})

list("raw" = wat.and.falff.lang.mets, "as_corrected"=df10) %>%
  write_rds("data/derivatives/R-func/task-falff-euc-or-dotprod-to-others-and-themselves.rds",compress="gz")
################################################################################
################################################################################
################################################################################
################################################################################
## get task resids
good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))
source("../../../workbench/customized-functions/make-nifti-array.R")
registerDoMC(7)
wat.afni.res <- foreach(ii=1:length(good.psvc)) %dopar% {
  id.n <- good.psvc[ii]
  file.n <- paste0("data/derivatives/func/AFNI-group-level-maps/word_association__word/resid_to_lang/res_",id.n,".nii.gz")
  df <- make_img_array(file.n);colnames(df)[7]<-id.n
  return(df)
}
t3 <- wat.afni.res %>% reduce(inner_join,by=c("x","y","z","mni_x","mni_y","mni_z"))
t4 <- read_rds("../../../refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds") %>% inner_join(t3)
t4 %>% write_rds("data/derivatives/R-func/psvc-subject-level-wat-task-afni-resids.rds",compress="gz")
## check dx correlations
t5 <- rbind(t4 %>% drop_na(roi_name) %>% group_by(roi_name) %>% summarise_at(.vars = vars(starts_with("2E")), .funs = mean) %>%
              mutate(roi_name=sub("17Networks_","", roi_name)) %>%
              pivot_longer(cols = -roi_name, names_to = "te_id") %>% rename(name=roi_name),
            t4 %>% drop_na(cerebr_structure) %>% group_by(cerebr_structure) %>% summarise_at(.vars = vars(starts_with("2E")), .funs = mean) %>%
              pivot_longer(cols = -cerebr_structure, names_to = "te_id") %>% rename(name=cerebr_structure)) %>%
  pivot_wider(names_from = name, values_from = value, id_cols = te_id) %>%
  inner_join(read_csv("../shared_data/data/demo-full.csv") %>% select(te_id, ADHD_dx,ASD_dx) %>% mutate_at(.vars=vars(-c(te_id)),.funs=as.numeric))
corr.table(t5[,-1]) %>% filter(grepl("_dx",V1),!grepl("_dx",V2)) %>%
  ggplot(aes(V1,V2,fill=r,label=ifelse(FDR<0.05,"*",ifelse(pval<0.05,".",""))))+
  geom_tile()+my.guides+redblu.col.gradient.2(label="r")+
  geom_text()+null_labs+bw.theme
ggsave2("figs/AFNI-WAT-resids-to-dx.png",6,34)
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
