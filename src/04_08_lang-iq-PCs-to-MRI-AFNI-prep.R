################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
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
demo <- read_csv("../shared_data/data/demo-full.csv")
good.r1 <- c(paste0("2E_00", c(1:6,8,9)),
             paste0("2E_0", c(10:13,16:20,22,24,25,27:30,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
             paste0("2E_", c(100,102,104:106,108,109,112,115,124,126,118,131,133,134)))
good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))
approach <- tribble(~name, ~lang_input, ~cog_input, ~cog_source,
                    "features_lang_IQ", "features", "features", "IQ",
                    "features_lang_cog", "features", "features", "IQ-and-NIH",
                    "PCs_lang_IQ", "PCs", "PCs", "IQ",
                    "PCs_lang_cog", "PCs", "PCs", "IQ-and-NIH")
approach.index <- 1; approach.name <- paste0(approach$name[approach.index])
fig.dir <- paste0(project.dir,"/../shared_data/figs/0126/rgcca_results/",approach.name,"/")
data.dir <- paste0(project.dir,"/../shared_data/data/derivatives/rgcca-results/",approach.name,"/")

results <- read_rds(paste0(data.dir,"lang-iq-shared-and-specific-components.rds"))
################################################################################
################################################################################
lang.iq.components <- (results$components)[,-c(2,3)] %>%
  rename_at(.vars = vars(contains("lang_specific")),.funs=function(x)sub("lang_specific","language-specific",x))%>%
  rename_at(.vars = vars(contains("iq_specific")),.funs=function(x)sub("iq_specific","cognition-specific",x))%>%
  rename("shared language-cognition CC1"=lang_iq_CC1_sum)
################################################################################
################################################################################
################################################################################
falff.t <- lang.iq.components[,c(1:4,11)] %>% filter(te_id %in% good.r1) %>% rename(Subj=te_id) %>% 
  rename_all(.funs = function(x) str_replace_all(x,"\\.","_")) %>%
  left_join(demo %>% mutate(Subj=te_id,sex=as.numeric(sex=="Male"))%>%select(Subj,sex,MRI_age))%>%
  mutate(InputFile = paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/",
                            Subj,"/run-3/REST1/L_sub-",Subj,"_REST1_MNI-reg-fALFF.nii.gz"))%>%
  mutate_at(.vars = vars(-c(Subj,sex,InputFile)),.funs=function(x)scale(x)[,1])%>%
  mutate(masked_file = paste0("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/masked/",Subj,"_masked.nii.gz"),
         final_file = paste0(project.dir,"/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/final/",Subj,".nii.gz"))
  
falff.t %>% select(1:7, InputFile=final_file) %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/AFNI-table.txt",sep="\t",row.names = F,quote = F)
falff.t$Subj %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/subj-list.txt",sep="\t",row.names = F,quote = F,col.names = F)
falff.t[,c(5,6)] %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/age-sex.txt",sep="\t",row.names = F,quote = F)
falff.t$InputFile %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/raw-files.txt",sep="\t",row.names = F,quote = F,col.names = F)
falff.t$masked_file %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/masked-files.txt",sep="\t",row.names = F,quote = F,col.names = F)
falff.t$final_file %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/fALFF/final-files.txt",sep="\t",row.names = F,quote = F,col.names = F)


reho.t <- lang.iq.components[,c(1,4:6)] %>%
  filter(te_id %in% good.r1) %>% rename(Subj=te_id) %>% rename_all(.funs = function(x) str_replace_all(x,"\\.","_")) %>%
  left_join(demo %>% mutate(Subj=te_id,sex=as.numeric(sex=="Male"))%>%select(Subj,sex,MRI_age))%>%
  mutate(InputFile = paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/",
                            Subj,"/run-3/REST1/L_sub-",Subj,"_REST1_MNI-reg-ReHo.nii.gz"))%>%
  mutate_at(.vars = vars(-c(Subj,sex,InputFile)),.funs=function(x)scale(x)[,1])%>%
  mutate(masked_file = paste0("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/masked/",Subj,"_masked.nii.gz"),
         final_file = paste0(project.dir,"/data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/final/",Subj,".nii.gz"))

reho.t %>% select(1:6, InputFile=final_file) %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/AFNI-table.txt",sep="\t",row.names = F,quote = F)
reho.t$Subj %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/subj-list.txt",sep="\t",row.names = F,quote = F,col.names = F)
reho.t[,c(5,6)] %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/age-sex.txt",sep="\t",row.names = F,quote = F)
reho.t$InputFile %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/raw-files.txt",sep="\t",row.names = F,quote = F,col.names = F)
reho.t$masked_file %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/masked-files.txt",sep="\t",row.names = F,quote = F,col.names = F)
reho.t$final_file %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/ReHo/final-files.txt",sep="\t",row.names = F,quote = F,col.names = F)


## WAT resids

wat.res <- read_rds("../shared_data/data/derivatives/lang-iq-shared-and-specific.rds")$components[,c(1,4:6)] %>%
  left_join(demo %>% mutate(sex=as.numeric(sex=="Male"))%>%select(te_id,sex,MRI_age))%>%
  rename(Subj=te_id) %>% filter(Subj %in% good.psvc) %>% 
  mutate(InputFile=paste0(project.dir,"/data/derivatives/func/AFNI-group-level-maps/word_association__word/resid_to_lang/res_",
                                          Subj,".nii.gz")) %>% mutate_at(.vars = vars(c(shared_PC1, lang_specific, iq_specific, MRI_age)),
                                                                         .funs=scale)
wat.res %>% write.table("data/derivatives/func/AFNI-group-level-maps/lang-iq-PCs/task-resids/table.txt",sep="\t",row.names = F,quote = F)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# ## testing here instead
# falff.raw <- read_rds("data/derivatives/R-func/REST-fALFF-voxel-wise.rds")
# falff.clean <- falff.raw %>% mutate(vox=paste(mni_x,mni_y,mni_z,sep="_"))%>%
#   relocate(vox)%>%select(-colnames(falff.raw)[1:6])%>%column_to_rownames("vox")%>%
#   rename_all(.funs =function(x) sub("sub__","",x)) %>% select(lang.iq.components$te_id) %>%
#   drop_na()
# 
# ## fit a model with the 3 PCs, and get the resids
# source("/wdata/msmuhammad/workbench/customized-functions/ruv-glmnet.R")
# registerDoMC(cores = 30)
# falff.ruv <- foreach(ii = 1:nrow(falff.clean)) %dopar% {
#   vv <- rownames(falff.clean)[ii]
#   vox <- t(falff.clean[ii,])
#   df <- data.frame(vox) %>% rownames_to_column("te_id") %>% rename(val=2) %>%
#     inner_join(lang.iq.components[,c(1,4:6)]) %>% 
#     column_to_rownames("te_id") %>% mutate_all(.funs = function(x)scale(x)[,1])%>%
#     as.matrix()
#   list("resid"=as.data.frame(ruv_glmnet(y = df[,1], x=df[,-1])[,1])%>%rename(resid=1))
# }
# falff.ruv %>% write_rds("data/derivatives/R-func/REST-fALFF-voxel-wise-to-PCs-RUV-resid.rds",compress = "gz")
# falff.ruv <- read_rds("data/derivatives/R-func/REST-fALFF-voxel-wise-to-PCs-RUV-resid.rds")
# 
# ## get a PC from the resids
# falff.resid <- do.call(cbind,lapply(falff.ruv, function(z) z$resid))
# falff.resid.pca <- prcomp(scale(falff.resid))
# plot(falff.resid.pca)
# falff.resid.pc1 <- as.data.frame(scale(falff.resid.pca$x[,1]))
# 
# ## re-fit the model with the 3PCs + RUV-PC1
# falff.ruv.corrected <- foreach(ii = 1:nrow(falff.clean),.combine = rbind) %dopar% {
#   vv <- rownames(falff.clean)[ii]
#   vox <- t(falff.clean[ii,])
#   df <- data.frame(vox) %>% rownames_to_column("te_id") %>% rename(val=2) %>%
#     inner_join(lang.iq.components[,c(1,4:6)]) %>% 
#     column_to_rownames("te_id") %>% mutate(RUV_PC1=falff.resid.pc1)%>%
#     mutate_all(.funs = function(x)scale(x)[,1])
#   v.model <- glm(df[,1]~., data=df[,-1])
#   coefs_table(v.model) %>% filter(!grepl("Intercept",x)) %>%
#     rename("t"=`t val.`) %>% pivot_longer(cols=c(t,pval))%>%
#     mutate(nn=paste0(x,"__",name),vox=vv) %>%
#     pivot_wider(names_from = nn,values_from = value,id_cols = vox)
# } %>% separate(vox, into=c("mni_x","mni_y","mni_z"),sep="_")
# falff.ruv.corrected %>% write_rds("data/derivatives/R-func/REST-fALFF-voxel-wise-to-PCs-RUV-corrected-stats.rds",compress = "gz")
# falff.ruv.corrected <- read_rds("data/derivatives/R-func/REST-fALFF-voxel-wise-to-PCs-RUV-corrected-stats.rds")

################################################################################
################################################################################
## convert stats to nifti maps
# source("/wdata/msmuhammad/workbench/customized-functions/make-array-to-nifti.R")
# ref <- read_rds("/wdata/msmuhammad/refs/labeled-MNI/2mm/voxels-w-labels-R_V2.rds")[,1:6]
# ## shared
# 
# vars <- c("shared_PC1", "lang_specific","iq_specific")
# 
# for(var in vars){
#   df <- falff.ruv.corrected%>%select(starts_with("mni"),paste0(var,"__t"),paste0(var,"__pval"))%>%
#     mutate_all(.funs=as.numeric)%>%
#     inner_join(ref)%>%select(x,y,z,t=4,pval=5)
#   make_nifti_from_df(df%>%select(x,y,z,intensity=t)) %>%
#     writeNIfTI(paste0("data/derivatives/R-func/shared-lang-iq-PCs-to-fALFF-RUV/",var,"_tstat"))
#   make_nifti_from_df(df%>%select(x,y,z,intensity=pval)) %>%
#     writeNIfTI(paste0("data/derivatives/R-func/shared-lang-iq-PCs-to-fALFF-RUV/",var,"_pval"))
#   make_nifti_from_df(df%>%filter(pval<0.005)%>%select(x,y,z,intensity=t)) %>%
#     writeNIfTI(paste0("data/derivatives/R-func/shared-lang-iq-PCs-to-fALFF-RUV/",var,"_tstat_pval-01"))
# }

################################################################################
################################################################################
################################################################################
################################################################################

