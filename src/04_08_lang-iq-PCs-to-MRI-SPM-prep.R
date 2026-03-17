################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
showtext::showtext_auto()
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
approach.index <- 2; approach.name <- paste0(approach$name[approach.index])
fig.dir <- paste0(project.dir,"/../shared_data/figs/0326/rgcca_results/",approach.name,"/")
data.dir <- paste0(project.dir,"/../shared_data/data/derivatives/rgcca-results/",approach.name,"/")

results <- read_rds(paste0(data.dir,"lang-iq-shared-and-specific-components.rds"))
################################################################################
################################################################################
lang.iq.components <- (results$components)[,-c(2,3)] %>%
  rename_at(.vars = vars(contains("lang_specific")),.funs=function(x)sub("lang_specific","language-specific",x))%>%
  rename_at(.vars = vars(contains("iq_specific")),.funs=function(x)sub("iq_specific","cognition-specific",x))%>%
  rename("shared language-cognition CC1"=lang_iq_CC1_sum)
################################################################################
out.dir <- paste0(project.dir, "/data/derivatives/func/SPM/lang-cog-pcs/",approach.name)
system(paste0("mkdir -p ",out.dir,"/fALFF"))
system(paste0("mkdir -p ",out.dir,"/ReHo"))
system(paste0("mkdir -p ",out.dir,"/resids"))
################################################################################
################################################################################
## fALFF
lang.iq.components[,c(1:4,11)] %>%
  left_join(demo %>% mutate(sex=as.numeric(sex=="Male"))%>%select(te_id,sex,MRI_age))%>%
  rename(Subj=te_id) %>% filter(Subj %in% good.r1) %>% 
  mutate(InputFile=paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/",
                          Subj,"/run-3/REST1/L_sub-",Subj,"_REST1_MNI-reg-fALFF.nii.gz")) %>% 
  mutate_at(.vars = vars(c(colnames(lang.iq.components)[c(2:4,11)], MRI_age)),.funs=scale) %>% 
  write.table(paste0(out.dir,"/fALFF/spm-table.txt"),sep="\t",row.names = F,quote = F)

## ReHo
lang.iq.components[,c(1:4,11)] %>%
  left_join(demo %>% mutate(sex=as.numeric(sex=="Male"))%>%select(te_id,sex,MRI_age))%>%
  rename(Subj=te_id) %>% filter(Subj %in% good.r1) %>% 
  mutate(InputFile=paste0("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/",
                          Subj,"/run-3/REST1/K_sub-",Subj,"_REST1_MNI-reg-ReHo.nii.gz")) %>% 
  mutate_at(.vars = vars(c(colnames(lang.iq.components)[c(2:4,11)], MRI_age)),.funs=scale) %>% 
  write.table(paste0(out.dir,"/ReHo/spm-table.txt"),sep="\t",row.names = F,quote = F)

## WAT resids
lang.iq.components[,c(1:4,11)] %>%
  left_join(demo %>% mutate(sex=as.numeric(sex=="Male"))%>%select(te_id,sex,MRI_age))%>%
  rename(Subj=te_id) %>% filter(Subj %in% good.psvc) %>% 
  mutate(InputFile=paste0(project.dir,"/data/derivatives/func/AFNI-group-level-maps/word_association__word/resid_to_lang/res_",Subj,".nii.gz")) %>% 
  mutate_at(.vars = vars(c(colnames(lang.iq.components)[c(2:4,11)], MRI_age)),.funs=scale) %>% 
  write.table(paste0(out.dir,"/resids/spm-table.txt"),sep="\t",row.names = F,quote = F)
################################################################################
################################################################################
################################################################################
################################################################################
### write command to visualize SPM results



groups <- c("fALFF","ReHo","resids")
patterns <- paste0(c("shared_CC1", "lang_specific_PC1", "lang_specific_PC2", "cog_specific_PC1"), "_pos")

cmds <- c("conda activate ENA")
for (gg in 1:length(groups)) {
  group <- groups[gg]
  group.dir <- paste0(out.dir,"/", group,"/spm_results/")
  
  for (pp in 1:length(patterns)) {
    patt <- patterns[pp]
    g.p.files <- data.frame(file = list.files(group.dir, pattern = patt)) %>%
      mutate(metric = sub("_000.*","",file), full_path = paste0(group.dir,file)) %>% 
      filter(!grepl("_neg|con_|beta",file))
    
    cmd <- paste0("python /wdata/msmuhammad/workbench/customized-functions/nilearn-viz-function__final.py ",
                  "--t_img ", g.p.files$full_path[g.p.files$metric=="spmT"]," ",
                  "--p_img ", g.p.files$full_path[g.p.files$metric=="p_uncorr"]," ",
                  "--p_thresh 0.005 ",
                  "--t_thresh 0 ",
                  "--plot_mode full ",
                  "--title ",patt," ",
                  "--outdir " ,out.dir,"/",group," ",
                  "--stat_label ",patt)
    cmds <- c(cmds,cmd)
  }
}
cmds %>% write_lines(paste0(out.dir,"/viz-commands.sh"))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

