################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"), "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/")
setwd(project.dir)
################################################################################
################################################################################
## get the niftis you want
maps.meta <- data.frame(file = list.files("data/derivatives/R-func/PS-VC-group-level-task-estimates-from-R/", 
                                          pattern = "\\.nii\\.gz")) %>%
  mutate(task = sub("-map.*","",file),
         value = case_when(grepl("_Est",task)~"estimate",
                           grepl("_pval",task)~"pval",
                           grepl("_statistic",task)~"statistic"),
         task = sub("_Est|_pval|_statistic","",task))
################################################################################
################################################################################
## prepare a csv file for all maps to run
config.path <- paste0(project.dir, "data/derivatives/NCT/pval-config")

output.df <- maps.meta %>% filter(value == "pval") %>%
  mutate(file_path = paste0(project.dir, "data/derivatives/R-func/PS-VC-group-level-task-estimates-from-R/",
                            file),
         config = config.path,
         output_path = paste0(project.dir, "data/derivatives/NCT/R-glm/", 
                              task, "__", value)) %>%
  select(file_path, config, output_path, task)
write_csv(output.df, "data/derivatives/NCT/NCT-files.csv")

## do it on AFNI group-level 3dMEMEA FDR results
output.df2 <- data.frame(task=unique(maps.meta$task))%>%
  mutate(file_path=paste0(project.dir,"data/derivatives/func/AFNI-group-level-maps/",task,"/",task,"_MEMA_group-FDR.nii.gz"),
         config = paste0(project.dir,"data/derivatives/NCT/FDR-config"),
         output_path=paste0(project.dir,"data/derivatives/NCT/AFNI-3dMEMA/",task,"__FDR")) %>%
  select(colnames(output.df)) %>%
  write_csv("data/derivatives/NCT/NCT-files_AFNI.csv")

## do it on word-association-word with language metrics AFNI FDR results
output.df3 <- data.frame(var = c("group","onset","word-count","thinking-time-mean"),
                         config = paste0(project.dir,"data/derivatives/NCT/FDR-config"),
                         task="word_association__word") %>%
  mutate(file_path = paste0(project.dir,"data/derivatives/func/AFNI-group-level-maps/",task,"/plus-language-metrics/",task,"_MEMA_",var,"-Coef.nii.gz"),
         output_path = paste0(project.dir,"data/derivatives/NCT/AFNI-3dMEMA/",task,"-",var,"__FDR")) %>%
  select(colnames(output.df)) %>%
  write_csv("data/derivatives/NCT/NCT-files_AFNI-WAT.csv")
################################################################################
################################################################################
## Now, run the python script
## the python script runs on the NCT conda environment
################################################################################
################################################################################
################################################################################
################################################################################
registerDoMC(cores = 3)
nct.res <- foreach(ff = 1:nrow(output.df2), .combine = rbind) %dopar% {
  read_csv(paste0(output.df2[ff,3], "/network_correspondence.csv")) %>% select(-1) %>%
    mutate(task = output.df2[ff,4])
}

# atlases <- data.frame(group = c("MG360J12","AS400Y17","WS90_14","TY17","XS268_8","AL20","TL12","EG17"),
#                       atlas_name = c("Mattew Glasser -- 2016 360-ROI with Ji2019 12 Cole-Anticevic networks",
#                                      "Alex Schaefer-- 2018 400-ROI with Kong2021 17 networks",
#                                      "William Shirer -- 2012 90-ROI 14 networks",
#                                      "Thomas Yeo -- 17 networks",
#                                      "Xilin Shen -- 2013 268-ROI with 8 networks",
#                                      "Angela Laird -- 2011 20-node ICA maps",
#                                      "Tim Laumann -- 2015 12 networks (Power2011)",
#                                      "Evan Gordon -- 2017 17 networks"))
atlases <- data.frame(group = c("MG360J12","AS400Y17","WS90_14","TY17","XS268_8","AL20","TL12","EG17"),
                      atlas_name = c("Glasser -- 360-ROI 12 networks",
                                     "Schaefer-- 400-ROI 17 networks",
                                     "Shirer -- 90-ROI 14 networks",
                                     "Yeo -- 17 networks",
                                     "Shen -- 268-ROI 8 networks",
                                     "Laird -- 20-node ICA maps",
                                     "Laumann -- 12 networks",
                                     "Gordon -- 17 networks"),
                      id=c(1:8))


# Plot as circular lollipop chart
nct.res %>% inner_join(atlases) %>%mutate(group=factor(group,levels=atlases$group))%>%
  filter(task %in% unique(maps.meta$task)[c(5,6,9:11)]) %>%
  arrange(task,group,dice) %>%
  group_by(task) %>% 
  mutate(id = row_number(),
         angle = 90 - 360 * (id - 0.5) / max(id),
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle),
         label = ifelse(p_value < 0.05, name, "")) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(id), y = dice, color = atlas_name)) +
  geom_hline(yintercept = c(0.05,0.15,0.25,0.35,0.45), linetype=1, color ="grey", alpha=0.3)+
  geom_segment(aes(xend = as.factor(id), y = 0, yend = dice), linewidth = 0.6) +  
  geom_point(size=1) +
  geom_text(aes(x=as.factor(id), label = label, angle = angle, hjust = hjust,y = 0.45), 
            size = 2, inherit.aes = FALSE) +
  coord_polar() + ylim(-0.1, 0.6) +
  scale_color_manual(values = palette.1, name = "",guide = guide_legend(ncol =1)) +
  facet_wrap(~task, scales = "free", ncol = 1, strip.position = "left") +
  bw.theme + theme(legend.position = "bottom",panel.border = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+null_labs
ggsave2("figs/AFNI-group-NCT-major.png", 4, 20)


nct.res %>% left_join(atlases) %>%
  filter(task == "word_association__word") %>% select(-task) %>%
  arrange(group,dice) %>%
  mutate(id = row_number(),angle = 90 - 360 * (id - 0.5) / max(id),hjust = ifelse(angle < -90, 1, 0),angle = ifelse(angle < -90, angle + 180, angle),label = ifelse(p_value < 0.05, name, "")) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(id), y = dice, color = atlas_name)) +
  geom_hline(yintercept = c(0.05,0.15,0.25,0.35,0.45), linetype=1, color ="grey", alpha=0.3)+
  geom_segment(aes(xend = as.factor(id), y = 0, yend = dice), linewidth = 0.6) +  
  geom_point(size=1) +
  geom_text(aes(x=as.factor(id), label = label, angle = angle, hjust = hjust,y = 0.45), 
            size = 2.5, inherit.aes = FALSE) +
  geom_text(aes(xx,yy, label = yy),size=1.5,color="black",inherit.aes = F,
            data = data.frame(xx=1,yy=c(0.05,0.15,0.25,0.35,0.45)))+
  coord_polar() + ylim(-0.1, 0.6) +
  scale_color_manual(values = palette.1, name = "",guide = guide_legend(ncol =1)) +
  bw.theme + theme(legend.position = "right",panel.border = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+null_labs
ggsave2("figs/AFNI-group-NCT-word_association-word.pdf", 7, 5)

nct.res %>% left_join(atlases) %>%
  filter(task == "word_association__word") %>% select(-task) %>%
  arrange(group,p_value) %>% select(atlas_name, network = name, dice,p_value)%>%
  mutate_at(.vars = vars(c(dice,p_value)), .funs = function(x) signif(x,4)) %>%
  write_tsv("data/derivatives/NCT/AFNI-3dMEMA/word_association__word__FDR/network_correspondence.tsv")

################################################################################
################################################################################
################################################################################
################################################################################
## word association w language metrics
nct.res3 <- foreach(ff = 1:nrow(output.df3), .combine = rbind) %dopar% {
  read_csv(paste0(output.df3[ff,3], "/network_correspondence.csv")) %>% select(-1) %>%
    mutate(task = sub(".*_word-","",sub("__FDR","",output.df3[ff,3])))
}
nct.res3 %>% inner_join(atlases) %>%mutate(group=factor(group,levels=atlases$group))%>%
  arrange(task,group,dice) %>%
  group_by(task) %>% 
  mutate(id = row_number(),angle = 90 - 360 * (id - 0.5) / max(id),hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle),label = ifelse(p_value < 0.05, name, "")) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(id), y = dice, color = atlas_name)) +
  geom_hline(yintercept = c(0.1,0.2,0.3,0.4), linetype=1, color ="grey", alpha=0.3)+
  geom_segment(aes(xend = as.factor(id), y = 0, yend = dice), linewidth = 0.6) +  
  geom_point(size=1) +
  geom_text(aes(x=as.factor(id), label = label, angle = angle, hjust = hjust,y = 0.4), 
            size = 2, inherit.aes = FALSE) +
  coord_polar() + ylim(-0.1, 0.7) +
  scale_color_manual(values = palette.1, name = "",guide = guide_legend(ncol =1)) +
  facet_wrap(~task, scales = "free", ncol = 1, strip.position = "left") +
  bw.theme + theme(legend.position = "bottom",panel.border = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+null_labs
ggsave2("figs/AFNI-group-NCT-WAT-plus-language.png", 4, 16)
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
