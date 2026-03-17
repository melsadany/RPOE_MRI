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
## get mitochondrial maps
mit.maps <- read_rds("/Dedicated/jmichaelson-wdata/msmuhammad/data/mitochondria-brain-maps/derivatives/all-data-2mm.rds")

## get fALFF maps
falff.maps <- read_rds("data/derivatives/R-func/REST-fALFF-voxel-wise.rds")

m1.m2 <- read_rds("../shared_data/data/m1m2.rds")
demo <- read_csv("../shared_data/data/demo-full.csv")
categ.df <- read_rds("../shared_data/data/vars-categories.rds")
lang.regions <- read_rds("../shared_data/data/language-curated-regions.rds")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## corr
all <- inner_join(mit.maps,falff.maps)
## entire map
all.c <- corr.table(all%>%select(colnames(mit.maps)[7:12]), all %>% select(starts_with("sub")))
all.c %>%
  filter(V1 %in% colnames(mit.maps), !V2 %in% colnames(mit.maps)) %>%
  ggplot(aes(V1,V2,fill=r))+
  geom_tile()+redblu.col.gradient.2() + my.guides +
  bw.theme + labs(x="",y="",caption = "correlation between fALFF map with mito-maps")
ggsave2("figs/mito-falff-participant.png",4,15)
mit.cor <- all.c %>%
  filter(V1 %in% colnames(mit.maps), !V2 %in% colnames(mit.maps)) %>%
  pivot_wider(names_from = V1, values_from = r, id_cols = V2) %>%
  rename(te_id=V2) %>% mutate(te_id=sub("sub__","",te_id))
################################################################################
################################################################################
################################################################################
all.2 <- demo %>% mutate(sexM=as.numeric(sex=="Male"),
                         right_handed=as.numeric(handedness=="R"),
                         ADHD_dx=as.numeric(ADHD_dx),
                         ASD_dx=as.numeric(ASD_dx)) %>%
  select(te_id, MRI_age, sexM, right_handed, ADHD_dx, ASD_dx) %>%
  inner_join(m1.m2[,-2]) %>%
  inner_join(mit.cor)
corr.table(all.2 %>% select(-colnames(mit.cor)),
           all.2 %>% select(colnames(mit.cor)[-1])) %>%
  filter(V1 %in% colnames(mit.cor), !V2 %in% colnames(mit.cor)) %>%
  mutate(var=V2) %>% left_join(categ.df) %>%
  mutate(var_clean = ifelse(is.na(var_clean), 
                            str_replace_all(sub("_age_correc.*","",var),"_"," "),
                            var_clean),
         cat1 = ifelse(is.na(cat1),"",cat1)) %>%
  ggplot(aes(V1,var_clean,fill=r,label=ifelse(pval<0.05,"*","")))+
  geom_tile()+redblu.col.gradient.2() + my.guides +
  geom_text() + ggh4x::facet_grid2(rows = vars(cat1),scales = "free", space = "free") +
  bw.theme + labs(x="",y="",caption = paste0("participants' fALFF maps were correlated with mito-maps\n",
                                             "then, correlated with other vars"))
ggsave2("figs/mito-falff-participant-to-cogn-and-demo.png",6,10)
################################################################################
################################################################################
################################################################################
## get the correlation between participants' fALFF and MitoD per region
c.s <- unique(all$cerebr_structure[!is.na(all$cerebr_structure)])
c.s.mit.corr <- foreach(ci = 1:length(c.s), .combine = rbind) %dopar% {
  c.s.n <- c.s[ci]
  df <- all %>% filter(cerebr_structure == c.s.n)
  do.call(rbind,lapply(df %>% select(starts_with("sub")), function(x){
    apply(df %>% select(CI,CII,CIV,MitoD,MRC,TRC), 2, function(y){
      cor(y,x)
    })
  })) %>% as.data.frame() %>%
    rownames_to_column("te_id") %>% mutate(te_id = sub("sub__","",te_id),
                                           cerebr_structure = c.s.n)
}
# now correlate these with others?
c.s.all.2 <- demo %>% mutate(sexM=as.numeric(sex=="Male"),
                             right_handed=as.numeric(handedness=="R"),
                             ADHD_dx=as.numeric(ADHD_dx),
                             ASD_dx=as.numeric(ASD_dx)) %>%
  select(te_id, MRI_age, sexM, right_handed, ADHD_dx, ASD_dx) %>%
  inner_join(m1.m2[,-2]) %>%
  inner_join(c.s.mit.corr %>% pivot_longer(cols = c(-te_id, -cerebr_structure)) %>%
               mutate(name2 = paste0(name,"____", cerebr_structure)) %>%
               pivot_wider(names_from = name2, values_from = value, id_cols = te_id))
corr.table(c.s.all.2 %>% select(-starts_with(colnames(c.s.mit.corr)[2:7]),-te_id),
           c.s.all.2 %>% select(starts_with(colnames(c.s.mit.corr)[2:7]))) %>%
  mutate(measure = sub("___.*","",V1), cerebr_structure = sub(".*___","",V1),var=V2) %>%
  filter(cerebr_structure %in% unique(c.s.mit.corr$cerebr_structure),
         !sub("__.*","",V2) %in% colnames(c.s.mit.corr)[2:7],
         !grepl("Ventricle|CSF|lobule",cerebr_structure)) %>%
  left_join(categ.df) %>%
  mutate(var_clean = ifelse(is.na(var_clean), str_replace_all(var,"_"," "),var_clean),
         cat1 = ifelse(is.na(cat1),"",cat1),
         hemisphere = sub("_.*","",cerebr_structure),
         var_clean = factor(var_clean, levels = c(unique(categ.df$var_clean), 
                                                  "ASD dx", "ADHD dx",
                                                  "MRI age", "sexM",
                                                  "right handed"))) %>%
  filter(measure %in% c("MitoD","MRC")) %>%
  # filter(cerebr_structure %in% lang.regions$y_clean) %>%
  ggplot(aes(var_clean,sub(".H_","", cerebr_structure),fill=r,
             label=ifelse(pval<0.05,"*","")))+
  geom_tile()+redblu.col.gradient.2() + my.guides +
  geom_text() + 
  ggh4x::facet_nested(cols = vars(measure,cat1), rows = vars(hemisphere),
                      scales = "free", space = "free",
                      nest_line = element_line(linewidth=0.6), solo_line = T) +
  bw.theme + theme(axis.text.x.bottom = element_text(angle=90,hjust=1,vjust=0.5)) +
  labs(x="",y="",caption = paste0("participants' fALFF maps were correlated with mito-maps on an ROI-level\n",
                                  "then, correlated with other vars"))
ggsave2("figs/mito-falff-participant-to-cogn-and-demo_roi-level.png",18,20)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## psvc example
psvc.design.r <- readxl::read_xlsx("../shared_data/data/RPOE_meta.xlsx",sheet = "PSVC-fMRI-metadata")  %>% 
  drop_na(task) %>% mutate(start_frame= start_frame+1, end_frame=end_frame+1) %>%
  select(order_in_sequence, task,exp_condition, category, specific,start_frame, end_frame)
psvc.design <- psvc.design.r %>%
  rowwise() %>% mutate(frame = list(seq(start_frame, end_frame))) %>%
  unnest(frame) %>% select(task, exp_condition, category, specific, frame)
psvc.act <- read_rds("data/derivatives/func/USE_THIS/JM/PS-VC/B_sub-JM_PS-VC_MNI-reg-voxel-ts-R.rds")
tt <- inner_join(mit.maps, psvc.act)
image(cor(tt[,-c(1:23)],tt$MRC))

## entire brain map
mit.corr <- data.frame(frame = colnames(tt[,-c(1:18)])) %>%
  mutate(frame = parse_number(frame),
         CI = cor(tt[,-c(1:18)],tt$CI),
         CII = cor(tt[,-c(1:18)],tt$CII),
         CIV = cor(tt[,-c(1:18)],tt$CIV),
         MitoD = cor(tt[,-c(1:18)],tt$MitoD),
         MRC = cor(tt[,-c(1:18)],tt$MRC),
         TRC = cor(tt[,-c(1:18)],tt$TRC))
left_join(psvc.design, mit.corr) %>%
  pivot_longer(cols = c(CI,CII,CIV,MitoD, TRC, MRC), values_to = "val") %>%
  filter(name=="MRC") %>%
  ggplot(aes(frame,val,color=name)) +
  geom_line() +
  scale_color_manual(values = palette.1,name="measure") +
  geom_vline(aes(xintercept = start_frame), data = psvc.design.r,alpha=0.1)+
  geom_text(aes(start_frame,0.26,label=paste0(task,"_",category)), color = "black",
            data = psvc.design.r,angle=90,hjust=1,size=2.5) +
  bw.theme + labs(y="r")
ggsave2("figs/processing-tmp/PSVC-mito-corr-JM.png",24,4)

## roi-level
psvc.mit.corr <- foreach(ci=1:length(c.s), .combine = rbind) %dopar% {
  c.s.n <- c.s[ci]
  df <- tt %>% filter(cerebr_structure == c.s.n)
  do.call(rbind,lapply(df %>% select(starts_with("Frame")), function(x){
    apply(df %>% select(CI,CII,CIV,MitoD,MRC,TRC), 2, function(y){
      cor(y,x)
    })
  })) %>% as.data.frame() %>%
    rownames_to_column("frame") %>% 
    mutate(frame = parse_number(frame),cerebr_structure = c.s.n)
}
left_join(psvc.design, psvc.mit.corr) %>%
  pivot_longer(cols = c(CI,CII,CIV,MitoD, TRC, MRC), values_to = "val") %>%
  filter(name %in% c("MRC")) %>%
  filter(cerebr_structure %in% lang.regions$y_clean) %>%
  mutate(hemisphere = sub("_.*", "", cerebr_structure),
         struct2 = sub(".H_","",cerebr_structure)) %>%
  filter(task == "word_association") %>%
  ggplot(aes(frame,val,color=hemisphere)) +
  geom_line() +
  scale_color_manual(values = redblu.col.2[c(2,1)],name="measure") +
  geom_vline(aes(xintercept = start_frame), 
             data = psvc.design.r %>% filter(task == "word_association"),alpha=0.1)+
  # geom_text(aes(start_frame,0.26,label=paste0(task,"_",category)), color = "black",
  #           data = psvc.design.r,angle=90,hjust=1,size=2.5) +
  ggh4x::facet_grid2(rows = vars(struct2), scales="free") +
  bw.theme + labs(y="r")

################################################################################
################################################################################
psvc.act.avg <- psvc.act %>% group_by(cerebr_structure) %>%
  summarise_at(.vars = vars(starts_with("Frame")), .funs = function(x) median(x)) %>%
  pivot_longer(cols = starts_with("Frame"), names_to = "frame", values_to = "int") %>%
  mutate(frame = parse_number(frame))
left_join(psvc.design.r, psvc.act.avg) %>%
  ggplot(aes(frame,int)) +
  geom_line(aes(group=cerebr_structure)) +
  geom_vline(aes(xintercept = start_frame), data = psvc.design,alpha=0.1)+
  geom_text(aes(start_frame,0.26,label=paste0(task,"_",category)), color = "black",
            data = psvc.design,angle=90,hjust=1,size=2.5) +
  bw.theme + labs(y="r")
ggsave2("figs/processing-tmp/PSVC-mito-corr-JM.png",24,4)

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