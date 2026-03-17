################################################################################
#                               PSVC fMRI task design                          #
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(fmri)
library(ggpubr)
################################################################################
# script is for extracting the BOLD signal for the fMRI task
# mainly trying to get activation per ROI for each stimulus
################################################################################
project.dir <- paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                      "/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
# read task metadata/design
psvc.design.r <- readxl::read_xlsx("../shared_data/data/RPOE_meta.xlsx",
                                   sheet = "PSVC-fMRI-metadata") %>% 
  drop_na(task) %>%
  mutate(start_frame= start_frame+1, end_frame=end_frame+1) %>%
  # drop the first 5 seconds/frames
  filter(order_in_sequence!=0) %>%
  mutate(start_frame=start_frame-5,end_frame=end_frame-5)
################################################################################
################################################################################
# create the design matrix for the task using the HRF model

# fixed params
TR <- 1
scans <- 1215
ttt <- psvc.design.r %>% filter(task %in% c("instructions","baseline","PS_samediff", "RAN", "semantic_coherence", "word_association"))
ntrials <- dim(ttt)[1]

# identify the different stimuli
ind.instructions <- (1:ntrials)[ttt$task=="instructions"]
ind.baseline <- (1:ntrials)[ttt$task=="baseline"]
ind.same <- (1:ntrials)[ttt$task=="PS_samediff"&ttt$exp_condition=="same"]
ind.diff <- (1:ntrials)[ttt$exp_condition=="diff"]
ind.samediff <- (1:ntrials)[ttt$task=="PS_samediff"]
ind.ran <- (1:ntrials)[ttt$exp_condition=="RAN"]
ind.coh <- (1:ntrials)[ttt$exp_condition=="coherent"]
ind.incoh <- (1:ntrials)[ttt$exp_condition=="incoherent"]
ind.cohincoh <- (1:ntrials)[ttt$exp_condition%in%c("coherent","incoherent")]
ind.word <- (1:ntrials)[ttt$exp_condition=="word"]
ind.num <- (1:ntrials)[ttt$exp_condition=="number"]

onsets <- ttt$start_frame
duration <- ttt$duration

# make the HRF for each stimulus
hrf.instructions <- fmri.stimulus(scans, onsets[ind.instructions], duration[ind.instructions], TR = TR, times = T)
hrf.baseline <- fmri.stimulus(scans, onsets[ind.baseline], duration[ind.baseline], TR = TR, times = T)
hrf.same <- fmri.stimulus(scans, onsets[ind.same], duration[ind.same], TR = TR, times = T)
hrf.diff <- fmri.stimulus(scans, onsets[ind.diff], duration[ind.diff], TR = TR, times = T)
hrf.samediff <- fmri.stimulus(scans, onsets[ind.samediff], duration[ind.samediff], TR = TR, times = T)
hrf.ran <- fmri.stimulus(scans, onsets[ind.ran], duration[ind.ran], TR = TR, times = T)
hrf.coh <- fmri.stimulus(scans, onsets[ind.coh], duration[ind.coh], TR = TR, times = T)
hrf.incoh <- fmri.stimulus(scans, onsets[ind.incoh], duration[ind.incoh], TR = TR, times = T)
hrf.cohincoh <- fmri.stimulus(scans, onsets[ind.cohincoh], duration[ind.cohincoh], TR = TR, times = T)
hrf.word <- fmri.stimulus(scans, onsets[ind.word], duration[ind.word], TR = TR, times = T)
hrf.num <- fmri.stimulus(scans, onsets[ind.num], duration[ind.num], TR = TR, times = T)

hrf.all <- cbind(hrf.instructions,hrf.baseline, hrf.same, hrf.diff, hrf.ran, hrf.coh, hrf.incoh, hrf.word, hrf.num)
hrf.major <- cbind(hrf.instructions, hrf.baseline, hrf.samediff, hrf.ran, hrf.cohincoh, hrf.word, hrf.num)

# make the design
# if you have any covariates, you can add them here
xdesign.r <- fmri.design(hrf.all, order = 2)
xdesign <- as.data.frame(xdesign.r)
colnames(xdesign)[1:9] <- paste0("hrf_", c("instructions","baseline",
                                           "same", "diff", "RAN", "coherent", 
                                           "incoherent", "word", "number"))
write_csv(xdesign, "data/derivatives/R-func/PSVC-design-minor_del5.csv")

xdesign.m <- fmri.design(hrf.major, order = 2)
xdesign2 <- as.data.frame(xdesign.m)
colnames(xdesign2)[1:7] <- paste0("hrf_", c("instructions","baseline", "samediff", "RAN", "cohincoh", "word", "num"))
write_csv(xdesign2, "data/derivatives/R-func/PSVC-design-major_del5.csv")

# make a plot for the conditions you have
xdesign %>% 
  rownames_to_column("frame") %>%
  pivot_longer(cols = c(2:10), names_to = "condition", values_to = "activation") %>%
  mutate(frame = as.numeric(frame)) %>%
  ggplot(aes(x=frame, y = activation, color = condition)) +
  geom_line(show.legend = F) +
  scale_color_manual(values = palette.1) +
  facet_wrap(~condition) + bw.theme
ggsave2("figs/processing-tmp/PSVC-HRF-design-matrix-2.png", 12, 8)
## check correlation between these HRFs
corr.table(xdesign %>% select(1:8),
           xdesign %>% select(-c(1:8)),
           method = "spearman") %>%
  filter(V1 != V2, grepl("hrf", V1), grepl("hrf", V2)) %>%
  mutate(FDR = p.adjust(pval, method = "fdr"),
         V1=sub("hrf_","",V1),V2=sub("hrf_","",V2),
         V1 = factor(V1, levels = sub("hrf_","",colnames(xdesign))),
         V2 = factor(V2, levels = sub("hrf_","",colnames(xdesign)))) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, paste0("**, ", round(r, 2)), 
                                                  ifelse(pval < 0.05, paste0("., ", round(r,2)), ""))))+
  geom_tile() + geom_text(color = "white") +
  redblu.col.gradient.2() + my.guides +
  bw.theme+ labs(x="", y="", title = "correlation between HRF of the tasks")
ggsave2("figs/corr-of-PSVC-HRFs-minor.png",8,8)
corr.table(xdesign2 %>% select(1:3),
           xdesign2 %>% select(-c(1:3)),
           method = "spearman") %>%
  filter(V1 != V2, grepl("hrf", V1), grepl("hrf", V2)) %>%
  mutate(FDR = p.adjust(pval, method = "fdr"),
         V1=sub("hrf_","",V1),V2=sub("hrf_","",V2),
         V1 = factor(V1, levels = sub("hrf_","",colnames(xdesign2))),
         V2 = factor(V2, levels = sub("hrf_","",colnames(xdesign2)))) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, paste0("**, ", round(r, 2)), 
                                                  ifelse(pval < 0.05, paste0("., ", round(r,2)), ""))))+
  geom_tile() + geom_text(color = "white") +
  redblu.col.gradient.2() + my.guides +
  bw.theme + labs(x="", y="", title = "correlation between HRF of the tasks")
ggsave2("figs/corr-of-PSVC-HRFs-major.png",8,8)



### cleaner task design plot
xdesign2 %>% 
  rownames_to_column("frame") %>%
  pivot_longer(cols = c(2:8), names_to = "condition", values_to = "activation") %>%
  mutate(frame = as.numeric(frame),
         condition = case_when(condition == "hrf_samediff" ~ "same/different task",
                               condition == "hrf_cohincoh" ~ "coherenet/incoherent task",
                               condition == "hrf_baseline" ~ "baseline",
                               condition == "hrf_instructions" ~ "instructions",
                               condition == "hrf_RAN" ~ "RAN task",
                               condition == "hrf_word" ~ "word association task",
                               condition == "hrf_num" ~ "number reading task")) %>%
  ggplot(aes(x=frame, y = activation, color = condition)) +
  geom_line() +
  scale_color_manual(values = palette.1, name = "task") +
  labs(y = "HRF", x = "frame (time in secs)") +
  bw.theme
ggsave2("figs/processing-tmp/PSVC-HRF-design-matrix-major.png", 8,6)


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
# detailed HRFs
ttt <- psvc.design.r %>% 
  filter(task %in% c("instructions","baseline","PS_samediff", "RAN", "semantic_coherence", "word_association"))
ntrials <- dim(ttt)[1]

tasks.v3 <- ttt %>%
  mutate(task_v3 = paste(task, exp_condition, category, sep = "_"))
tasks.v3.1 <- tasks.v3 %>% distinct(task_v3)
scans <- 1215
onsets <- tasks.v3$start_frame
duration <- tasks.v3$duration

registerDoMC(5)
tasks.v3.design <- foreach(i = 1:nrow(tasks.v3.1), .combine = cbind) %dopar% {
  t <- tasks.v3.1$task_v3[i]
  ind.t <- (1:ntrials)[tasks.v3$task_v3==t]
  hrf.t <- fmri.stimulus(scans, 
                         onsets[ind.t], 
                         duration[ind.t], 
                         TR = TR, 
                         times = T)
  return(hrf.t)
}

# make the design
# if you have any covariates, you can add them here
xdesign.3 <- fmri.design(tasks.v3.design, order = 2)
xdesign.32 <- as.data.frame(xdesign.3)
colnames(xdesign.32)[1:nrow(tasks.v3.1)] <- str_replace_all(paste0("hrf_", tasks.v3.1$task_v3),"-","_")
write_csv(xdesign.32, "data/derivatives/R-func/PSVC-design-minor-detailed_del5.csv")

xdesign.32 %>%
  rownames_to_column("frame") %>%
  pivot_longer(cols = c(2:nrow(tasks.v3.1)), names_to = "condition", values_to = "activation") %>%
  mutate(frame = as.numeric(frame)) %>%
  ggplot(aes(x=frame, y = activation, color = condition)) +
  geom_line(show.legend = F) +
  # scale_color_manual(values = palette.1, name = "task") +
  labs(y = "HRF", x = "frame (time in secs)") +
  bw.theme

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# AFNI task TS
ff <- list.files("data/derivatives/PSVC-task-TS/del5", pattern = "\\.txt")
system("mkdir -p data/derivatives/PSVC-task-TS/del5/afni/")
foreach(ii = c(1,3:length(ff))) %dopar% {
  f <- ff[ii]
  read_delim(paste0("data/derivatives/PSVC-task-TS/del5/", f), col_names = F) %>%
    select(1) %>% unlist() %>% as.numeric() %>%
    write_lines(paste0("data/derivatives/PSVC-task-TS/del5/afni/", f), sep = " ")
}
read_delim(paste0("data/derivatives/PSVC-task-TS/del5/instructions.txt"), col_names = F) %>%
  mutate(cc = paste0(X1,"*",X2)) %>% select(cc) %>% unlist() %>%
  write_lines(paste0("data/derivatives/PSVC-task-TS/del5/afni/instructions.txt"), sep = " ")
################################################################################
################################################################################
################################################################################
