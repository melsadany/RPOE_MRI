################################################################################
#                       correlate fMRI log data and IQ tests                   #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri"
setwd(project.dir)
################################################################################
################################################################################
psvc.design.r <- readxl::read_xlsx("../shared_data/data/RPOE_meta.xlsx",
                                   sheet = "PSVC-fMRI-metadata") %>% 
  drop_na(task) %>%
  mutate(start_frame= start_frame+1, end_frame=end_frame+1)
################################################################################
################################################################################
# read fmri processed data
p.list <- list.dirs("/Dedicated/jmichaelson-sdata/MRI/RPOE", recursive = F)
registerDoMC(cores = 20)
mri.meta <- foreach (i = 1:(length(p.list)-1), .combine = rbind) %dopar% {
  f.path <- p.list[i]
  pid <- sub("/Dedicated/jmichaelson-sdata/MRI/RPOE/", "", f.path)
  f.path2 <- paste0(project.dir, "/data/derivatives/MRI-log/mdata/", pid, "_mdata.tsv")
  if (file.exists(f.path2)) {
    t <- read_tsv(f.path2) %>%
      mutate(te_id = pid) %>%
      mutate(keypress = readr::parse_number(keypress)) %>%
      mutate(hand = ifelse(keypress %in% c(6:9), "R", 
                           ifelse(keypress %in% c(1:5), "L", NA))) %>%
      mutate(converted_keypress = ifelse(hand == "L", 
                                         keypress+5, keypress),
             keypress_rel_time = ifelse(keypress_rel_time>0, keypress_rel_time, NA),
             order_in_sequence = psvc.design.r$order_in_sequence) %>%
      relocate(order_in_sequence, .after = file)
    return(t)
  }else {
    return(NULL)
  }
}
################################################################################
################################################################################
task.template <- mri.meta %>% 
  mutate(template_answer = ifelse(exp_condition == "same", 7, 
                                 ifelse(exp_condition == "diff", 8, 
                                        ifelse(exp_condition == "coherent", 7,
                                               ifelse(exp_condition == "incoherent", 8,
                                                      NA))))) %>%
  select(task, exp_condition, category, specific, file, order_in_sequence, 
         duration, template_answer) %>%
  slice_head(n=209)
task.total <- task.template %>%
  mutate(category = ifelse(task == "semantic_coherence", NA, category)) %>%
  group_by(task, exp_condition, category) %>%
  dplyr::summarise(task_total = n())
################################################################################
################################################################################

# add the task answer and a logic indicating if answered correctly or not
# some participants pressed the wrong key during some tasks
# you need to correct those

# 2E_095 pressed key 6 instead of 7 for some tasks
# 2E_040: 
#     baseline: pressed 8 instead of 7
#     file images_2/prompt_16_same.png: pressed 8 instead of 7
#     file images_2/prompt_33_diff.png: pressed 9 instead of 8
#     file images_2/prompt_15_diff.png: pressed 9 instead of 8
#     file images_2/prompt_13_diff.png: pressed 9 instead of 8
#     file images_2/prompt_48_same.png: pressed 8 instead of 7
#     file images_2/prompt_18_same.png: pressed 8 instead of 7
#     I would recommend dropping this one!


psvc.log <- mri.meta %>%
  left_join(task.template %>% distinct(file, template_answer) %>% drop_na()) %>%
  filter(!te_id %in% c("2E_040")) %>%
  mutate(converted_keypress = ifelse(te_id == "2E_095" & converted_keypress == 6, 7, converted_keypress),
         answered_correctly = case_when((!is.na(template_answer))&is.na(converted_keypress) ~ F,
                                        (!is.na(template_answer))&(converted_keypress != template_answer) ~ F,
                                        (!is.na(template_answer))&(converted_keypress == template_answer) ~ T),
         drop = case_when(is.na(answered_correctly) ~ F,
                          answered_correctly == F ~ T,
                          answered_correctly == T ~ F))
write_rds(psvc.log, "data/derivatives/MRI-log/combined-log-files.rds")

################################################################################
# how many were correct per participant
p.scores <- mri.meta %>%
  left_join(task.template) %>%
  distinct() %>%
  filter(!is.na(template_answer)) %>%
  select(te_id, task, exp_condition, category,file, converted_keypress, template_answer) %>%
  mutate(answered_correctly = ifelse(is.na(converted_keypress), F,
                                     ifelse(converted_keypress == template_answer, T, F)),
         category = ifelse(task == "semantic_coherence", NA, category),
         exp_condition = paste0(exp_condition, "_", category)) %>%
  group_by(te_id, task, exp_condition) %>%
  dplyr::summarise(total = sum(answered_correctly)) %>%
  left_join(task.total, relationship = "many-to-many") %>% 
  mutate(task_exp = paste0(task, "_",exp_condition, "_score")) %>% 
  pivot_wider(names_from = task_exp, values_from = total, id_cols = te_id) %>% ungroup() %>%
  rowwise() %>%
  mutate(total_correct = sum(c_across(-c(contains("baseline"), te_id))))

# get the relative time needed to define coh/non-coh and same/diff for the ones answered
p.scores.t <- mri.meta %>%
  filter(task %in% c("PS_samediff", "semantic_coherence")) %>%
  mutate(abs_keypress_t = keypress_rel_time - rel_time,
         category = ifelse(task == "semantic_coherence", NA, category),
         exp_condition = paste0(exp_condition, "_", category)) %>%
  select(te_id, task, exp_condition, abs_keypress_t) %>%
  drop_na() %>%
  group_by(te_id, task, exp_condition) %>%
  dplyr::summarise(avg_time_for_same_coh = mean(abs_keypress_t)) %>%
  mutate(task_2 = paste0(task, "_", exp_condition, "_avg_t")) %>%
  pivot_wider(names_from = task_2, values_from = avg_time_for_same_coh, id_cols = "te_id") %>% ungroup() %>%
  rowwise() %>%
  mutate(avg_time_sd_ci = mean(c_across(-c(te_id)), na.rm = T))
################################################################################
# how many words were thought of
word.count <- mri.meta %>%
  filter(task == "word_association", exp_condition == "word") %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(te_id) %>%
  dplyr::summarise(total_word_association_words = sum(count))
################################################################################
# get the relative time or how long it takes them to read the number in mind
# number.presses <- mri.meta %>%
#   filter(task == "word_association", exp_condition == "number") %>%
#   mutate(press = ifelse(is.na(keypress), F, T)) %>%
#   select(te_id, press) %>%
#   group_by(te_id) %>%
#   dplyr::summarise(total_number_press = sum(press))
number.presses.time <- mri.meta %>%
  filter(task == "word_association", exp_condition == "number") %>%
  mutate(abs_keypress_t = keypress_rel_time - rel_time) %>%
  select(te_id, abs_keypress_t) %>%
  drop_na() %>%
  group_by(te_id) %>%
  dplyr::summarise(avg_time_number_reading = mean(abs_keypress_t, na.omit = T))
# same thing for RAN
# participant 2E_032 had no RAN response
RAN.presses.time <- mri.meta %>%
  filter(task == "RAN") %>%
  mutate(abs_keypress_t = keypress_rel_time - rel_time) %>%
  select(te_id, abs_keypress_t) %>%
  drop_na() %>%
  group_by(te_id) %>%
  dplyr::summarise(avg_time_RAN_reading = mean(abs_keypress_t, na.omit = T))
################################################################################
# read tests data
m1.m2 <- read_rds("../shared_data/data/m1m2.rds")
m1.m2.s <- read_rds("../shared_data/data/m1m2-sex-corrected.rds")
################################################################################
################################################################################
m0 <- inner_join(p.scores, p.scores.t) %>%
  inner_join(word.count) %>%
  inner_join(number.presses.time) %>%
  full_join(RAN.presses.time) # one participant w no response
write_rds(m0, "data/derivatives/R-func/participants-PSVC-log-results.rds")
################################################################################
################################################################################
# correlation between fmri scores and iq/nih
m123 <- inner_join(inner_join(p.scores, p.scores.t), 
                   m1.m2.s) %>%
  inner_join(word.count) %>%
  inner_join(number.presses.time) %>%
  inner_join(RAN.presses.time) 

m123%>%
  rename(pattern_comparison = pattern_comparison_PS_age_corrected_standard_score) %>% 
  pivot_longer(cols = c(pattern_comparison,
                        # PSI_composite_score,
                        Coding, Symbol_Search), 
               names_to = "IQ_task", values_to = "IQ_score") %>%
  pivot_longer(cols = c(total_correct,
                        total_word_association_words,
                        avg_time_sd_ci,
                        avg_time_number_reading,
                        avg_time_RAN_reading),
               names_to = "PSVC_task", values_to = "PSVC_score") %>%
  group_by(PSVC_task) %>% 
  mutate(PSVC_score = scale(PSVC_score, scale = T, center = T)[,1]) %>% 
  ungroup() %>%
  mutate(PSVC_task = case_when(PSVC_task == "total_correct" ~ "count of correct answers",
                               PSVC_task == "total_word_association_words" ~ "count of words generated",
                               PSVC_task == "avg_time_sd_ci" ~ "avg time for answers",
                               PSVC_task == "avg_time_number_reading" ~ "avg time for reading single numbers",
                               PSVC_task == "avg_time_RAN_reading" ~ "avg time for fininshing RAN")) %>%
  ggplot(aes(y = IQ_score, x = PSVC_score)) +
  geom_point(shape = 1, alpha = 0.6) +
  geom_smooth(method = "lm", color = six.colors[3]) + 
  ggpubr::stat_cor(color = "red", method = "spearman") +
  ggh4x::facet_grid2(cols = vars(PSVC_task), rows = vars(IQ_task),
                     scales = "free", independent = T) +
  labs(y = "Z-scaled IQ score (residualized for sex)",
       x = "Z-scaled feature value",
       caption = paste0("n(samples): ", length(unique(m123$te_id)))) +
  bw.theme
  # theme_linedraw() +
  # theme(strip.background = element_rect(fill = "white", color = "white"),
  #       strip.text = element_text(color = "black"))
ggsave(filename = "figs/corr-psvc-fmri-performance-w-PS.png",
       width = 13, height = 8, units = "in", dpi = 320, bg = "white")




corr.table(m123 %>% select(colnames(m1.m2.s)[c(3:24)]),
           m123 %>% select(-colnames(m1.m2.s))) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  filter(V1 %in% colnames(m1.m2.s), !V2 %in% colnames(m1.m2.s),
         ! grepl("baseline", V2)) %>%
  mutate(V1 = sub("_age_corrected_standard_score", "_NIH", V1),
         V1 = factor(V1, levels = unique(V1)),
         V2 = factor(V2, levels = unique(V2)),
         cat2 = ifelse(grepl("NIH", V1), "NIH-TB", "IQ"),
         V1 = sub("_NIH", "", V1),
         V1 = factor(V1, levels = unique(V1)),
         cat1 = ifelse(grepl("total", V2)|grepl("_score", V2), "total",
                      ifelse(grepl("time", V2)|grepl("_t", V2), "time", "other")),
         V2 = sub("_NA", "", V2)) %>%
  ggplot(aes(x=V1, y=V2, fill=r, label=ifelse(FDR<0.05, "**", ifelse(pval<0.05, ".", "")))) +
  geom_tile() +
  geom_text(size=3, color = "white") +
  ggh4x::facet_grid2(cols = vars(cat2), rows = vars(cat1), scales = "free", space = "free") +
  redblu.col.gradient.2() + my.guides + null_labs + bw.theme +
  labs(caption = paste0("n(samples): ", nrow(m123))) +
  theme(axis.text.x.bottom = element_text(angle=90, hjust=1, vjust=0.5))
ggsave2("figs/corr-psvc-fmri-performance-w-IQ.png",width = 10, height = 8)
################################################################################


################################################################################


################################################################################