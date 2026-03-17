################################################################################
#                               processing fMRI log data                       #
################################################################################
# this is to process the Psydat files for the PSVC fMRI task
################################################################################
rm(list = ls())
gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri"
setwd(project.dir)
################################################################################
################################################################################
good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))
################################################################################
################################################################################
# read fMRI tabulated data? 
# using Jake's function to process the file here
format_2e_fmri_metadata = function(fn,nframes=1220){ 
  dat = read.table(fn,sep="\t",header=F,stringsAsFactors=F) 
  d1 = dat[grepl("EXP",dat[,2]),] 
  wa.df <- d1[grepl("images/WA",d1$V3)&grepl("New trial",d1$V3),] %>% rename(task_start = V1)
  d1 = d1[grepl("New trial",d1[,3]),] 
  mdata = t(sapply(strsplit(d1[[3]],split="\\, |\\: "),'[',c(4,6,8,10,12,14))) 
  mdata = data.frame(abs_time=d1[[1]], 
                     rel_time=d1[[1]]-d1[[1]][1], 
                     end_time=d1[[1]]-d1[[1]][1]+as.numeric(mdata[,6]), 
                     task=mdata[,2], 
                     exp_condition=mdata[,3], 
                     category=mdata[,4], 
                     specific=mdata[,5], 
                     duration=as.numeric(mdata[,6]), 
                     file=mdata[,1],stringsAsFactors=F) 
  ipt = dat[grepl("DATA",dat[,2]),] 
  ipt = ipt[!grepl("Keypress: 5",ipt[,3],fixed=T),] #what's keypress 5?
  o = outer(ipt[,1],mdata$abs_time,'-') 
  o[o<0] = Inf  
  kp = apply(o,2,which.min) 
  mdata$keypress = ipt[kp,3] 
  mdata$keypress_rel_time = ipt[kp,1] - mdata$abs_time[1] 
  mdata$keypress[mdata$task%in%c("neutral","instructions") | mdata$keypress_rel_time>mdata$end_time] = NA 
  mdata$keypress_rel_time[mdata$task%in%c("neutral","instructions") | mdata$keypress_rel_time>mdata$end_time] = NA 
  rt = mdata$keypress_rel_time - mdata$rel_time # is that how long it takes to press per task?
  # boxplot(tapply(rt,mdata$category,function(x) x),las=2) 
  # table(mdata$category) 
  #nframes = dim(m)[4] #I assumed it's dim[4] of fMRI
  tr = 1.0 
  tt = seq(0,tr*nframes,length.out=nframes) 
  ii = sapply(tt,function(x) rev(which(x-mdata$rel_time>=0))[1]) 
  fdat = data.frame(time=tt,mdata[ii,],stringsAsFactors=F) 
  expr = rep("instructions",nrow(fdat)) 
  expr[fdat$prompt%in%"baseline"] = "baseline" 
  expr[fdat$prompt%in%c("coherent","incoherent")] = "verbal" 
  expr[fdat$prompt%in%c("same","diff")] = "PS" 
  fdat$experiment=as.factor(expr) 
  prompt = sapply(unique(fdat$category),function(x) fdat$category%in%x) 
  
  # Muhammad's section
  # tried to extract how many presses were done per word in word association task
  first <- mdata[171:209,]
  first <- first[order(first$abs_time), ]
  # Function to find the max abs_time from first dataframe before the given abs_time in second dataframe
  find_max_time <- function(abs_time) {
    max_time <- max(first$abs_time[first$abs_time < abs_time])
    ifelse(is.na(max_time), NA, max_time)
  }
  counts <- ipt %>%
    rename(abs_time = V1, keypress=V3) %>%
    filter(abs_time > min(first$abs_time))
  # Add task_b4_start column to the counts dataframe using sapply/mapply
  counts$task_b4_start <- sapply(counts$abs_time, find_max_time)
  tmp <- inner_join(counts %>% select(abs_time=task_b4_start), mdata %>% select(-keypress))  %>%
    group_by(task,exp_condition, category, specific,abs_time) %>%
    dplyr::summarise(count=n())
  # View(left_join(mdata,tmp))
  mdata <- left_join(mdata,tmp)
  
  
  wa.df2<- wa.df %>% mutate(task_end = task_start + 5, task_file = sub(".*file: ","",sub("png.*","png",V3)))%>%select(task_start,task_end,task_file)%>%
    as.data.table()
  setkey(wa.df2,task_start,task_end)
  wa.resp.0 <- ipt %>% filter(V1>min(wa.df$task_start),V1<max(wa.df$task_start))%>%select(keypress_time=V1,key=V3)%>%as.data.table()
  
  wa.keypresses <- foverlaps(x=wa.resp.0[, .(keypress_time, key_start = keypress_time, key_end = keypress_time)],
                             y=wa.df2[, .(task_start, task_end, task_file)],
                             by.x = c("key_start","key_end"),by.y=c("task_start","task_end"),type="within",nomatch = 0L)
  
  return(list(mdata = mdata, fdat = fdat, prompt = prompt,wa_keypresses=wa.keypresses)) 
}
# # EXAMPLE USE  
# fn = "/sdata/MRI/RPOE/2E_045/scan/metadata/2E_045_semanticMap_2023-07-26_11h05.52.863.log" 
# l = format_2e_fmri_metadata(fn, nframes = 1000) 
# questions for Jake
# what's keypress 7,8,4, and 5
# where's the multiple keypresses for the word association task
################################################################################
################################################################################
################################################################################
# apply the function on users that have done the fMRI-PSVC
system(paste0("mkdir -p ", project.dir, "/data/derivatives/MRI-log/mdata"))
p.of.int <- readxl::read_xlsx("../shared_data/data/RPOE_participants_metadata.xlsx", sheet = 1) %>%
  select(devGenes_id, te_id) %>% filter(te_id %in% good.psvc) %>%
  drop_na() %>%
  # 2E_046 incomplete/no resting/sleeping
  # 2E_031 no fMRI
  # 2E_032 has a different version of the task. 
  # 2E_053 has no data recorded
  # 2E_036 no fMRI
  # 2E_037 no t-fMRI/ no other in-person data
  # 2E_047 no fMRI / braces
  # 2E_049 no fMRI
  # 2E_059 needs different processing for possible key switch
  # please drop them
  mutate(mri_valid = ifelse(te_id %in% paste0("2E_0", c(32,46,53,36,37,47,49)), F, T)) %>%
  filter(mri_valid==T) %>%
  mutate(f = paste0("/Dedicated/jmichaelson-sdata/MRI/RPOE/", te_id)) %>%
  filter(te_id %in% paste0("2E_1",c(31,33,34))) # filter this to the new participants
p.list <- p.of.int$f
# participant 2E_095 switched keys in the scanner, so handle them differently


####
psvc.meta <- readxl::read_xlsx("../shared_data/data/RPOE_meta.xlsx", sheet = "PSVC-fMRI-metadata")%>% drop_na()
registerDoMC(cores = 4)
foreach (i = 1:length(p.list)) %dopar% {
  f.path <- p.list[i]
  pid <- sub("/Dedicated/jmichaelson-sdata/MRI/RPOE/", "", f.path)
  file <- list.files(paste0(f.path, "/scan/metadata"), pattern = "log", full.names = T)
  
  if (length(file)!=0) { 
    out <- format_2e_fmri_metadata(fn = file)
    write_rds(out, paste0(project.dir, "/data/derivatives/MRI-log/", pid, "_fmri-all.rds"))
    write_tsv(out[["mdata"]], paste0(project.dir, "/data/derivatives/MRI-log/mdata/", pid, "_mdata.tsv"))
    
    # write timeseries tables for correct responses
    system(paste0("mkdir -p ", "data/derivatives/MRI-log/timings/", pid))
    out$mdata %>% cbind(psvc.meta %>% select(ends_with("frame"))) %>%
      filter(parse_number(keypress)==7,task == "PS_samediff", exp_condition == "same") %>%
      select(start_frame, duration) %>% mutate(tt=1) %>%
      write_tsv(paste0("data/derivatives/MRI-log/timings/", pid,"/same-corr.txt"), col_names = F)
    out$mdata %>% cbind(psvc.meta %>% select(ends_with("frame"))) %>%
      filter(parse_number(keypress)==8, task == "PS_samediff", exp_condition == "diff") %>%
      select(start_frame, duration) %>% mutate(tt=1)%>%
      write_tsv(paste0("data/derivatives/MRI-log/timings/", pid,"/diff-corr.txt"), col_names = F)
    out$mdata %>% cbind(psvc.meta %>% select(ends_with("frame"))) %>%
      filter(parse_number(keypress)==7, task == "semantic_coherence", exp_condition == "coherent") %>%
      select(start_frame, duration) %>% mutate(tt=1)%>%
      write_tsv(paste0("data/derivatives/MRI-log/timings/", pid,"/coh-corr.txt"), col_names = F)
    out$mdata %>% cbind(psvc.meta %>% select(ends_with("frame"))) %>%
      filter(parse_number(keypress)==8, task == "semantic_coherence", exp_condition == "incoherent") %>%
      select(start_frame, duration) %>% mutate(tt=1)%>%
      write_tsv(paste0("data/derivatives/MRI-log/timings/", pid,"/incoh-corr.txt"), col_names = F)
  }
}

## get word association task response metrics
wa.all <- foreach(i=1:nrow(p.of.int),.combine = rbind)%dopar%{
  id.n <- p.of.int$te_id[i]
  out.f <- paste0(project.dir, "/data/derivatives/MRI-log/", id.n, "_fmri-all.rds")
  if (file.exists(out.f)) {
    out <- read_rds(out.f)
    out$wa_keypresses%>%ungroup() %>% group_by(task_file)%>%
      summarise(responses_count=n(),
              thinking_time=list(diff(key_start)),
              onset=key_start[1]-task_start,
              thinking_time_mean = ifelse(responses_count>1,mean(unlist(thinking_time)),NA))%>%
      ungroup()%>%distinct()%>%
      mutate(te_id = id.n)
  } else{
    return(NULL)
  }
}
wa.all %>% group_by(te_id) %>% summarise(total_word_count = sum(responses_count,na.rm=T),
                                         thinking_time_mean = mean(thinking_time_mean,na.rm=T),
                                         onset = mean(onset,na.rm=T))%>%
  write_rds("data/derivatives/MRI-log/word-association-metrics.rds",compress="gz")

##

# participant 2E_053 had something weird/different in their log file
# So, I'll try to extract their keypresses from their csv file
log <- "/sdata/MRI/RPOE/2E_053/scan/metadata/2E_053_semanticMap_2023-09-01_15h46.12.339.log" 
csv <- read_csv("/sdata/MRI/RPOE/2E_053/scan/metadata/2E_053_semanticMap_2023-09-01_15h46.12.339.csv")
tmp <- csv %>% 
  # select(task, key_resp.keys, key_resp.rt) %>%
  mutate(key_resp.keys = strsplit(gsub("[\\[\\]']", "", as.character(key_resp.keys)), ", "),
         key_resp.rt = strsplit(gsub("[\\[\\]']", "", key_resp.rt), ", ")) %>%
  unnest(key_resp.keys, key_resp.rt) %>%
  mutate(keypress = paste0("Keypress: ", parse_number(key_resp.keys)),
         rt = parse_number(key_resp.rt),
         abs_time = key_resp.started,
         keypress_rel_time = key_resp.started+rt) %>%
  filter(!grepl("Keypress: 5",keypress,fixed=T))
# apparently the file has keypresses for numbers only, other than that nothing was recorded
################################################################################

################################################################################


################################################################################




