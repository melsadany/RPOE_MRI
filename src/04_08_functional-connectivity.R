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
## identify participants
good.r1 <- c(paste0("2E_00", c(1:6,8,9)),
             paste0("2E_0", c(10:13,16:20,22,24,25,27:30,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
             paste0("2E_", c(100,102,104:106,108,109,112,115,124,126,118,131,133,134)))
good.psvc <- c(paste0("2E_0", c(22,29,23,31:35,38:45,48,50:57,66,70,75,84,85,90,95:99)),
               paste0("2E_", c(100,102,104,105,106,108,109,112,115,124,126,118,131,133,134)))

roi.meta <- read_rds("../../../refs/Schaefer2018/Parcellations/MNI/roi-meta.rds")
################################################################################
################################################################################
################################################################################
################################################################################
## function for connectivity calculation
calculate_connectivity <- function(df, compcorr=F) {
  
  if (compcorr) {
    df2 <- as.matrix(df[,-1])
    # 1. Identify high-variance voxels (potential noise sources)
    variance_by_voxel <- apply(df2, 1, var)
    variance_threshold <- quantile(variance_by_voxel, probs = 0.95)  # Top 5%
    noise_voxels <- which(variance_by_voxel > variance_threshold)
    # 2. Extract time courses from these high-variance voxels
    noise_timecourses <- df2[noise_voxels, ]  # Noise voxels × time
    # 3. Perform PCA on the noise timecourses
    # Transpose to time × noise_voxels
    noise_pca <- prcomp(t(noise_timecourses), center = TRUE, scale. = FALSE)
    # 4. Select number of components (usually based on variance explained)
    noise_variance <- noise_pca$sdev^2 / sum(noise_pca$sdev^2)
    cumulative_noise <- cumsum(noise_variance)
    # Typically keep enough components to explain 50-70% of noise variance
    n_compcor <- which(cumulative_noise > 0.5)[1]  # First to exceed 50%
    cat(sprintf("Using %d CompCor components (explain %.1f%% variance)\n",n_compcor, cumulative_noise[n_compcor] * 100))
    # 5. Extract the CompCor regressors
    compcor_regressors <- noise_pca$x[, 1:n_compcor]
    # 6. Regress these out from ALL voxels
    fmri_compcor_cleaned <- foreach(voxel=1:nrow(df2),.combine = rbind) %dopar% {
      z_from_lm(y = df2[voxel, ],x = compcor_regressors)
    }
    # average per ROI
    roi_ts <- cbind(df[,1],fmri_compcor_cleaned) %>% group_by(roi_name) %>%
      summarise_all(.funs = mean) %>% column_to_rownames("roi_name")
    # 7. Calculate functional connectivity
    fc_compcor <- cor(t(roi_ts), use = "pairwise.complete.obs",method = "pearson")
    # 8. Apply Fisher's z-transform for normality
    connectivity_z <- atanh(fc_compcor)
    # 9. reframe and keep unique connections
    roi_names <- colnames(connectivity_z)
    # 10. Get upper triangle indices
    upper_tri_indices <- upper.tri(connectivity_z, diag = F)
    # 11. Create connection names
    connection_pairs <- expand.grid(ROI1 = roi_names, ROI2 = roi_names)
    connection_pairs$index <- paste(connection_pairs$ROI1, connection_pairs$ROI2, sep = "--")
    # 12. Filter to upper triangle
    upper_conn <- connection_pairs[upper_tri_indices, ]
    upper_conn$value <- connectivity_z[upper_tri_indices]
    
  } else {
    # 1. Aggregate to ROI level
    roi_ts <- df %>%
      group_by(roi_name) %>%
      summarise(across(starts_with("Frame"), mean, na.rm = T)) %>%
      column_to_rownames("roi_name") %>%
      t()  # Time points as rows
    # 2. Scale the time series
    roi_ts_scaled <- scale(roi_ts, center = TRUE, scale = TRUE)
    # 3. Calculate correlations
    connectivity <- cor(roi_ts_scaled, use = "pairwise.complete.obs",method = "pearson")
    # 4. Apply Fisher's z-transform for normality
    connectivity_z <- atanh(connectivity)
    # 5. reframe and keep unique connections
    roi_names <- colnames(connectivity_z)
    # Get upper triangle indices
    upper_tri_indices <- upper.tri(connectivity_z, diag = F)
    # Create connection names
    connection_pairs <- expand.grid(ROI1 = roi_names, ROI2 = roi_names)
    connection_pairs$index <- paste(connection_pairs$ROI1, connection_pairs$ROI2, sep = "--")
    # Filter to upper triangle
    upper_conn <- connection_pairs[upper_tri_indices, ]
    upper_conn$value <- connectivity_z[upper_tri_indices]
  }
  
  return(upper_conn)
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## resting-state
registerDoMC(10)
rest.connectivity <- foreach(ii = 1:length(good.r1)) %dopar% {
  id.n <- good.r1[ii]
  file.n <- paste0(project.dir,"/data/derivatives/func/USE_THIS/",id.n,"/REST1/B_sub-",id.n,"_REST1_MNI-reg-voxel-ts-R.rds")
  if (!file.exists(file.n)) {
    return(NULL)
  }
  id.df <- read_rds(file.n) %>%
    filter(!is.na(roi_name)) %>%
    select(-c(1:11, paste0("Frame_", c(1:5)))) %>%
    mutate(roi_name = sub("17Networks_","",roi_name))
  id.df.clean <- id.df[!((rowSums(id.df[,-1])<1000)),]
  id.conn <- calculate_connectivity(df=id.df.clean, T)
  colnames(id.conn)[4] <- id.n
  write_rds(id.conn,paste0(project.dir,"/data/derivatives/func/USE_THIS/",id.n,"/REST1/B_sub-",id.n,"_REST1_MNI-reg-voxel-connectivity-R.rds"),compress = "gz")
  return(id.conn)
}
t1 <- rest.connectivity %>% reduce(full_join, by=c("ROI1","ROI2","index"))
################################################################################
################################################################################
################################################################################
################################################################################
## PSVC
wat.frames <- read_table("data/derivatives/PSVC-task-TS/all/word_association__word.txt",col_names = F) %>%
  rowwise() %>% mutate(frames = list(c(X1:(X2-1+X1)))) %>% select(frames) %>% unnest(frames)


registerDoMC(10)
task.connectivity <- foreach(ii = 1:length(good.psvc)) %dopar% {
  id.n <- good.psvc[ii]
  file.n <- paste0(project.dir,"/data/derivatives/func/USE_THIS/",id.n,"/PS-VC/B_sub-",id.n,"_PS-VC_MNI-reg-voxel-ts-R.rds")
  file.o <- paste0(project.dir,"/data/derivatives/func/USE_THIS/",id.n,"/PS-VC/B_sub-",id.n,"_PS-VC_MNI-reg-voxel-WAT-connectivity-R.rds")
  if (file.exists(file.o)) {
    id.conn <- read_rds(file.o)
  } else {
    if (!file.exists(file.n)) {
      return(NULL)
    }
    id.df <- read_rds(file.n) %>%
      filter(!is.na(roi_name)) %>%
      # select(-c(1:11, paste0("Frame_", c(1:5)))) %>%
      select(roi_name, paste0("Frame_", unique(wat.frames$frames))) %>%
      mutate(roi_name = sub("17Networks_","",roi_name))
    id.df.clean <- id.df[!((rowSums(id.df[,-1])<1000)),]
    id.conn <- calculate_connectivity(df=id.df.clean,T)
    colnames(id.conn)[4] <- id.n
    write_rds(id.conn,file.o,compress = "gz")
  }
  
  return(id.conn)
}
t2 <- task.connectivity %>% reduce(full_join, by=c("ROI1","ROI2","index"))
################################################################################
################################################################################
################################################################################
################################################################################
list("resting-state"=t1, "task"=t2) %>%
  write_rds("data/derivatives/R-func/functional-connectivity-rest-and-task-pc-corr.rds",compress = "gz")
################################################################################
