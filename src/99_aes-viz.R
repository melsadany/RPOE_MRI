################################################################################
################################################################################
rm(list = ls())
gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(ggseg)
library(ggseg3d)
library(plotly)
library(ggsegSchaefer)
library(ciftiTools)
library(gifti)
################################################################################
################################################################################
project.dir <- paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                      "/jmichaelson-wdata/msmuhammad/projects/RPOE/mri")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# identify tasks
tasks <- sub("\\.txt", "", list.files("data/derivatives/PSVC-task-TS/afni/"))[-2]

major.t <- tasks[c(6,7,31,56,46)]
minor.t.1 <- tasks[c(2:5,13,18)]
minor.t.2 <- tasks[c(32:45,47:55,57)]
minor.t.3 <- tasks[c(8:12,14:17,19:30)]
tasks.of.int <- c(major.t, minor.t.1)
################################################################################
################################################################################
# statistics map
ciftiTools.setOption("wb_path", "../../../workbench/connectome-wb/workbench")
data.folder <- paste0("data/derivatives/func/AFNI-group-level-maps")


registerDoMC(cores = 4)
foreach(tt = 1:length(tasks.of.int)) %dopar% {
  task.n <- tasks.of.int[tt]
  
  # task.nii <- readNIfTI(paste0("data/derivatives/func/AFNI-group-level-maps/",
  #                              task.n, "/", task.n, "_MEMA_group-Coef.nii"))
  
  l.surf <- read_gifti("../../../refs/mni-freesurfer/2mm/surf/lh.white.surf.gii")
  r.surf <- read_gifti("../../../refs/mni-freesurfer/2mm/surf/rh.white.surf.gii")
  
  rh.gi <- read_gifti(paste0(data.folder, "/", task.n, "/", task.n, "_MEMA_group-Coef_rh.stats.shape.gii"))
  lh.gi <- read_gifti(paste0(data.folder, "/", task.n, "/", task.n, "_MEMA_group-Coef_lh.stats.shape.gii"))
  sc.gi <- read_gifti(paste0(data.folder, "/", task.n, "/", task.n, "_MEMA_group-Coef_sc.stats.shape.gii"))
  
  rh_data <- as.numeric(unlist(rh.gi$data))
  lh_data <- as.numeric(unlist(lh.gi$data))
  sc_data <- as.numeric(unlist(sc.gi$data))
  
  # Reshape subcortical data into a matrix or array
  # Assuming it represents a single frame, reshape to 1-row matrix
  sc_matrix <- matrix(sc_data, nrow = 1)
  # Create placeholder subcortical labels (matching the data dimensions)
  sc_labels <- matrix(seq_len(ncol(sc_matrix)), nrow = 1)
  # Create a mask for the subcortical volume
  sc_mask <- matrix(1, nrow = 1, ncol = ncol(sc_matrix))  # All voxels are valid
  
  task.xii <- as_cifti(cortexL = lh_data,cortexR = rh_data, 
                       surfL = l.surf,
                       surfR = r.surf)
  
  
  
  
  
  xii <- read_xifti(paste0("data/derivatives/func/AFNI-group-level-maps/",
                           task.n, "/", task.n, "_MEMA_group-Coef_cortex.dtseries.nii"))
  xii <- read_cifti(ciftiTools.files()$cifti[1])
  ## Optional manipulations such as smoothing:
  xii <- smooth_xifti(xii,surf_FWHM = 8)
  ## Plot data:
  view_cifti_surface(xii,
                     surfL = 'midthickness',
                     surfR= 'midthickness',
                     colors  = viridis::inferno(n = 500),
                     legend_embed = F,zlim = c(-Inf, Inf),
                     hemisphere = 'both')
  
}



# xii <- read_xifti("data/derivatives/func/AFNI-group-level-maps/baseline/baseline_MEMA_group-Coef_cortex.dscalar.nii") 
xii <- read_cifti(ciftiTools.files()$cifti[4])
## Optional manipulations such as smoothing:
# xii <- smooth_xifti(xii,surf_FWHM = 8)
## Plot data:
view_cifti_surface(xii,
                   # surfL = 'midthickness',
                   # surfR= 'midthickness',
                   # colors  = viridis::inferno(n = 500),
                   # legend_embed = F,zlim = c(-Inf, Inf),
                   hemisphere = 'both')

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
# task fMRI?
library(gganimate)
# select 1 participant
d01 <- psvc.bold %>%
  filter(te_id == "2E_023") %>%
  mutate(value = scale(value, scale = T, center = T)[,1],
         label = paste0(tolower(hemisphere), "_", `ROI Name`))

schaefer7_100_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>%
  mutate(atlas_label=label) %>%
  left_join(d01 %>% 
              filter(frame %in% c(1:3)) %>%
              # mutate(frame = as.factor(frame)) %>%
              select(label, value, frame)) %>%
  drop_na() %>%
  select(-label) %>%
  rename(label=atlas_label) %>%
  ggplot() +
  geom_brain(atlas = schaefer7_100,
             position = position_brain(hemi ~ side),
             aes(fill = value)) +
  scale_fill_gradient2(low = "#5e5a6d",
                       high = "red") +
  my.guides +
  theme(axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_blank()) +
  transition_time(frame, range = c(1,3))
  
  
  
#
  
  
  
  
  ggseg3d(atlas = schaefer7_100_3d,
          surface = list("lh" = schaefer7_100_3d$surf[1],
                         "rh" = schaefer7_100_3d$surf[2]),
          hemisphere = c("right"),frame = "frame",
          colour = "value", text = "value",
          palette = c("forestgreen" = min(d01$value), 
                      "white" = 0, 
                      "firebrick" = max(d01$value))) %>%
  remove_axes() %>% 
  # add_glassbrain() %>%
  layout(scene = list(
    camera = list(eye = list(x = 0, y = 0, z = -2),
                  up = list(x = 0, y = -1, z = 0),
                  center = list(x = 0, y = 0, z = 0))),
    title = paste0("brain")) %>%
  config(toImageButtonOptions = list(
    format = "svg",
    filename = "brain",
    width = 800,
    height = 700)) %>%
  animation_button(x = 1, y = 0, 
                   xanchor = "right", 
                   yanchor = "bottom")

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