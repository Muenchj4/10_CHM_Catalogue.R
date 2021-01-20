#------------------------------------------------------------------------------
# Type: control script 
# Name: 10_CHM_Catalog.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  script creates a canopy height model from generic Lidar 
#              las data using the lidR package. In addition the script cuts the area to 
#              a defined extent using the lidR catalog concept.
#              Furthermore the data is tiled into a handy format even for poor memory 
#              and a canopy height model is calculated
# Data: regular las LiDAR data sets 
# Output: Canopy heightmodel RDS file of the resulting catalog
# Copyright: Chris Reudenbach, Thomas Nauss 2017,2020, GPL (>= 3)
#------------------------------------------------------------------------------

## clean your environment
rm(list = ls()) 
gc()

# 0 - load packages
#-----------------------------

library("future")

options(future.rng.onMisue = "ignore")
# 1 - source files
#-----------------
source(file.path(envimaR::alternativeEnvi(root_folder = "D:/Benutzer/Muench/edu/mpg-envinsys-plygrnd",
                                          alt_env_id = "COMPUTERNAME",
                                          alt_env_value = "PCRZP",
                                          alt_env_root_folder = "F:/BEN/edu/mpg-envinsys-plygrnd"),
                 "/src/000_setup.R"))
unlink(envrmt$path_tmp)

# 2 - define variables
#---------------------

# switch if lasclip is called
lasclip = TRUE
plott=TRUE
set.seed(1000)
#Saftflusshalbmond
coord = c(xmin,ymin,xmax,ymax)

# setup future for parallel
future::plan(multisession, workers = 4L)
set_lidr_threads(4L)
# 3 - start code 
#-----------------

#---- if not clipped yet
if (lasclip){
  ctg <- lidR::readLAScatalog(envrmt$path_lidar_level0)
  projection(ctg) <-25832
  lidR::opt_chunk_size(ctg) = 400
  
  lidR::opt_chunk_buffer(ctg) <- 20
  lidR::opt_output_files(ctg) <- paste0(envrmt$path_tmp,"{ID}_cut") # add output filname template
  ctg@output_options$drivers$Raster$param$overwrite <- TRUE
  
  crop_aoimof_ctg = lidR::clip_rectangle(ctg, 
                                         xleft = coord[[1]], 
                                         ybottom = coord[[2]], 
                                         xright = coord[[3]], 
                                         ytop = coord[[4]])
  projection(crop_aoimof_ctg) <-25832
  lidR::opt_chunk_size(crop_aoimof_ctg) = 650
  lidR::opt_chunk_buffer(crop_aoimof_ctg) <- 50
  saveRDS(crop_aoimof_ctg,file= file.path(envrmt$path_level1,"crop_aoimof.rds"))
  
}