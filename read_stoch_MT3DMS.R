###############################################################################
# R script to calculate the Central Mass Moments CMM from a stochastic run of MT3DMS
# NOT in a markdown because there are currently some issues with markdown/knitr
# when in a code chunk the directory is changed for reading all *.UCN files
# G.Bier January, 2024
###############################################################################
rm(list = ls())
source("MMfunctions.R")  #file containing several functions to calculate CMM's

# Before running this code assign the proper folder at which the subfolders of
# the stochastic run is located (i.e. folder numbered 001,002,...010)
# assign this folder "between quotes" to variable sto_dir.
# Next the UCN files will be read from the subsequent folders (001,002,..010)
#
# Assign the number of time steps which need to be considered.
# It can be found in the sub panel in the explorer window (lower left on screen)
# when a solute/species is selected 
# This number shodiruld be assigned to variable NTTS

sto_dir = "stoch_144_72_k7_5_std7_5_sim250_MT3DMS"
NTTS = 22

cur_dir = getwd()  #get currrent directory to go back to this after running

setwd(paste0(cur_dir,"/",sto_dir,"/")) #setting the workin directory for reading the subsequent folders
print(paste("Folder to read UCN files : ",getwd()))
sto_subdirs = getwd()

stoch_MM_df = c() # open a data container for the results

###########calculating the Volume per cell
delx = 0.25#1.0 #10. 
dely = 0.25#1.0 #10.
delz = 0.25#1.0 #10.
porosity =0.4
volume = delx*dely*delz*porosity
x_origin = 0.0
y_origin = 0.0
###########calculating the Volume per cell

# loop (i) through all subdirectories to read the MT3D001.UCN file
# and to determine M0, M1x, M1y, M2x and M2y for each time step and subfolder
# nr.sim = to calculate the number of simulations in the sto_dir
nr.sim = length(list.dirs(path = sto_subdirs,recursive = F))
# start.time is a variable to keep track of the running time
start.time = Sys.time()
for (i in 1: nr.sim)
{
  if (i < 10 ) setwd(paste0(sto_subdirs,"/00",i))
  if( i >= 10) setwd(paste0(sto_subdirs,"/0",i))
  if( i >= 100) setwd(paste0(sto_subdirs,"/",i))
  trans_conc = readucn("MT3D001",NLAY = 1,NTTS = NTTS )
  ##pts is a vector containing all times steps and required for the calculation of the transient mass moments.
  pts = unique((trans_conc$TIME))
  ## below the dimension of the transport model
  nrCols = max(trans_conc$COL)
  nrRows = max(trans_conc$ROW)
  
  M0 = Calc_M0(trans_conc)
  cat(paste("calculating M0 for run number ", i),"\n")
  M1x = Calc_M1_x(trans_conc,M0)
  cat(paste("calculating M1x for run number ", i),"\n")
  M1y = Calc_M1_y(trans_conc,M0_time = M0)
  cat(paste("calculating M1y for run number ", i),"\n")
  M2x = Calc_M2_x(trans_conc,M0,M1x)
  cat(paste("calculating M2x for run number ", i),"\n")
  M2y = Calc_M2_y(trans_conc,M0,M1y)
  cat(paste("calculating M2y for run number ", i),"\n")
  COVxy = Calc_M2_xy(trans_conc,M0_time = M0,M1_x_time = M1x,M1_y_time = M1y)
  cat(paste("calculating COVXY for run number ", i),"\n")
  stoch_MM_df = cbind(stoch_MM_df,M0,M1x,M1y,M2x,M2y,COVxy)
  one.calculation = Sys.time()-start.time
  print(paste("Current calculation time :",Sys.time()-start.time))
}

setwd(cur_dir) #go back to the base directory where this file is located
stoch_MM_df = cbind(pts,stoch_MM_df) #add the time steps to the data
stoch_MM_data =  data.frame(stoch_MM_df) #create a data.frame from this
#write.csv(file = "stoch_MM_data.csv",x = stoch_MM_data) #write the data file
write.csv(file = paste0(sto_dir,".csv"),x = stoch_MM_data) #write the data file
