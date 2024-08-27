###############################################################################
# R script to calculate the (Central) Mass Moments CMM from a stochastic run of MT3DMS
# NOT in a markdown because there are currently some issues with markdown/knitr
# when in a code chunk the directory is changed for reading all *.UCN files
# G.Bier February-August, 2024
###############################################################################
rm(list = ls())
library(dplyr) # using the %>% function in readucn()
library(raster) # for creating a grid again from the ucn data
# Before running this code assign the proper folder at which the subfolders of
# the stochastic run is located (i.e. folder numbered 001,002,...010)
# assign this folder "between quotes" to variable sto_dir.
# Next the UCN files will be read from the subsequent folders (001,002,..010)
#
# Assign the number of time steps which need to be considered. 
# It can be found in the sub panel in the explorer window (lower left on screen) of GMS
# when a solute/species is selected 
# This number should be assigned to variable NTTS
# If NTTS is larger than the actual time steps the program will 
# throw an error regarding the readBin() function

# the main directory where the tochastic results can be found in the 001,002 etc. folders
sto_dir = "stoch_mm/stochastic_ass2_MT3DMS"
# the number of transport times steps
NTTS = 22

cur_dir = getwd()  #get currrent directory to go back to this after running

setwd(paste0(cur_dir,"/",sto_dir,"/")) #setting the working directory for reading the subsequent folders
print(paste("Folder to read UCN files : ",getwd()))
sto_subdirs = getwd()

stoch_MM_df = c() # open a data container for the results

###########calculating the Volume per cell
dx = 1.0#0.25#1.0 #10. 
dy =1.0#0.25#1.0 #10.
dz =1.0#0.25#1.0 #10.
porosity =0.4
vol = dx*dy*dz*porosity
###########calculating the Volume per cell


#####function for reading the binary UCN file containing concentrations(x,y,z,t)
readucn <- function(ucrootname, NLAY, NTTS){
  # NLAY IS THE NUMBER OF LAYERS IN THE MODEL
  # NTTS IS THE NUMBER OF TRANSPORT TIME STEPS IN THE MODEL
  ucname <- paste(ucrootname, ".ucn", sep = "")    
  to.read <- file(ucname, "rb")    
  TRANS <- c()
  STP <- c()
  PER <- c()
  TIME <- c()
  TEXT <- c()
  LAY <- c()
  CONC <- c()
  dat <- c()
  readblock <- function(){
    TRANS <- readBin(to.read, integer(), n = 1)    
    STP <- readBin(to.read, integer(), n = 1)
    PER <- readBin(to.read, integer(), n = 1)
    TIME <- readBin(to.read, double(), size = 4, n = 1)
    TEXT <- readChar(to.read, 16)
    NC <- readBin(to.read, integer(), n = 1)
    NR <- readBin(to.read, integer(), n = 1)
    #print(NC) #GB
    #print(NR) #GB
    LAY <- readBin(to.read, integer(), n = 1)
    CONC <- readBin(to.read, double(), size = 4, n = NR * NC, endian = "little")
    out <- list(TRANS, STP, PER, TIME, TEXT, NC, NR, LAY, CONC)
    return(out)            
  }    
  for(Q in 1:NTTS){
    for(K in 1:NLAY){
      dat[[length(dat) + 1]] <- readblock()
      #print(paste("loading time step :",Q))
    }
  }
  close(to.read)
  TRANS <- sapply(dat, "[[", 1)
  STP <- sapply(dat, "[[", 2)
  PER <- sapply(dat, "[[", 3)
  TIME <- sapply(dat, "[[", 4)
  LAY <- sapply(dat, "[[", 8)
  NC <- dat[[1]][6] %>% as.integer
  NR <- dat[[1]][7] %>% as.integer
  CONC <- sapply(dat, "[[", 9)
  # UCN <- tibble::data_frame( "data_frame was deprecated in tibble 1.1.0. GB-WUR 8/4/21
  UCN <- tibble::tibble(
    TRANS = rep(TRANS, each = (NC * NR)) %>% as.integer(),
    STP = rep(STP, each = (NC * NR)) %>% as.integer(),
    PER = rep(PER, each = (NC * NR)) %>% as.integer(), 
    TIME = rep(TIME, each = (NC * NR)) %>% as.double(),  
    LAY = rep(LAY, each = (NC * NR)) %>% as.integer(),
    ROW = rep(rep(rep(1:NR, each = NC), NLAY), NTTS) %>% as.integer(), 
    COL = rep(rep(rep(seq(1, NC, 1), NR), NLAY), NTTS) %>% as.integer(),
    CONC = CONC %>% as.double()
  )
  rm(TRANS)
  rm(STP)
  rm(CONC)
  rm(LAY)
  rm(NR)
  rm(NC)
  rm(TEXT)
  rm(TIME)
  rm(PER)    
  rm(dat)
  gc()
  return(UCN)        
}

###setting up some data containers for the CMM's


# loop (i) through all subdirectories to read the MT3D001.UCN file
# and to determine M0, M1x, M1y, M2x, M2y and COV(x,y) for each time step and subfolder
# nr.sim = to calculate the number of simulations in the sto_dir
nr.sim = length(list.dirs(path = sto_subdirs,recursive = F))
# start.time is a variable to keep track of the running time
start.time = Sys.time()
for (i in 1: nr.sim)
{
  if (i < 10 ) setwd(paste0(sto_subdirs,"/00",i))
  if( i >= 10 && i < 100 ) setwd(paste0(sto_subdirs,"/0",i))
  if( i >= 100) setwd(paste0(sto_subdirs,"/",i))
  ucn_data = readucn("MT3D001",NLAY = 1,NTTS = NTTS )
  ##time.steps is a vector containing all times steps and required for the calculation of the transient mass moments.
  time.steps = unique(ucn_data$TIME)  #which should be the same as NTTS when all stress periods are read
  ## below the dimension of the transport model
  nr.cols = max(ucn_data$COL)
  nr.rows = max(ucn_data$ROW)
  M0 = c()
  M1x = c()
  M1y = c()
  M2x = c()
  M2y = c()
  Covxy = c()
  ##loping over the data per time step
  for (t in 1: length(time.steps))
  {
    cur.t = time.steps[t] #numbering of the time steps
    conc_timestep = raster(nrow = nr.rows, ncol = nr.cols) #setup a raster with proper dimensions
    current.conc = ucn_data$CONC[ucn_data$TIME == cur.t] # assign(?) the concentration data for the current time step
    values(conc_timestep) = current.conc  #put the current concentrations in the raster
    
    ##calculating the total mass for this time step
    ##  BE AWARE that the origin of the coordinates of the model is in the lower left
    ## BUT the ROW NUMBERING STARTS from the UPPER ROW. 
    ## This applies for both the UCN row-column data as for the raster() function
    #
    #Calculating total mass by simply using the sum(function) for the raster with concentration X the volume
    M0 = c(M0,sum(current.conc)*vol)
    # since the C is allocated on a regular grid each row has one y-coordinate value
    # as well as each column will have one x-coordinate value, the sum of the concentrations can 
    # directly be calculated row-wise as well as column wise
    col_sums = apply(as.matrix(conc_timestep),MARGIN = 2, FUN = sum)
    row_sums = apply(as.matrix(conc_timestep),MARGIN = 1, FUN = sum)
    # For M1x; each column number is multiplied with dx-dx/2 to determine the center of the grid cells
    # this is then multiplied with the column concentrations in col_sums time the grid volume
    # for M1y the order of the row numbers is reversed to from nr_rows to 1 since the row numbering starts at the top 
    # of the model but the coordinates has its origin in the lower left corner
    currentM1x = 1/M0[t] * sum(((1:nr.cols)*dx -dx/2) * col_sums[1:nr.cols]) * vol
    currentM1y = 1/M0[t] * sum(((nr.rows:1)*dy-dy/2) * row_sums[1:nr.rows]) * vol
    M1x = c(M1x,currentM1x)
    M1y = c(M1y,currentM1y)
    # For M2x the previously calculated M1x is used to calculate (x_i - M1x)^2
    # For M2y like with M1y the numbering of rows starts from the top but the y-coordinate from the bottom
    # therefor the y-coordinate is calculated as (nr.row:1)*dy - dy/2
    currentM2x =  1/M0[t] * vol * sum(col_sums[1:nr.cols] * (((1:nr.cols)*dx - dx/2) - currentM1x)^2)
    M2x = c(M2x,currentM2x)
    currentM2y = 1/M0[t] * vol * sum(row_sums[1:nr.rows] * (((nr.rows:1)*dy - dy/2) - currentM1y)^2)
    M2y = c(M2y,currentM2y)
    sum_M2xy = 0.0
    
    # Here the covariance of concentration distribution is calculated
    # 1/M0 * vol * (x_i - M1x) * (y_i) - M1y) for each time step
    # It consists of a row loop and column loop. So in row order the data is 
    # evaluated. The rows are from the top to the bottom but the y-coordinate
    # starts at the bottom going up.Therefor the row numbering starts
    # from the last row which is the lower coordinate
    for (r in nr.rows:1)
    {
      for (c in 1:nr.cols)
      {
        currentM2xy = conc_timestep[r,c]*((c*dx - dx/2) - M1x[t]) * ((r*dy - dy/2) - M1y[t])
        sum_M2xy = sum_M2xy + currentM2xy
      }
    }
    current.M2xy = 1/M0[t] * vol * sum_M2xy
    Covxy = c(Covxy,current.M2xy)
  }
  
  stoch_MM_df = cbind(stoch_MM_df,M0,M1x,M1y,M2x,M2y,Covxy)
  ## some lines to track elapsed and remaining time
  avg.elapsed.per.run = as.numeric(difftime(Sys.time(),start.time,units = "secs")/i)
  print(paste("Current run :", i, "Elapsed time :",difftime(Sys.time(),start.time), " remaining time :", (nr.sim - i)*avg.elapsed.per.run/60," minutes"))
}


setwd(cur_dir) #go back to the base directory where this file is located
pts = time.steps
stoch_MM_df = cbind(pts,stoch_MM_df) #add the time steps to the data

stoch_MM_data =  data.frame(stoch_MM_df) #create a data.frame from this
#write.csv(file = "stoch_MM_data.csv",x = stoch_MM_data) #write the data file
write.csv(file = paste0(sto_dir,".csv"),x = stoch_MM_data) #write the data file
##writing some data of the models regarding their dimensions
dim.data = data.frame(dx,nr.cols,dy,nr.rows,dz,porosity)
write.table(file = paste0(sto_dir,".dim"),x = dim.data)
