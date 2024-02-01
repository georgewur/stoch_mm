#################################################################
##source file which contains functions to calculation of 
#various spatial central moments
#For the full description of the workings and the arguments
#look at the markdown
# G. Bier January 2024
#################################################################
library(magrittr)
#function to read UNC files
#' Read MT3DMS .ucn File
#'
#' This function reads in a binary ucn file and creates a data frame
#' composed of the following elements:
#' \describe{
#' \item{TRANS}{Transport Time Step}
#' \item{STP}{Flow Time Step}
#' \item{PER}{Stress Period}
#' \item{TIME}{Elapsed Time}
#' \item{LAY}{Model Layer}
#' \item{ROW}{Model Row}
#' \item{COL}{Model Column}
#' \item{CONC}{Concentration}
#' }
#' @param ucname This is the name of the ucn file
#' @param NLAY This is the number of layers assigned to the model
#' @param NTTS This is the number of transport time steps that are printed to the .ucn file. This can be obtained from the .btn file. 
#' This will cause an error if the value assigned is greater than the correct value. Future versions will need to be developed so that NTTS is not
#' needed for this function to operate properly.
#' @export
#' @examples
#' readucn("MT3DMS001S", NLAY = 8, NTTS = 22)

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
      print(paste("loading time step :",Q))
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

# function extract_sample_data()
# This function extracts the concentration time series of the selected (sample) locations.\
# It requires `sampler_rows` and `sampler_cols` to determine which values should be extracted.\
# It does not require any arguments.\
# The output is a data set containing time series of the concentration at the selected locations.
extract_sample_data = function()
{
  sample_ROW = filter(Trans_conc, ROW %in% sampler_rows)
  sample_ROW_COL = filter(sample_ROW, COL %in% sampler_cols)
  return(sample_ROW_COL)
}

####  calc_M0(data_set)
# This function calculates the total mass of the distribution series of `data_set`.\
# The argument: `data_set` could be the whole data set: `Trans_conc` or only the sampler data: `Sample_conc`\
# The output is a time series of the masses of the `data_set`.
Calc_M0 = function(data_set)
{
  M0_time = c()
  for (i in 1:length(pts))
  {
    time_step = pts[i]
    M0 = volume*sum(subset(data_set, TIME == time_step, select = CONC))
    M0_time = c(M0_time,M0)
  }
  return(M0_time)
}

####  calc_M1_x(data_set, M0_time, smpl = FALSE)
# This function calculates the first mass moment (mass midpoint) in the x-direction in time of the `data_set`.\
# `M0_time` is the time series of the mass of the selected `data_set`.\
# `smpl = FALSE` is Boolean indicating whether the `data_set` contains only data for the sampler locations (`TRUE`). Default is set to `FALSE`, expecting the whole concentration distribution.
Calc_M1_x = function(data_set, M0_time,smpl = FALSE)
{
  if (smpl == TRUE)
  {
    nrCols = length(sampler_cols)
    col_numbers = sampler_cols
  }else{
    col_numbers = rep(1:nrCols) # the actual column numbers to process
  }
  M1_pts = c()
  for (i in 1:length(pts))
  {
    M1_cols = c()
    time_step = pts[i] #cycle through all time steps
    for (c in 1:nrCols) #cycle through all columns
    {
      Sum_C_col = sum(subset(data_set, TIME == time_step & COL == col_numbers[c], select = CONC))
      #print(Sum_C_col)
      M1_c = volume*(x_origin + col_numbers[c]*delx - delx/2)*Sum_C_col
      M1_cols = c(M1_cols,M1_c)
      #print(M1_c)
    }
    M1_x = 1/M0_time[i]*sum(M1_cols)
   # cat(paste("M1_x for time step :", i),'\n')
    #print(M1_x)
    M1_pts = c(M1_pts,M1_x)
  }
  return(M1_pts)
}


####  calc_M2_x(data_set, M0_time, M1_x\_time, smpl = FALSE)
# This function calculates the second central mass moment (spreading of mass as VAR(x)) in the x-direction in time of the `data_set`.\
# `M0_time` is the time series of the mass of the selected `data_set`.\
# `M1_x_time` is the first mass moment, mass midpoint, in the x-direction which required to determine the standard deviation (spread) of the plume. `smpl = FALSE` is Boolean indicating wheter the `data_set` contains only data for the sampler locations (`TRUE`). Default is set to `FALSE`, expecting the whole concentration distribution.
Calc_M2_x = function(data_set, M0_time, M1_x_time, smpl = FALSE)
{
  if (smpl == TRUE)
  {
    nrCols = length(sampler_cols)
    col_numbers = sampler_cols
  }else{
    col_numbers = rep(1:nrCols)
  }
  M2_pts = c()
  for (i in 1:length(pts))
  {
    M2_cols = c()
    time_step = pts[i]
    for (c in 1:nrCols)
    {
      Sum_C_col = sum(subset(data_set, TIME == time_step & COL == col_numbers[c], select = CONC))
      M2_c = volume*Sum_C_col*((x_origin + col_numbers[c]*delx - delx/2) - M1_x_time[i])^2
      M2_cols = c(M2_cols,M2_c)
    }
    #M2_x = sqrt(1/M0_time[i]*sum(M2_cols)) #this would calculate sqrt(VAR(x))
    M2_x = 1/M0_time[i]*sum(M2_cols)
    #print(i)
    M2_pts = c(M2_pts,M2_x)
  }
  return(M2_pts)
}


####  calc_M1_y(data_set, M0_time, smpl = FALSE)
# This function calculates the first mass moment (mass midpoint) in the y-direction in time of the `data_set`.\
# `M0_time` is the time series of the mass of the selected `data_set`.\
# `smpl = FALSE` is Boolean indicating wheter the `data_set` contains only data for the sampler locations (`TRUE`). Default is set to `FALSE`, expecting the whole concentration distribution.
Calc_M1_y = function(data_set, M0_time, smpl = FALSE)
  #Be aware that row numbering starts at the top and not at the bottom
{
  if (smpl == TRUE)
  {
    loc_nrRows = length(sampler_rows)
    row_numbers = sampler_rows
  }else{
    row_numbers = rep(1:nrRows)
    loc_nrRows = nrRows
  }
  M1_pts = c()
  for (i in 1:length(pts))
  {
    M1_rows = c()
    time_step = pts[i]  
    for (r in 1:loc_nrRows)
    {
      Sum_C_row = sum(subset(data_set, TIME == time_step & ROW == row_numbers[r], select = CONC))
      #      M1_r = volume*(y_origin + row_numbers[r] - dely/2) * Sum_C_row
      M1_r = volume*(y_origin + (nrRows - row_numbers[r])*dely + dely/2) * Sum_C_row
      M1_rows = c(M1_rows,M1_r)
    }
    M1_y = 1/M0_time[i]*sum(M1_rows)
    #print(i)
    M1_pts = c(M1_pts,M1_y)
  }
  return(M1_pts)
}

#### calc_M2_y(data_set, M0_time, M1_y\_time, smpl = FALSE)
# This function calculates the second central mass moment (spreading of mass as VAR(y)) in the y-direction in time of the `data_set`.\
# `M0_time` is the time series of the mass of the selected `data_set`.\
# `M1_y_time` is the first mass moment, mass midpoint, in the y-direction which required to determine the standard deviation (spread) of the plume. `smpl = FALSE` is Boolean indicating wheter the `data_set` contains only data for the sampler locations (`TRUE`). Default is set to `FALSE`, expecting the whole concentration distribution. 
Calc_M2_y = function(data_set,M0_time,M1_y_time, smpl=FALSE)
{
  if (smpl == TRUE)
  {
    loc_nrRows = length(sampler_rows)
    row_numbers = sampler_rows
  }else{
    row_numbers = rep(1:nrRows)
    loc_nrRows = nrRows
  }
  M2_pts = c()
  for (i in 1:length(pts))
  {
    M2_rows = c()
    time_step = pts[i]  
    for (r in 1:loc_nrRows)
    {
      Sum_C_row = sum(subset(data_set, TIME == time_step & ROW == row_numbers[r], select = CONC))
      # M2_r = volume*Sum_C_row*((y_origin + row_numbers[r] - dely/2)- M1_y_time[i])^2
      M2_r = volume*Sum_C_row*((y_origin + (nrRows - row_numbers[r])*dely + dely/2)- M1_y_time[i])^2
      M2_rows = c(M2_rows,M2_r)
    }
    #M2_y = sqrt(1/M0_time[i]*sum(M2_rows)) THIS WOULD NOT CALCULATE VAR(Y) BUT SQRT(VAR(Y))
    M2_y = 1/M0_time[i]*sum(M2_rows)
    #print(i)
    M2_pts = c(M2_pts,M2_y)
  }
  return(M2_pts)
}

#### Calc_M2_xy(data_set, M0_time, M1_x_time,M1_y_time, smpl = FALSE)
# This function calculated the covariance between x and y in time of `data_set` 
# `M0_time` is the time series of the mass of the selected `data_set`.\
# `M1_y_time` is the first mass moment, mass midpoint, in the y-direction which required to determine the standard deviation (spread) of the plume. `smpl = FALSE` is Boolean indicating wheter the `data_set` contains only data for the sampler locations (`TRUE`). Default is set to `FALSE`, expecting the whole concentration distribution. 
Calc_M2_xy = function(data_set, M0_time, M1_x_time, M1_y_time, smpl = FALSE)
  #to calculate the covariance the procedure below cycles through the grid in a column wise fashion
  #at every cell of the grid, the concentration is retrieved from the data set
  #and is the position in x and y known from which the covariance can be calculated
  #Cov(x,y) = Volume*C_row_col*(row -avg_row)*(col - avg_col)
{
  if (smpl == TRUE)
  {
    nrCols = length(sampler_cols)
    col_numbers = sampler_cols
    nrRows = length(sampler_rows)
    row_numbers = sampler_rows
  }else{
    col_numbers = rep(1:nrCols)
    row_numbers = rep(1:nrRows)
  }
  M2_xy_pts = c()
  
  for (timestep in 1:length(pts))
  {
    #select all Ci from one time step
    C_timestep = subset(data_set, TIME == pts[timestep], select = CONC)
    #print(cat("time:",timestep," conc: ",C_timestep$CONC))
    # Sum_M2_xy is a container for M2_xy for this time step
    Sum_M2_xy = 0
    for (row in 1:nrRows)#row numbers starts at the top!
    {
      for (col in 1:nrCols)
      {
        y_i = y_origin + (nrRows - row_numbers[row])*dely + dely/2
        x_i = x_origin + (col_numbers[col]*delx - delx/2)
        C_i = C_timestep$CONC[(row-1) * nrCols + col]
        M2_xy = volume * C_i * ((x_i - M1_x_time[timestep]) * (y_i - M1_y_time[timestep]))
        #print(cat('xi ',x_i,"yi",y_i,"Ci",C_i,"M2_xy",M2_xy))
        #print(cat('col ',col,"row",row,"Ci",C_i,"M2_xy",M2_xy))
        Sum_M2_xy = Sum_M2_xy + M2_xy
      }
    }
    #a simple check if counters row and col results in the same number to C_timestep
    if (length(C_timestep$CONC) != nrRows * nrCols)
    {
      print("Something went wrong determining COV(x,y)")
      #break
    }
    COVXY = 1/M0_time[timestep] * Sum_M2_xy
    #print(cat("time step :",timestep, "COV(x,y) :", COVXY))
    M2_xy_pts = c(M2_xy_pts,COVXY)
  }
  return(M2_xy_pts)
}

