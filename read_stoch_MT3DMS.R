###############################################################################
# R script to calculate the MM from a stochastic run of MT3DMS
# NOT in a markdown because there are currently some issues with markdown/knitr
# when in a code chunk the directory is changed for reading all *.UCN files
# G.Bier January, 2024
###############################################################################
rm(list = ls())
source("MMfunctions.R")

# First read all *.UCN files in the STOCHASTIC MT3DMS-FOLDER.
# For this give the folder name (case sensitive) in de "string"variable `sto_dir`.
# Next the UCN files will be read in the subsequent folders (001,002,..010)
sto_dir = "stoch_mm/test_stoch_ass2_MT3DMS"
cur_dir = getwd()

setwd(paste0(cur_dir,"/",sto_dir))
print(paste("Folder to read UCN files : ",getwd()))
sto_subdirs = getwd()

stoch_MM_df = c() #data.frame(colnames("time","stoch.nr","M0","M1x","My1","M2x","M2y","COVxy")) #"time","stoch.nr","M0","M1x","My1","M2x","M2y","COVxy")

###########calculating the Volume per cell
delx = 1.0 #10. 
dely = 1.0 #10.
delz = 1.0 #10.
porosity =0.4
volume = delx*dely*delz*porosity
x_origin = 0.0
y_origin = 0.0
###########calculating the Volume per cell


for (i in 1: 10)
{
  if (i < 10 ) setwd(paste0(sto_subdirs,"/00",i))
  if( i == 10) setwd(paste0(sto_subdirs,"/010"))
  trans_conc = readucn("MT3D001",NLAY = 1,NTTS = 74 )
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
}

setwd(cur_dir)
stoch_MM_df = cbind(pts,stoch_MM_df)
stoch_MM_data =  data.frame(stoch_MM_df)
write.csv(file = "stoch_MM_data.csv",x = stoch_MM_data)

# M0 = Calc_M0(trans_conc)
# M1x = Calc_M1_x(trans_conc,M0)
# M1y = Calc_M1_y(trans_conc,M0_time = M0)
# M2x = Calc_M2_x(trans_conc,M0,M1x)
# M2y = Calc_M2_y(trans_conc,M0,M1y)
# COVxy = Calc_M2_xy(trans_conc,M0_time = M0,M1_x_time = M1x,M1_y_time = M1y)



# library("ellipse")
# 
# for (i in 1:length(pts))#requires gifski and ellipse packages!!
# {
#   corxy = COVxy[i]/(M2x[i] * M2y[i])
#   M  = matrix(c(1,corxy,corxy,1),ncol = 2, byrow=TRUE) #the correlation matrix seems to be required for ellipse
#   
#   plot(M1x[1:i],M1y[1:i],type="l",xlim=c(2,40),asp=1,col="blue",lwd=2,
#        main=paste("time=",pts[i]),xlab="x",ylab="y")
#   points(M1x[i],M1y[i],pch=20,col="blue",cex=2)
#   polygon(ellipse(M,scale=c(M2x[i],M2y[i]),centre=c(M1x[i],M1y[i]),level=0.68,fill=TRUE),col=rgb(0,0,1,0.2),border=NA)
#   grid(col="black")
# }


#go back to the original directory

