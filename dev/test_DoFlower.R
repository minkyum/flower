library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)



###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
# numSite <- 6; cc <- 50



###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/flower/PLCF_Parameters.json')
source(params$setup$rFunctions)

# Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)


###############################
##
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks # the number of cells in a chunk: nearly 67,234 

# writeRaster(imgNum, filename=paste0(params$setup$outDir,strSite,'/base_image_num.tif'), format="GTiff", overwrite=TRUE)


# # Points
# pt1 <- shapefile('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/shp/Ansan_2023_16_36_inside flowering.shp')
# pt2 <- shapefile('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/shp/Ansan_2023_16_36_outside flowering.shp')
# pt3 <- shapefile('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/shp/Ansan_2023_16_36_outside non-flowering.shp')
# pt1 <- spTransform(pt1,crs(imgBase))
# pt2 <- spTransform(pt2,crs(imgBase))
# pt3 <- spTransform(pt3,crs(imgBase))
# 
# pixNums1 <- extract(imgNum,pt1)
# pixNums2 <- extract(imgNum,pt2)
# pixNums3 <- extract(imgNum,pt3)

pts <- shapefile('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/shp/seoul_flower.shp')
pts <- spTransform(pts,crs(imgBase))
pixNums <- extract(imgNum,pts)



# ###############################
# ##
# registerDoMC(params$setup$numCores)
# registerDoMC()
# 
# 
# # for(j in 1:2){
# #   if(j==1){
# #     pixNums <- c(pixNums1,pixNums2)
# #   }else{
# #     pixNums <- pixNums3
# #   }
#   
#   foreach(ppp=1:length(pixNums)) %dopar% {
#   # for(ppp in 1:length(pixNums)){
#     pixNum  <- pixNums[ppp]
#     
#     ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
#     file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)
#     
#     load(file)
#     
#     ##
#     blue  <- band1[pixNum%%chunk,]
#     green <- band2[pixNum%%chunk,]
#     red   <- band3[pixNum%%chunk,]
#     nir   <- band4[pixNum%%chunk,]
#     phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr
#     
#     blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000
#     
#     pheno_pars <- params$phenology_parameters
#     qa_pars    <- params$qa_parameters
#     
#     
#     i1 <- 2.5*(nir - red) / (nir + 2.4*red + 1)         # EVI2
#     i2 <- (nir - green) / (nir + green)                 # GNDVI
#     i3 <- red/nir                                       # Red/NIR      # for whiteness
#     i4 <- (red+green+blue)/((green/blue)*(red-blue+1))  # EBI          # for whiteness
#     
#     
#     subdates <- dates[dates > as.Date('2022-12-31') & dates < as.Date('2024-01-01')]
#     subvi1    <- i1[dates > as.Date('2022-12-31') & dates < as.Date('2024-01-01')]
#     subvi2    <- i2[dates > as.Date('2022-12-31') & dates < as.Date('2024-01-01')]
#     subvi3    <- i3[dates > as.Date('2022-12-31') & dates < as.Date('2024-01-01')]
#     subvi4    <- i4[dates > as.Date('2022-12-31') & dates < as.Date('2024-01-01')]
#     
#     # save data
#     # save(subdates,subvi1,subvi2,subvi3,subvi4,
#     #      file=paste0('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/dat/ts_airborn_',j,'_',sprintf('%03d',ppp),'.rda'))
#     save(subdates,subvi1,subvi2,subvi3,subvi4,
#          file=paste0('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/dat/ts_seoul_',sprintf('%03d',ppp),'.rda'))
#     
#     ########################################
#     setwd('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/fig')
#     # png(filename=paste0('ts_',j,'_',sprintf('%03d',ppp),'.png'),width=12,height=18,units='in',res=300)
#     png(filename=paste0('ts_seoul_',sprintf('%03d',ppp),'.png'),width=12,height=18,units='in',res=300)
#     
#     par(mfrow=c(4,1),oma=c(2,1,0,0),mar=c(4,5,1,2),mgp=c(2.8,1.2,0))
#     plot(subdates,subvi1,ylim=c(-0.1,0.85),pch=21,bg='darkgreen',cex=1.8,
#          xlab='',ylab='EVI2',cex.lab=2,cex.axis=1.5)
#     plot(subdates,subvi2,ylim=c(-0.1,0.85),pch=21,bg='green',cex=1.8,
#          xlab='',ylab='GNDVI',cex.lab=2,cex.axis=1.5)
#     plot(subdates,subvi3,ylim=c(-0.1,0.85),pch=21,bg='blue',cex=1.8,
#          xlab='',ylab='SR',cex.lab=2,cex.axis=1.5)
#     plot(subdates,subvi4,ylim=c(-0.1,0.85),pch=21,bg='red',cex=1.8,
#          xlab='',ylab='EBI',cex.lab=2,cex.axis=1.5)
#     
#     dev.off()
#     
#     print(ppp)
#   }
# # }
# 
# 
# # ########################################
# # path  <- '/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/dat'
# # file1 <- list.files(path=path,pattern=glob2rx('ts_airborn_1_*.rda'),full.names=T)
# # file2 <- list.files(path=path,pattern=glob2rx('ts_airborn_2_*.rda'),full.names=T)
# # 
# # load(file1[1])
# # 
# # dat1 <- matrix(NA,length(file1),length(subdates))
# # dat2 <- matrix(NA,length(file2),length(subdates))
# # for(i in 1:length(file1)){
# #   load(file1[i])
# #   dat1[i,] <- subvi4
# # }
# # for(i in 1:length(file2)){
# #   load(file2[i])
# #   dat2[i,] <- subvi4
# # }
# # 
# # d1 <- apply(dat1,2,median,na.rm=T)
# # d2 <- apply(dat2,2,median,na.rm=T)
# # plot(d1,col='red')
# # points(d2,col='blue')
# # 
# # dat <- cbind(subdates,d1)
# # dat <- na.omit(dat)
# # 
# # spl <- smooth.spline(dat[,1],dat[,2], spar=0.4)
# # smoothVI <- predict(spl,as.numeric(as.Date('2023-01-01')):as.numeric(as.Date('2023-12-31')))$y
# # # screen and fill values less than the the dormant value
# # xSmooth[xSmooth < dormant_value] <- dormant_value

  
  
  
  
  
  
########################################  

ppp <- 2

pixNum  <- pixNums[ppp]
  
ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)
  
##
blue  <- band1[pixNum%%chunk,]
green <- band2[pixNum%%chunk,]
red   <- band3[pixNum%%chunk,]
nir   <- band4[pixNum%%chunk,]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

# blue  <- band1[i,]
# green <- band2[i,]
# red   <- band3[i,]
# nir   <- band4[i,]


###
DoFlower <- function(blue, green, red, nir, dates, phenYrs, params){

  #### Despike, calculate dormant value, fill negative VI values with dormant value
  log <- try({

    pheno_pars <- params$phenology_parameters
    
    blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000

    vi  <- 2.5*(nir - red) / (nir + 2.4*red + 1)        # EVI2
    vi1 <- (red+green+blue)/((green/blue)*(red-blue+1)) # Enhanced Bloom Index
    vi2 <- red / nir                                    # Simple Ratio
    vi3 <- (nir - green) / (nir + green)                # GNDVI
    
    
    # Spikes check, and remove
    spikes <- CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars)
    blue[spikes] <- NA; green[spikes] <- NA; red[spikes] <- NA; nir[spikes] <- NA
    
    vi[spikes] <- NA; vi1[spikes] <- NA; vi2[spikes] <- NA; vi3[spikes] <- NA

    #
    splineStart <- as.Date(paste0(phenYrs,'-03-01'))
    numDaysFit  <- 180
    splineEnd   <- splineStart+(numDaysFit-1)
    numYrs <- length(phenYrs)

  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,length(phenYrs)))}


  outAll <- c()
  for(y in 1:numYrs){
    log <- try({

      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")

      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)
      dateSub   <- dates[dateRange]
      inWin     <- pred_dates > splineStart[y] + 30 & pred_dates < splineEnd[y] - 30
      
      viSub <- vi[dateRange]
      viSub1 <- vi1[dateRange] + vi2[dateRange] + (1-vi3[dateRange])
      
      
      #Get index of pixels with good values
      ind <- !is.na(viSub)  
      # smooth with a spline to get continuous daily series
      spline_fit  <- smooth.spline(dateSub[ind],  viSub[ind], spar=pheno_pars$splineSpar)
      spline_fit1 <- smooth.spline(dateSub[ind], viSub1[ind], spar=pheno_pars$splineSpar)
      
      # Smooted VI
      smoothed_vi  <- predict( spline_fit, as.numeric(pred_dates))$y
      smoothed_vi1 <- predict(spline_fit1, as.numeric(pred_dates))$y
      
      # Check vegetated or not
      amp <- max(smoothed_vi) - min(smoothed_vi)
      if (amp < 0.3) {outAll <- c(outAll,NA);next}
      
      # Find peak
      peaks <- FindPeaks(smoothed_vi1)
      if (all(is.na(peaks))) {outAll <- c(outAll,NA);next}
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi1, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,NA);next}
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inWin[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      if (length(full_segs)==0) {outAll <- c(outAll,NA);next}
      
      
      # Find the point of maximum increase date
      date_max_ebi <- full_segs[[1]][2]
      out <- date_max_ebi + 59
      
    },silent=TRUE) #End of the try block
    
    if(inherits(log, "try-error")){
      outAll <- c(outAll,NA)
    }else if(length(out)==0){
      outAll <- c(outAll,NA)
    }else{
      outAll <- c(outAll,out)
    }
  }

  return(outAll)
}

plot(smoothed_vi)
plot(smoothed_vi1)
abline(v=peaks)
