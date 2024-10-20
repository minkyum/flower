library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 6



########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/flower/PLCF_Parameters.json')
source(params$setup$rFunctions)

phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_flower')
print(ckPheDir)

files <- list.files(path=ckPheDir,pattern=glob2rx('*.rda'),full.names=T)
print(length(files))

#
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks


if(length(files)==200){
   
  pheDir <- paste0(params$setup$outDir,strSite,'/phe')
  if (!dir.exists(pheDir)) {dir.create(pheDir)}
  
  Fmat <- matrix(NA,numPix,length(phenYrs))
  for(i in 1:numChunks){
    cc <- sprintf('%03d',i)
    cfile <- paste0(ckPheDir,'/chunk_flower_',cc,'.rda') 
    log <- try(load(cfile),silent=F)
    if (inherits(log, 'try-error')) next 
    
    if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
    }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
    chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
    Fmat[chunkStart:chunkEnd,] <- f_mat
  }
  
  for(yToDo in 1:length(phenYrs)){
    frast <- setValues(imgBase,(Fmat[,yToDo]))
    writeRaster(frast,filename=paste0(pheDir,'/flower_',phenYrs[yToDo],'.tif'), format="GTiff", overwrite=TRUE)
  }
  
}



  
