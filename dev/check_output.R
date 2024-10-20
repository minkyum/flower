library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 9



###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/flower/PLCF_Parameters.json')
source(params$setup$rFunctions)


## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/phe')
print(ckDir)

files <- list.files(path=ckDir,pattern=glob2rx('*.tif'),full.names=T)
print(length(files))

map <- raster(files[1])
poly <- shapefile('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/polygon/target_polygon/Robinia_Seoul_polygons.shp')
poly <- spTransform(poly,crs(map))

maps <- vector('list',length(files))
for(i in 1:length(files)){
  map       <- raster(files[i])
  maps[[i]] <- mask(map,poly)
}

########################################
setwd('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/fig')
# png(filename=paste0('ts_',j,'_',sprintf('%03d',ppp),'.png'),width=12,height=18,units='in',res=300)
png(filename=paste0('ts_seoul_',sprintf('%03d',ppp),'.png'),width=12,height=18,units='in',res=300)

par(mfrow=c(6,1),oma=c(2,1,0,0),mar=c(4,5,1,2),mgp=c(2.8,1.2,0))
for(i in 1:6){
  hist(maps[[i]])  
}

dev.off()

hist(map_mask)


  
