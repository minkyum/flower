###############################
library(rjson)

params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/flower/PLCF_Parameters.json')
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)

print(strSite)


###############################
## 03_Get flowering dates
setwd(paste0(params$setup$logDir,'03'))
for(numSite in c(3,4,8)){
  nn <- sprintf('%03d',numSite)
  for(cc in 1:params$setup$numChunks){
    system(paste('qsub -V -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_03.sh ', nn, cc, sep=''))
  }
}


## 04_Create map
setwd(paste0(params$setup$logDir,'04'))
for(numSite in c(3,4,8)){
  system(paste('qsub -V -pe omp 2 -l h_rt=02:00:00 ',params$setup$rScripts,'run_script_04.sh ',numSite,sep=''))  
}



