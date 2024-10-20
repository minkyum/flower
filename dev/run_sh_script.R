###############################
setwd('/projectnb/modislc/users/mkmoon/sukyungkim/flowering_planet/tsplot/seoul/')
for(numSite in c(6)){
  system(paste('qsub -V -pe omp 28 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/flower/dev/run_script_test_DoFlower.sh ', numSite, sep=''))
}
