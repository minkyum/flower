#!/bin/bash

# module load python3/3.7.7
# module load gdal/3.1.2
# module load R 

echo Submitting $1
R --vanilla < /usr3/graduate/mkmoon/GitHub/flower/dev/test_DoFlower.R $1


