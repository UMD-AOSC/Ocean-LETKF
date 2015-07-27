#!/bin/bash

## Setting up to run:  
#On a disk with large memory storage:
#Create an output directory, e.g. "/lustre/uid/OUTPUT/tmp_XYZ"
#Link the output directory to a local tmp directory, e.g. "ln -fs /lustre/uid/OUTPUT/tmp_XYZ tmp"
#Create a log directory, e.g. "/lustre/uid/OUTPUT/tmp_XYZ/log"
#Line the log directory to a local log directory, e.g. "ln -fs /lustre/uid/OUTPUT/tmp_XYZ/log log"
#Put an initial experiment directory in tmp/
#  e.g. "cp -rl /lustre/uid/OUTPUT/tmp_ABC/1990122700 tmp/"
#  note: it must contain all model RESTART directories (for each ensemble member) and letkf and/or hybrid analyses
#Also put the static files in the tmp directory
#  e.g. "cp -rl /lustre/uid/OUPUT/tmp_ABC/INIT/INPUT tmp/INIT/INPUT"
#Install rocoto if not already installed.
#Update the submit_job file with the correct directory and file names
#Update the xml script with the local drectories
#Run ./submit_job
#update the crontab with the calling line from submit_job, and the desired time interval: "crontab -e"

# Setup the 
CWD=`pwd` #/autofs/mnt/ncrc-svm1_home1/Steve.Penny/letkf/Ocean-LETKF
cd ..
root=`pwd`
cd $CWD
RUN=$root/run
LOG=$RUN/log

mkdir -p $LOG/workflow
mkdir -p $LOG/model_prep
mkdir -p $LOG/model
mkdir -p $LOG/letkf_prep
mkdir -p $LOG/letkf
