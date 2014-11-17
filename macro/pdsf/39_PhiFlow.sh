#!/bin/bash
date

#. ./39_PhiFlow.sh

if [ $# -eq 0 ]
  then
#    counter=0
    Name="_39GeV_Lambda_ME_"
    suffix=".root"
    for((counter=0;counter<=31;counter=counter+1))
#    for((counter=23;counter<=28;counter=counter+1)) #test
    do
      cp ./run.csh ./run$Name$counter.csh

      echo -n "root4star -b -q -x 'PhiFlow.C(" >> run$Name$counter.csh
###############################energy###################################
#      echo -n 0',' >> run$Name$counter.csh  # 200GeV
      echo -n 1',' >> run$Name$counter.csh  # 39GeV
#      echo -n 2',' >> run$Name$counter.csh  # 27GeV
###############################energy###################################

###############################X_flag###################################
#      echo -n 0',' >> run$Name$counter.csh # Same Event 
      echo -n 1',' >> run$Name$counter.csh  # Mixed Event
###############################X_flag###################################

      echo -n $counter',' >> run$Name$counter.csh # List

############################start_event#################################
      echo -n 0',' >> run$Name$counter.csh  # start_event
############################start_event#################################

#############################stop_event#################################
      echo -n 1000000000',' >> run$Name$counter.csh  # stop_event
#      echo -n 100024',' >> run$Name$counter.csh  # stop_event: test mode
#############################stop_event#################################

##############################Partilce##################################
#      echo 0')'"'" >> run$Name$counter.csh  # phi meson
      echo 1')'"'" >> run$Name$counter.csh  # Lambda
#      echo 2')'"'" >> run$Name$counter.csh  # anti-Lambda
#      echo 3')'"'" >> run$Name$counter.csh  # K0s
##############################Partilce##################################

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /project/projectdirs/star/xusun/OutPut/AuAu39GeV/Log/Lambda/flow_Lambda/job$Name$counter.log -e /project/projectdirs/star/xusun/OutPut/AuAu39GeV/Log/Lambda/flow_Lambda/job$Name$counter.err ./run$Name$counter.csh

      mv run$Name$counter.csh /project/projectdirs/star/xusun/OutPut/AuAu39GeV/Script/PhiFlow/
    done

  else
    echo "Wrong number of parameters"
fi
