#!/bin/bash
date

#. ./Flow.sh

if [ $# -eq 0 ]
  then
#    counter=0
    Name="_PhiFlow_Default_"
    suffix="_SE"
    for((counter=0;counter<=9;counter=counter+1))
    do
      cp ./run.csh ./run$Name$counter$suffix.csh

      echo -n "root4star -b -q -x 'PhiFlow.C(" >> run$Name$counter$suffix.csh
###############################energy###################################
#      echo -n 0',' >> run$Name$counter$suffix.csh  # 7GeV
#      echo -n 1',' >> run$Name$counter$suffix.csh  # 11GeV
#      echo -n 2',' >> run$Name$counter$suffix.csh  # 19GeV
#      echo -n 3',' >> run$Name$counter$suffix.csh  # 27GeV
      echo -n 4',' >> run$Name$counter$suffix.csh  # 39GeV
#      echo -n 5',' >> run$Name$counter$suffix.csh  # 62GeV
#      echo -n 6',' >> run$Name$counter$suffix.csh  # 200GeV
###############################energy###################################


##############################AMPT_Mode#################################
      echo -n 0',' >> run$Name$counter$suffix.csh  # Default
#      echo -n 1',' >> run$Name$counter$suffix.csh  # String Melting
##############################AMPT_Mode#################################

################################List####################################
      echo -n $counter',' >> run$Name$counter$suffix.csh # List
################################List####################################


############################start_event#################################
      echo -n 0',' >> run$Name$counter$suffix.csh  # start_event
############################start_event#################################

#############################stop_event#################################
#      echo -n 1000000000',' >> run$Name$counter$suffix.csh  # stop_event
      echo -n 1024',' >> run$Name$counter$suffix.csh  # stop_event: test mode
#############################stop_event#################################

#############################Mixed_Event#################################
      echo 0')'"'" >> run$Name$counter$suffix.csh  # Same Event
#      echo 1')'"'" >> run$Name$counter$suffix.csh  # Mixed Event
#############################Mixed_Event#################################

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /project/projectdirs/star/xusun/OutPut/AMPT_Default/Log/run$Name$counter$suffix.log -e /project/projectdirs/star/xusun/OutPut/AMPT_Default/Log/run$Name$counter$suffix.err ./run$Name$counter$suffix.csh

      mv run$Name$counter$suffix.csh /project/projectdirs/star/xusun/OutPut/AMPT_Default/Script/
    done

  else
    echo "Wrong number of parameters"
fi
