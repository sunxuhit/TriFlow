#!/bin/bash
date

#. ./Flow.sh

if [ $# -eq 0 ]
  then
#    counter=0
    Name="_PhiFlow_StringMelting_"
    suffix=".root"
    for((counter=0;counter<=9;counter=counter+1))
    do
      cp ./run.csh ./run$Name$counter.csh

      echo -n "root4star -b -q -x 'PhiFlow.C(" >> run$Name$counter.csh
###############################energy###################################
#      echo -n 0',' >> run$Name$counter.csh  # 7GeV
#      echo -n 1',' >> run$Name$counter.csh  # 11GeV
#      echo -n 2',' >> run$Name$counter.csh  # 19GeV
#      echo -n 3',' >> run$Name$counter.csh  # 27GeV
      echo -n 4',' >> run$Name$counter.csh  # 39GeV
#      echo -n 5',' >> run$Name$counter.csh  # 62GeV
#      echo -n 6',' >> run$Name$counter.csh  # 200GeV
###############################energy###################################


##############################AMPT_Mode#################################
#      echo -n 0',' >> run$Name$counter.csh  # Default
      echo -n 1',' >> run$Name$counter.csh  # String Melting
##############################AMPT_Mode#################################

################################List####################################
      echo -n $counter',' >> run$Name$counter.csh # List
################################List####################################


############################start_event#################################
      echo -n 0',' >> run$Name$counter.csh  # start_event
############################start_event#################################

#############################stop_event#################################
#      echo -n 1000000000',' >> run$Name$counter.csh  # stop_event
      echo -n 1024',' >> run$Name$counter.csh  # stop_event: test mode
#############################stop_event#################################

#############################Mixed_Event#################################
      echo  0')'"'" >> run$Name$counter.csh  # Same Event
#      echo  1')'"'" >> run$Name$counter.csh  # Mixed Event
#############################Mixed_Event#################################

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /project/projectdirs/star/xusun/OutPut/AMPT_StringMelting/Log/run$Name$counter.log -e /project/projectdirs/star/xusun/OutPut/AMPT_StringMelting/Log/run$Name$counter.err ./run$Name$counter.csh

      mv run$Name$counter.csh /project/projectdirs/star/xusun/OutPut/AMPT_StringMelting/Script/
    done

  else
    echo "Wrong number of parameters"
fi
