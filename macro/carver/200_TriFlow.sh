#!/bin/bash
date

#. ./200_TriFlow.sh

if [ $# -eq 0 ]
  then
    Name="_200GeV_Lambda_ME_"
    suffix=".root"
    List="/project/projectdirs/star/xusun/OutPut/AuAu200GeV/List/run_list/200GeV_"
    suffixlist=".list"
    for((counter=1;counter<=1030;counter=counter+1))
#    for((counter=111;counter<=122;counter=counter+1))
    do
      cp ./200_run.pbs ./run$Name$counter.pbs
      sed -i "s/my_job/my_job_$counter/g" ./run$Name$counter.pbs 

      echo -n "root4star -b -q -x 'TriFlow.C(" >> run$Name$counter.pbs
      echo -n '"'$List$counter$suffixlist'",' >> run$Name$counter.pbs
      echo -n $counter',' >> run$Name$counter.pbs
###############################mode###################################
#      echo -n 0',' >> run$Name$counter.pbs  # fill ReCenterPar mode
#      echo -n 1',' >> run$Name$counter.pbs  # ReCenter ShiftPar mode
#      echo -n 2',' >> run$Name$counter.pbs  # Charged Flow mode
#      echo -n 3',' >> run$Name$counter.pbs  # Pion and Kaon mode
#      echo -n 4',' >> run$Name$counter.pbs  # Proton and Yields mode
#      echo -n 5',' >> run$Name$counter.pbs  # Phi mode
      echo -n 6',' >> run$Name$counter.pbs  # Lambda mode
###############################mode###################################

###############################energy###################################
      echo  0')'"'" >> run$Name$counter.pbs  # 200GeV
#      echo  1')'"'" >> run$Name$counter.pbs  # 39GeV
#      echo  2')'"'" >> run$Name$counter.pbs  # 27GeV
###############################energy###################################

      qsub run$Name$counter.pbs

      mv run$Name$counter.pbs /project/projectdirs/star/xusun/OutPut/AuAu200GeV/Script/TriFlow/
    done

  else
    echo "Wrong number of parameters"
fi
