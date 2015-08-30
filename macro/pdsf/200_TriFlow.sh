#!/bin/bash
date

#. ./200_TriFlow.sh

if [ $# -eq 0 ]
  then
    Name="_200GeV_Lambda0_SE_"
    suffix=".root"
    List="/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/List/run_list/200GeV_"
    suffixlist=".list"
    for((counter=1;counter<=1030;counter=counter+1))
#    for((counter=111;counter<=155;counter=counter+1))
    do
      cp ./run.csh ./run$Name$counter.csh

      echo -n "root4star -b -q -x 'TriFlow.C(" >> run$Name$counter.csh
      echo -n '"'$List$counter$suffixlist'",' >> run$Name$counter.csh
      echo -n $counter',' >> run$Name$counter.csh
###############################mode###################################
#      echo -n 0',' >> run$Name$counter.csh  # fill ReCenterPar mode
#      echo -n 1',' >> run$Name$counter.csh  # ReCenter ShiftPar mode
#      echo -n 2',' >> run$Name$counter.csh  # Charged Flow mode
#      echo -n 3',' >> run$Name$counter.csh  # Pion and Kaon mode
#      echo -n 4',' >> run$Name$counter.csh  # Proton and Yields mode
#      echo -n 5',' >> run$Name$counter.csh  # Phi mode
      echo -n 6',' >> run$Name$counter.csh  # Lambda mode
###############################mode###################################

###############################energy###################################
      echo -n 0',' >> run$Name$counter.csh  # 200GeV
#      echo -n 1',' >> run$Name$counter.csh  # 39GeV
#      echo -n 2',' >> run$Name$counter.csh  # 27GeV
###############################energy###################################

###############################flag_ME###################################
      echo -n 0')' >> run$Name$counter.csh  # Same Event
#      echo -n 1')' >> run$Name$counter.csh  # Mixed Event
###############################flag_ME###################################
      echo -n "' > /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/Log/Correction/run" >> run$Name$counter.csh
      echo -n $Name$counter >> run$Name$counter.csh
      echo ".log" >> run$Name$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu39GeV/Log/Correction/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu39GeV/Log/Correction/job$Name$counter.err ./run$Name$counter.csh

      mv run$Name$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/Script/TriFlow/ 
    done

  else
    echo "Wrong number of parameters"
fi
