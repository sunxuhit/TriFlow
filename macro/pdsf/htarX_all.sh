#!/bin/bash
date

#. ./htarX_all.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    List_SM=SE
    SM=_${List_SM}_
    Energy=200GeV
    InPutDir="/home/x/xusun/AuAu$Energy/$PID/file_${Energy}_${PID}${SM}"
    suffix=".tar"
    for((counter=0;counter<=2;counter=counter+1)) #test
    do
      cp ./run_hadd.csh ./run_htarX_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID" >> run_htarX_$PID$SM$counter.csh
      echo " " >> run_htarX_$PID$SM$counter.csh

      InPutName=$InPutName$counter$suffix

      echo -n "htar -xf $InPutName " >> run_htarX_$PID$SM$counter.csh

      echo " " >> run_htarX_$PID$SM$counter.csh
      echo " " >> run_htarX_$PID$SM$counter.csh
      echo "echo 'This is the end of htarX\!\!\!'" >> run_htarX_$PID$SM$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/htar/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/htar/job_$PID$SM$counter.err ./run_htarX_$PID$SM$counter.csh

      mv run_htarX_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/htar/
    done

  else
    echo "Wrong number of parameters"
fi
