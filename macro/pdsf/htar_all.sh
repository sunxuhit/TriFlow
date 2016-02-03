#!/bin/bash
date

#. ./htar_all.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    List_SM=SE
    SM=_${List_SM}_
    Energy=200GeV
    OutDir="/home/x/xusun/AuAu$Energy/$PID/file_${Energy}_${PID}${SM}"
    suffix=".tar"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/backup_List/backup_${List_SM}.list"
    counter=0
    for item in `cat $InPutList`
    do
      cp ./run_hadd.csh ./run_htar_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID" >> run_htar_$PID$SM$counter.csh
      echo " " >> run_htar_$PID$SM$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_htar_$PID$SM$counter.csh
      echo " " >> run_htar_$PID$SM$counter.csh
      echo -n "htar $OutName " >> run_htar_$PID$SM$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_htar_$PID$SM$counter.csh
      done

      echo " " >> run_htar_$PID$SM$counter.csh
      echo " " >> run_htar_$PID$SM$counter.csh
      echo "echo 'This is the end of htar\!\!\!'" >> run_htar_$PID$SM$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/htar/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/htar/job_$PID$SM$counter.err ./run_htar_$PID$SM$counter.csh

      mv run_htar_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/htar/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
