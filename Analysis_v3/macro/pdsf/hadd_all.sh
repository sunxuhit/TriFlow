#!/bin/bash
date

#. ./hadd_all.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    List_SM=SE
    SM=_${List_SM}_
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/merged_Yields$SM${Energy}_"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${List_SM}.list"
    counter=0
    for item in `cat $InPutList`
    do
      cp ./run_hadd.csh ./run_hadd_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID/flow_$PID" >> run_hadd_$PID$SM$counter.csh
      echo " " >> run_hadd_$PID$SM$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_$PID$SM$counter.csh
      echo " " >> run_hadd_$PID$SM$counter.csh
      echo -n "hadd $OutName " >> run_hadd_$PID$SM$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_hadd_$PID$SM$counter.csh
      done

      echo " " >> run_hadd_$PID$SM$counter.csh
      echo " " >> run_hadd_$PID$SM$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_$PID$SM$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=6:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.err ./run_hadd_$PID$SM$counter.csh

      mv run_hadd_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
