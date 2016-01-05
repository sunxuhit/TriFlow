#!/bin/bash
date

#. ./hadd_all.sh /project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/List/Phi_list/flow_List/flow_SE.list


if [ $# -eq 1 ]
  then
    PID=Phi
    SM=_SE_
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/merged_Yields$SM${Energy}_"
    suffix=".root"
    counter=0
    for item in `cat $1`
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

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.err ./run_hadd_$PID$SM$counter.csh

      mv run_hadd_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
