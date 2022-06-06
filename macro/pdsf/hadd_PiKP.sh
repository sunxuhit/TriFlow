#!/bin/bash
date

#. ./hadd_PiKP.sh 


if [ $# -eq 0 ]
  then
    PID=Proton
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Mass2_$PID/merged_file/merged_flow_${Energy}_"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${PID}.list"
    counter=0
    for item in `cat $InPutList`
    do
      cp ./run_hadd.csh ./run_hadd_${PID}_$counter.csh
      echo "cd ./AuAu$Energy/Mass2_$PID" >> run_hadd_${PID}_$counter.csh
      echo " " >> run_hadd_${PID}_$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_${PID}_$counter.csh
      echo " " >> run_hadd_${PID}_$counter.csh
      echo -n "hadd $OutName " >> run_hadd_${PID}_$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_hadd_${PID}_$counter.csh
      done

      echo " " >> run_hadd_${PID}_$counter.csh
      echo " " >> run_hadd_${PID}_$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_${PID}_$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=6:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${PID}_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${PID}_$counter.err ./run_hadd_${PID}_$counter.csh

      mv run_hadd_${PID}_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
