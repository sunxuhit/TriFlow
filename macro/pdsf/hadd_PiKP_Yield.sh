#!/bin/bash
date

#. ./hadd_PiKP_Yield.sh 


if [ $# -eq 0 ]
  then
    PID=nSigmaPion
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Mass2_$PID/Yields/merged_file/merged_yield_${Energy}_"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/yield_List/yield_${PID}.list"
    counter=0
    for item in `cat $InPutList`
    do
      cp ./run_hadd.csh ./run_hadd_${PID}_Yields_$counter.csh
      echo "cd ./AuAu$Energy/Mass2_$PID/Yields" >> run_hadd_${PID}_Yields_$counter.csh
      echo " " >> run_hadd_${PID}_Yields_$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_${PID}_Yields_$counter.csh
      echo " " >> run_hadd_${PID}_Yields_$counter.csh
      echo -n "hadd $OutName " >> run_hadd_${PID}_Yields_$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_hadd_${PID}_Yields_$counter.csh
      done

      echo " " >> run_hadd_${PID}_Yields_$counter.csh
      echo " " >> run_hadd_${PID}_Yields_$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_${PID}_Yields_$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=6:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${PID}_Yields_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${PID}_Yields_$counter.err ./run_hadd_${PID}_Yields_$counter.csh

      mv run_hadd_${PID}_Yields_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
