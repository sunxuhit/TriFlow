#!/bin/bash
date

#. ./hadd_PiKP_final.sh 


if [ $# -eq 0 ]
  then
    PID=Proton
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Mass2_$PID/merged_file/Flow_${Energy}"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${PID}_final.list"
      cp ./run_hadd.csh ./run_hadd_final_${PID}_$counter.csh
      echo "cd ./AuAu$Energy/Mass2_$PID/merged_file/" >> run_hadd_final_${PID}_$counter.csh
      echo " " >> run_hadd_final_${PID}_$counter.csh

      OutName=$OutDir$suffix

      echo "rm $OutName " >> run_hadd_final_${PID}_$counter.csh
      echo " " >> run_hadd_final_${PID}_$counter.csh
      echo -n "hadd $OutName " >> run_hadd_final_${PID}_$counter.csh

      for item in `cat $InPutList`
      do
	echo -n "$item " >> run_hadd_final_${PID}_$counter.csh
      done

      echo " " >> run_hadd_final_${PID}_$counter.csh
      echo " " >> run_hadd_final_${PID}_$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_final_${PID}_$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=2:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_final_${PID}_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_final_${PID}_$counter.err ./run_hadd_final_${PID}_$counter.csh

      mv run_hadd_final_${PID}_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
    else
      echo "Wrong number of parameters"
fi
