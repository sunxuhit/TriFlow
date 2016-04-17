#!/bin/bash
date

#. ./hadd_PiKP_final.sh 


if [ $# -eq 0 ]
  then
    PID=nSigmaPion
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Mass2_$PID/merged_file/Flow_${Energy}"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${PID}_final.list"
      cp ./run_hadd.csh ./run_hadd_final_${PID}.csh
      echo "cd ./AuAu$Energy/Mass2_$PID/merged_file/" >> run_hadd_final_${PID}.csh
      echo " " >> run_hadd_final_${PID}.csh

      OutName=$OutDir$suffix

      echo "rm $OutName " >> run_hadd_final_${PID}.csh
      echo " " >> run_hadd_final_${PID}.csh
      echo -n "hadd $OutName " >> run_hadd_final_${PID}.csh

      for item in `cat $InPutList`
      do
	echo -n "$item " >> run_hadd_final_${PID}.csh
      done

      echo " " >> run_hadd_final_${PID}.csh
      echo " " >> run_hadd_final_${PID}.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_final_${PID}.csh

      qsub -hard -l scratchfree=500,h_cpu=6:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_final_${PID}.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_final_${PID}.err ./run_hadd_final_${PID}.csh

      mv run_hadd_final_${PID}.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
    else
      echo "Wrong number of parameters"
fi
