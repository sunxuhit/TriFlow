#!/bin/bash
date

#. ./hadd_all_final.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    List_SM=SE
    SM=_${List_SM}
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/Yields${SM}_${Energy}"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${List_SM}_final.list"
      cp ./run_hadd.csh ./run_hadd_final_$PID$SM.csh
      echo "cd ./AuAu$Energy/$PID/flow_$PID/merged_file/" >> run_hadd_final_$PID$SM.csh
      echo " " >> run_hadd_final_$PID$SM.csh

      OutName=$OutDir$suffix

      echo "rm $OutName " >> run_hadd_final_$PID$SM.csh
      echo " " >> run_hadd_final_$PID$SM.csh
      echo -n "hadd $OutName " >> run_hadd_final_$PID$SM.csh

      for item in `cat $InPutList`
      do
	echo -n "$item " >> run_hadd_final_$PID$SM.csh
      done
      echo " " >> run_hadd_final_$PID$SM.csh
      echo " " >> run_hadd_final_$PID$SM.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_final_$PID$SM.csh

      qsub -hard -l scratchfree=500,h_cpu=2:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM.err ./run_hadd_final_$PID$SM.csh

      mv run_hadd_final_$PID$SM.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
  else
    echo "Wrong number of parameters"
fi
