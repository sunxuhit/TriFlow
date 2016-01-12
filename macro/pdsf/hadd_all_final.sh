#!/bin/bash
date

#. ./hadd_all_final.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    List_SM=SE
    SM=_${List_SM}_
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/MERGED_Yields$SM${Energy}_"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${List_SM}_final.list"
    counter=0
      cp ./run_hadd.csh ./run_hadd_final_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID/flow_$PID/merged_file/" >> run_hadd_final_$PID$SM$counter.csh
      echo " " >> run_hadd_final_$PID$SM$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_final_$PID$SM$counter.csh
      echo " " >> run_hadd_final_$PID$SM$counter.csh
      echo -n "hadd $OutName " >> run_hadd_final_$PID$SM$counter.csh

      for item in `cat $InPutList`
      do
	echo -n "$item" >> run_hadd_final_$PID$SM$counter.csh

	echo " " >> run_hadd_final_$PID$SM$counter.csh
	echo " " >> run_hadd_final_$PID$SM$counter.csh
	echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_final_$PID$SM$counter.csh
      done

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.err ./run_hadd_final_$PID$SM$counter.csh

      mv run_hadd_final_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"

  else
    echo "Wrong number of parameters"
fi
