#!/bin/bash
date

#. ./hadd_all_2nd.sh 


if [ $# -eq 0 ]
  then
    PID=Phi
    SM=_SE_
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/MERGED_Yields$SM${Energy}_"
    suffix=".root"
    InPutList="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/List/${PID}_list/flow_List/flow_${SM}_2nd.list"
    counter=0
    for item in `cat $1`
    do
      cp ./run_hadd.csh ./run_hadd_2nd_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID/flow_$PID/merged_file/" >> run_hadd_2nd_$PID$SM$counter.csh
      echo " " >> run_hadd_2nd_$PID$SM$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_2nd_$PID$SM$counter.csh
      echo " " >> run_hadd_2nd_$PID$SM$counter.csh
      echo -n "hadd $OutName " >> run_hadd_2nd_$PID$SM$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_hadd_2nd_$PID$SM$counter.csh
      done

      echo " " >> run_hadd_2nd_$PID$SM$counter.csh
      echo " " >> run_hadd_2nd_$PID$SM$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_2nd_$PID$SM$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_$PID$SM$counter.err ./run_hadd_2nd_$PID$SM$counter.csh

      mv run_hadd_2nd_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
