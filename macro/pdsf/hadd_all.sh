#!/bin/bash
date

#. ./hadd_all.sh


if [ $# -eq 0 ]
  then
    PID=Phi
    SM=_ME_
    Energy=200GeV
    n_stpes=9
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/$PID/flow_$PID/merged_file/merged_Yileds$SM${Energy}_"
    suffix=".root"
#    for((counter=0;counter<=52;counter=counter+1))
    for((counter=0;counter<=1;counter=counter+1)) #test
    do
      cp ./run_hadd.csh ./run_hadd_$PID$SM$counter.csh
      echo "cd ./AuAu$Energy/$PID/flow_$PID" >> run_hadd_$PID$SM$counter.csh
      echo " " >> run_hadd_$PID$SM$counter.csh

      OutName=$OutDir$counter$suffix

      echo -n "hadd $OutName " >> run_hadd_$PID$SM$counter.csh

      Order=0
      for((i_loop=0;i_loop<=$n_stpes;i_loop=i_loop+1)) # loop for Yileds files
      do
	let "Order=$counter*10+$i_loop"
	echo -n " Yields${SM}${Energy}_$Order$suffix" >> run_hadd_$PID$SM$counter.csh
      done
      echo " " >> run_hadd_$PID$SM$counter.csh

      qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job$Name$counter.err ./run_hadd_$PID$SM$counter.csh

      mv run_hadd_$PID$SM$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
    done

  else
    echo "Wrong number of parameters"
fi
