#!/bin/bash
date

#. ./39_Resolution.sh 39_Resolution_list

if [ $# -eq 1 ]
  then
    counter=0
    Name="_39GeV_"
    for INPUTFILE in `cat $1`
    do
      echo $INPUTFILE
      cat $INPUTFILE | while read line; do
	let "counter+=1"
        job_start=`awk 'BEGIN {split("'"$line"'",arr);print arr[1]}'`
        job_stop=`awk 'BEGIN {split("'"$line"'",arr);print arr[2]}'`

	cp ./run.csh ./run$Name$counter.csh

	echo -n "root4star -b -q -x 'Resolution.C(" >> run$Name$counter.csh
	echo -n $job_start','  >> run$Name$counter.csh  
	echo -n $job_stop',' >> run$Name$counter.csh  
#	echo  0')'"'" >> run$Name$counter.csh  #200GeV
	echo  1')'"'" >> run$Name$counter.csh  #39GeV 
#	echo  2')'"'" >> run$Name$counter.csh  #27GeV 

	qsub ./run$Name$counter.csh

	mv run$Name$counter.csh /project/projectdirs/star/xusun/OutPut/AuAu39GeV/Script/PhiFlow/
      done
    done

  else
    echo "Wrong number of parameters"
fi
