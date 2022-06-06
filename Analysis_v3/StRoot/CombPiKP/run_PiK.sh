#!/bin/bash
date

# . ./run_PiK.sh

if [ $# -eq 0 ] 
then 

    # mode
    # 0 = fit mode
    # 1 = test mode
    
    # mEnergy
    # 0 = 200GeV
    # 1 = 39GeV
    
    # mCharge
    # 0 = positive particles
    # 1 = negative particels

    # mCentrality
    # 0 = 0080
    # 1 = 0010
    # 2 = 1040
    # 3 = 4080
  
  for mode in {0..0}
  do
     for mEnergy in {0..0}
     do
        for mCharge in {0..1}
        do
	  for mCentrality in {0..0}
	  do
	    echo Send command to batch farm: root -l -b -q CombPID_Xu.C+\($mode,$mEnergy,$mCharge,$mCentrality\)

	    FILE="A"$mode"_B"$mEnergy"_C"$mCharge"_D"$mCentrality
	    echo FILE = $FILE
	    cp ./run.csh ./run$FILE.csh

	    echo -n "/u/huck/root/bin/root -l -b -q 'CombPID_Xu.C+(">>run$FILE.csh
	    echo -n $mode','$mEnergy','$mCharge','$mCentrality')' >>run$FILE.csh
	    echo -n "' > log/run_">>run$FILE.csh
	    echo -n $FILE>>run$FILE.csh
	    echo ".log">>run$FILE.csh

	    qsub -hard -l scratchfree=500,h_cpu=24:00:00,h_vmem=1.8G,eliza17io=1 -o log/job_$FILE.log -e log/job_$FILE.err ./run$FILE.csh

	    mv run$FILE.csh script/

	  done
        done
     done
  done
  echo All files produced
else
  echo "Wrong number of parameters"
fi
