### v3 analysis (STAR picoDst only):
> to extract flow signal with systematic errors, please use PiKFlow, ProtonFlow and V0Flow

#### Copy Code:
> git clone /global/project/projectdirs/starprod/rnc/xusun/gitrepo/Analysis_v3.git ./Analysis_v3 (Recommended, then you can get updates by simply doint git pull origin master:master)
> Or copy it directly from: /global/homes/x/xusun/STAR/Analysis_v3

#### TriFlow.C (use 39 GeV as an example)
> Fill Histogram for v3 analysis: 
> 
> root4star -b -q 'TriFlow.C("39GeV_111.list",111,1,0,0)'
> - first argument is inputlist
> - second argument is jobCounter
> - third argument is analysis mode: 0 for re-center correction, 1 for shift correction, 2 for charged hadron, 3 for pion and kaon, 4 for proton, 5 for phi meson, 6 for Lambda, 7 for antiLambda
> - fourth argument is energy: 0 for 200 GeV, 1 for 39 GeV, 2 for 27 GeV (on progress)
> - fifth argument is flag_ME: 0 for Same Event, 1 for Mixed Event
> 
> The scripts to submit jobs are in /macro/pdsf/ (the script to submit on carver is on progress) 
> . ./39_TriFlow.sh (you can change submit jobs in script)
> NOTICE: please CHANGE the working directory in run.csh
> 
> the main code is in /StRoot/StTriFlowMaker/
> 
> - StTriFlowConstants.h/StTriFlowConstants.cxx is the class for all the constants used in the analysis
> - StTriFlowCut.h/StTriFlowCut.cxx is the class for event selection and track selection and PID cut
> - StTriFlowCorrection.h/StTriFlowCorrection.cxx is the class for all correction include re-center correction, shift correction and event plane reconstruction
> - StTriFlowProManger.h/StTriFlowProManger.cxx is the class for all the TProfile which are used in re-center correction, shift correction and charged hadron v3 calculation
> - StTriFlowHistoManger.h/StTriFlowHistoManger.cxx is the class for all the THistogram which are used in pion, kaon proton analysis
> - StTriFlowV0.h/StTriFlowV0.cxx is the class for all v0 analysis (only phi meson right now)
> - StTriFlowMEKey.h is the class for key used in v0 analysis
> - StTriFlowMaker.h/StTriFlowMaker.cxx is the main maker to do v3 analysis. There are three main function in this class:
>   - Init: initialize the class used in different analysis mode
>   - Make: main analysis loop for all tracks in one event (event loop is in TriFlow.C)
>   - Finish: save everything

#### Resolution calculation: 
> Calculate resolution after Shift Correction before Charged Hadron flow calculation
> 
> root4star -b -q Resolution.C\(231,232,0\)
> first argument is start job
> second argument is stop job
> thrid argument is energy: 0 for 200 GeV, 1 for 39 GeV, 2 for 27 GeV(on progress)
> 
> to submit job into batch:
> . ./39_Resolution.sh 39_Resolution_list (you can change submit jobs in 39_Resolution_list)
> 
> the main code is in /StRoot/StTriFlowResolution/
> 
> StTriFlowConstants.h/StTriFlowConstants.cxx is the class for the constants used in analysis
> 
> StTriFlowHistManger.h/StTriFlowHistManger.cxx is the class for the THistogram for event plane distribution
> 
> StTriFlowResolution.h/StTriFlowResolution.cxx is the main class to do resolution calculation. It is similar to StTriFlowCorrection.h/StTriFlowCorrection.cxx in /StRoot/StTriFlowMaker/ but have one more function to calculation event plane resoltuion (for some reason I wrote a separate code to calculate event plane resolution, I will add another mode in StTriFlowMaker to calculate event plane resoltuion some day)

------------------------------------------------------------------------

#### V0 recontruction: (phi, Lambda, anti-Lambda)

> root4star -b -q PhiFlow.C\(1,0,0,1000,1024,1\)
> first argument is Energy: 0 for 200 GeV, 1 for 39 GeV
> second argument is flag for Mixed Event: 0 for same event, 1 for mixed event
> third argument is inputlist: depends on how many output TTree you have
> fourth argument is start event
> fifth argument is stop event
> sixth argument is Mode: 0 for phi, 1 for Lambda, 2 for anti-Lambda
> 
> to submit job into batch:
> . ./39_PhiFlow.sh
> 
> the main code is in /StRoot/StStrangenessAna/
> 
> StStrangenessCons.h/StStrangenessCons.cxx is the class for constants used in analysis
> 
> StStrangenessCorr.h/StStrangenessCorr.cxx is the class for corrections used in analysis
> 
> StStrangenessCut.h/StStrangenessCut.cxx is the class for cuts used in analysis
> 
> StStrangenessHistoManger.h/StStrangenessHistoManger.h is the class for Histogram used in analysis
> 
> StStrangenessAna.h/StStrangenessAna.cxx is the main analysis class to recontruct v0 particles
> 
> /StRoot/StStrangenessAna/macro/PhiFlow.C is the macro for phi flow calculation
> /StRoot/StStrangenessAna/macro/LambdaFlow.C is the macro for Lambda flow calculation
> /StRoot/StStrangenessAna/macro/antiLambdaFlow.C is the macro for antiLambda flow calculation
> NOTICE: please CHANGE the input files to your own v0 histogram and create some output directory

#### Output Directory sturcture

> Please check: /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu39GeV
> 
> You may need to generate your own List to submit the jobs and the examples are in /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu39GeV/List.
> NOTICE: all the corrections are run by run basis, therefore, it is good to submit jobs run by run when calculate the correction factors. But when you try to calculate v0 like Lambda and anti-Lambda, run by run calculation is not possible anymore as the calculate time will exceed 24 hours and pdsf will kill your jobs. Then you need to split your jobs. You can find all the examples in the directory above.

If you have any problem, please contact: xsun@hit.edu.cn | sunxuhit@gmail.com
