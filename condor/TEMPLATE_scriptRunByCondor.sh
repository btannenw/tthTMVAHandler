#!/bin/bash 

outputdir=`pwd`

cmssw_wd=__MYCMSSWDIR_
# Directory where the script calls "cmsenv".

jobdir=__JOBDIR___
# = Directory in which the command exists.


source /cvmfs/cms.cern.ch/cmsset_default.sh

cd ${cmssw_wd}
echo "## satoshi ## pwd = " `pwd` 

eval `scramv1 runtime -sh` # This is instead of cmsenv
echo "## satoshi ## CMSSW env done"

echo "## satoshi ## Moved to job directory."
cd ${jobdir}
echo "## satoshi ## pwd = " `pwd` 
echo "## satoshi ## ls = " `ls` 

# Need to be in the job directory in order to access some files with relative-path coded in the program.


#  "" is needed for nominal.
for systematics in "" SYSTEMATICSOPTION
do
${cmssw_wd}/${jobdir}/PROGRAMTORUN  INPUTOPTIONS -d ${outputdir} -o OUTPUTFILENAME SPECIALOPTIONS $systematics
done

cd ${outputdir}/
for f in *.root
do
xrdcp  ${f}  root://cmseos.fnal.gov//store/user/lpctthrun2/ttHTMVA/NICKNAMEOFJOBS/${f}
done

rm *.root 

echo "##########################"
echo "script ends"

#- - - If you want to submit CMSRUN job, you can do like this : 
#
#  cat ${pythonscript} | sed "s|OUTPUTDIR|${outputdir}|g" > ${outputdir}/job.py
#  cmsRun ${outputdir}/job.py

