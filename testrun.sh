#!/bin/sh


OUTPUTDIR=output
rm   -rf ${OUTPUTDIR}
mkdir -p ${OUTPUTDIR}


#for channel in ttbar tth DataEl DataMu wjets qcd_ht100 qcd_ht200 qcd_ht300 qcd_ht500 qcd_ht700 qcd_ht1000 qcd_ht1500 qcd_ht2000

for channel in ttbar tth DataElB DataElC DataElD DataMuB DataMuC DataMuD
do

INTUPLES=../ttHAnalyzer/condor/


EX_OPTION=" "

if [ `echo ${channel} | grep  Data ` ]
then
EX_OPTION=" -D "
#          \__ This tells "input is data file"
fi

inputoptions=""
for rootfile in ` ls ${INTUPLES}/output__${channel}_subjob*.root | head -n 10000 `
do
inputoptions="${inputoptions} -i $rootfile "
done

outputoption=" -o ${OUTPUTDIR}/${channel}"

./main $inputoptions  $outputoption  ${EX_OPTION} &


done

wait

echo "done"
