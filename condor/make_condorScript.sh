#!/bin/sh

# Prepare EOS.
source /etc/profile.d/eos_aliases.sh      # In case of LPC


rm CheckJobResult.sh
touch CheckJobResult.sh


CMSSW_DIR=${CMSSW_BASE}/src/


INPUT_DATA_DIR=/store/user/lpctthrun2/ttHAnalyzer/20170724_ObjSeleRevisedAndUpdated


nickname="2017_08_03_FakeEstAdded"

JOBDIR="ttHTMVAHandler"
# Job DIR : where the executable exists.


FakeEstimationON="YES"


# if false : treat as one ttbar-inclusive
# if  true : prepare subjobs badsed on additionalJetId.
TTBarAdditionalJetID="true"




function registerTTbarMC(){

if [ _$TTBarAdditionalJetID = "false" ] 
then

ijob=${ijob}+1 ;
channel[${ijob}]=${1}

else

prefix[1]=ttPlusB    
prefix[2]=ttPlus2B   
prefix[3]=ttPlusBBbar
prefix[4]=ttPlusCCbar
prefix[5]=ttOther    

for n in 1 2 3 4 5 
do
ijob=${ijob}+1 ;
channel[${ijob}]=${1}_${prefix[$n]}
done

fi

}

function getSpecialOption(){

channel_name=${1}

SPOP="";

# use only even enve 
if [[ $channel_name == *"tttosemilep"*     ]] ; then SPOP=" -O "  ; fi
if [[ $channel_name == *"ttto2l2nu"*       ]] ; then SPOP=" -O "  ; fi
if [[ $channel_name == "tth"*             ]] ; then SPOP=" -O "  ; fi

if [[ $channel_name == *"ttPlusB"     ]] ; then SPOP=${SPOP}" -t 1 " ; fi 
if [[ $channel_name == *"ttPlus2B"    ]] ; then SPOP=${SPOP}" -t 2 " ; fi 
if [[ $channel_name == *"ttPlusBBbar" ]] ; then SPOP=${SPOP}" -t 3 " ; fi 
if [[ $channel_name == *"ttPlusCCbar" ]] ; then SPOP=${SPOP}" -t 4 " ; fi 
if [[ $channel_name == *"ttOther"     ]] ; then SPOP=${SPOP}" -t 5 " ; fi 

echo ${SPOP}

}

function getSampleName(){

channel_name=${1}
if [[ $channel_name == *"ttPlusB"     ]] ; then echo ${channel_name%_ttPlusB}     ; return ; fi 
if [[ $channel_name == *"ttPlus2B"    ]] ; then echo ${channel_name%_ttPlus2B}    ; return ; fi 
if [[ $channel_name == *"ttPlusBBbar" ]] ; then echo ${channel_name%_ttPlusBBbar} ; return ; fi 
if [[ $channel_name == *"ttPlusCCbar" ]] ; then echo ${channel_name%_ttPlusCCbar} ; return ; fi 
if [[ $channel_name == *"ttOther"     ]] ; then echo ${channel_name%_ttOther}     ; return ; fi 

echo ${channel_name}

}








declare -i ijob
ijob=0

ijob=${ijob}+1 ;
channel[${ijob}]="DataElB"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElC"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElD"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElE"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElF1"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElF2"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElG"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElHv2"

ijob=${ijob}+1 ;
channel[${ijob}]="DataElHv3"



ijob=${ijob}+1 ;
channel[${ijob}]="DataMuB"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuC"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuD"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuE"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuF1"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuF2"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuG"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuHv2"

ijob=${ijob}+1 ;
channel[${ijob}]="DataMuHv3"






# ijob=${ijob}+1 ;
# channel[${ijob}]="ttbar"

# ijob=${ijob}+1 ;
# channel[${ijob}]="ttto2l2nu"
# ijob=${ijob}+1 ;
# channel[${ijob}]="tttosemilep"
# 

registerTTbarMC ttto2l2nu
registerTTbarMC tttosemilep


ijob=${ijob}+1 ;
channel[${ijob}]="wjetsincl"

ijob=${ijob}+1 ;
channel[${ijob}]="zjetsincl"


ijob=${ijob}+1 ;
channel[${ijob}]="ww"

ijob=${ijob}+1 ;
channel[${ijob}]="wz"

ijob=${ijob}+1 ;
channel[${ijob}]="zz"



ijob=${ijob}+1 ;
channel[${ijob}]="tth"


ijob=${ijob}+1 ;
channel[${ijob}]="ZjetLowMass"

ijob=${ijob}+1 ;
channel[${ijob}]="schan_both"

ijob=${ijob}+1 ;
channel[${ijob}]="tchan_top"

ijob=${ijob}+1 ;
channel[${ijob}]="tbarW"


ijob=${ijob}+1 ;
channel[${ijob}]="tchan_tbar"

ijob=${ijob}+1 ;
channel[${ijob}]="tW"



ijob=${ijob}+1
channel[${ijob}]="TTWJetsToLNu_1"


ijob=${ijob}+1
channel[${ijob}]="TTWJetsToLNu_2"

ijob=${ijob}+1
channel[${ijob}]="TTWJetsToQQ"


ijob=${ijob}+1
channel[${ijob}]="TTZToLLNuNu_1"

ijob=${ijob}+1
channel[${ijob}]="TTZToLLNuNu_2"

ijob=${ijob}+1
channel[${ijob}]="TTZToLLNuNu_3"

ijob=${ijob}+1
channel[${ijob}]="TTZToQQ"



# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT1000to1500"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT100to200"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT1500to2000"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT2000toInf"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT200to300"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT300to500"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT500to700"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT50to100"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="QCDHT700to1000"
# 
# 


# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT100to200"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT1200to2500"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT200to400"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT2500toInf"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT400to600"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT600to800"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT70to100"
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]="DYM50HT800to1200"
# 



# 1 - 18 : data el and mu.
# 1 : ttbar 
# 2 : ttbar LJ, DL
# 1 : w
# 1 : z
# 3 : di bosons
# 1 : ttH
# 1 : low Z mass
# 5 single top 
# 9 : QCD
#
# sum 42 

if [ $FakeEstimationON = "YES" ]
then

for i in `seq 1 ${ijob}`
do

echo checking $i 

ijob=${ijob}+1 ;
filepath[${ijob}]=${filepath[${i}]}

if [ `echo ${channel[$i]} | grep  Data ` ]
then
channel[${ijob}]=Fake${channel[${i}]#Data}
# Point 1 : Category name is FakeXXX where XXX comes from the original DataEl/Mu[Perido].
else 
channel[${ijob}]=Fake${channel[${i}]}
fi

done
# event loop 

fi
# end if FakeEstimationMode 





# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarmass1735
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarmass1715
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbartuneup
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbartunedown
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarhdampup
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarhdampdown
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarwidth0p2
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarwidth0p5
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarwithd4
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarwidth8
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarcolorflip
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=ttbarerdon
# 
# 
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted1ttto2l2nu
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted2ttto2l2nu
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted3ttto2l2nu
# 
# 
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted1tttosemilep
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted2tttosemilep
# 
# ijob=${ijob}+1 ;
# channel[${ijob}]=TopGenReweighted3tttosemilep
# 




registerTTbarMC ttbarmass1735
registerTTbarMC ttbarmass1715
registerTTbarMC ttbartuneup
registerTTbarMC ttbartunedown
registerTTbarMC ttbarhdampup
registerTTbarMC ttbarhdampdown
registerTTbarMC ttbarwidth0p2
registerTTbarMC ttbarwidth0p5
registerTTbarMC ttbarwithd4
registerTTbarMC ttbarwidth8
registerTTbarMC ttbarcolorflip
registerTTbarMC ttbarerdon
registerTTbarMC TopGenReweighted1ttto2l2nu
registerTTbarMC TopGenReweighted1tttosemilep


registerTTbarMC ttbarFSRDown
registerTTbarMC ttbarFSRUp
registerTTbarMC ttbarISRDown
registerTTbarMC ttbarISRUp
registerTTbarMC ttbarUEdown
registerTTbarMC ttbarUEup




#  66 = all but QCD.
# ttbars sytematics : total 12.
#                   ( 67 - 78 )
# ttbar LJ/DL x 3 TopGenWeight = 6
#                   ( 79 - 84

for i in `seq 1 ${ijob}`
do

#------------------------------------------
#--- Setting which you may need to change -----
#------------------------------------------

# Input file directory. 
echo "satoshi test "`getSampleName ${channel[$i]}`
indir=${INPUT_DATA_DIR}/`getSampleName ${channel[$i]}`

# Prefix of output root file. 
prefixoutputrootfile=output__${channel[$i]}

executable=./main

# Number of maximum input files
NMAXFILES=9999

# Number of files per one condor job
NFILESFORONEJOB=3

#--------------------------
#----- Setting end --------
#--------------------------



SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)


source /cvmfs/cms.cern.ch/cmsset_default.sh
eos root://cmseos.fnal.gov mkdir /store/user/lpctthrun2/ttHTMVA/
eos root://cmseos.fnal.gov mkdir /store/user/lpctthrun2/ttHTMVA/${nickname}
eos root://cmseos.fnal.gov mkdir /store/user/lpctthrun2/ttHTMVA/${nickname}/${channel[$i]}/




# --------
# Preparation of input files is done
#


inputs=""
declare -i iFiles
declare -i iJOB
declare -i iFilesForThisJob
iFiles=0
iJOB=1
iFilesForThisJob=0
inputs[$iJOB]="";

#  for subdir in ` eosls ${indir} `
#  do
#  
#  for subsubdir in ` eosls ${indir}/${subdir} `
#  do
#  
#  for thefile in ` eosls ${indir}/${subdir}/${subsubdir} | grep ".root" `
#  do

echo eosls ${indir}
eosls ${indir} | wc 

for thefile in ` eosls "${indir}" | grep ".root" `
do

echo ${thefile}

if [ $iFiles -lt ${NMAXFILES} ]
then

inputs[$iJOB]="${inputs[$iJOB]} -i root://cmseos.fnal.gov/${indir}/${thefile}"
iFilesForThisJob=${iFilesForThisJob}+1

if [ ${iFilesForThisJob} -ge $NFILESFORONEJOB ]
then
iJOB=${iJOB}+1
inputs[$iJOB]="";
iFilesForThisJob=0
fi

fi

iFiles=${iFiles}+1

done
# ls-Loop end.


#  done
#  done
#  done


#Adjust last set (It can be empty. If this is the case, ignore it)
if [ $iFilesForThisJob -eq 0 ]
then
iJOB=${iJOB}-1
fi

echo "###########"
echo ${filepath[${i}]}
echo "Number of found files : "${iFiles}
echo "Total NJob : "${iJOB}
echo "###########"

cat >> CheckJobResult.sh <<EOF
echo "# - - - - - - "
echo "Number of file  for ${channel[$i]}"
echo " - should be : ${iJOB}"
for syst in "" \
JERUP \
JERDOWN \
JESUP \
JESDOWN \
btag_LFUp \
btag_LFDown \
btag_HFUp \
btag_HFDown \
btag_HFStats1Up \
btag_HFStats1Down \
btag_HFStats2Up \
btag_HFStats2Down \
btag_LFStats1Up \
btag_LFStats1Down \
btag_LFStats2Up \
btag_LFStats2Down \
btag_CharmStat1up \
btag_CharmStat1down \
btag_CharmStat2up \
btag_CharmStat2down 
do 
nameforgrep=output\${syst}_${prefixoutputrootfile}_subjob
numberoffiles=\` ls *root | grep \${nameforgrep} | wc -l \`
if [ \$numberoffiles -eq ${iJOB} ]
then
echo -n "[GOOD] "
else
echo -n "[ BAD] "
fi
echo "  # of \${syst} : " \$numberoffiles
done

EOF


declare -i jJob
jJob=1
while [ $jJob -le ${iJOB} ]
do

outputoption="${prefixoutputrootfile}_subjob${jJob}"

scriptname=__script_${channel[$i]}_subjob${jJob}.sh
condorscript=__condorJobScript_${channel[$i]}_subjob${jJob}.py


EX_OPTION=""
SYST_OPTION=""
if [ `echo ${channel[$i]} | grep  DataMu ` ]
then
    EX_OPTION=" -M "
elif [ `echo ${channel[$i]} | grep  DataEl ` ]
then
    EX_OPTION=" -E "
elif [ `echo ${channel[$i]} | grep  FakeMu ` ]
then
    EX_OPTION=" -M -F "
    # Add FakeEstimationRelated uncertainty.
    for syst in `seq 101 104 `
    do
	SYST_OPTION="${SYST_OPTION} \" -S ${syst}\" "
    done
elif [ `echo ${channel[$i]} | grep  FakeEl ` ]
then
    EX_OPTION=" -E -F "
    # Add FakeEstimationRelated uncertainty.
    for syst in `seq 101 104 `
    do
	SYST_OPTION="${SYST_OPTION} \" -S ${syst}\" "
    done
else 

# if this is MC : 

if [[ ${channel[${i}]} = *"ttbarmass1735"*   || \
      ${channel[${i}]} = *"ttbarmass1715"*   || \
      ${channel[${i}]} = *"ttbartuneup"*     || \
      ${channel[${i}]} = *"ttbartunedown"*   || \
      ${channel[${i}]} = *"ttbarhdampup"*    || \
      ${channel[${i}]} = *"ttbarhdampdown"*  || \
      ${channel[${i}]} = *"ttbarwidth0p2"*   || \
      ${channel[${i}]} = *"ttbarwidth0p5"*   || \
      ${channel[${i}]} = *"ttbarwithd4"*     || \
      ${channel[${i}]} = *"ttbarwidth8"*     || \
      ${channel[${i}]} = *"ttbarcolorflip"*  || \
      ${channel[${i}]} = *"ttbarerdon"*      || \
      ${channel[${i}]} = *"ttbarFSR"*      || \
      ${channel[${i}]} = *"ttbarISR"*      || \
      ${channel[${i}]} = *"ttbarUE"*      || \
      ${channel[${i}]} = *"TopGenReweighted"*  ]]
then 
 echo "For ${channel[${i}]}, only nominal is process -- no systematic trees."
else 

for syst in `seq 5 78 `
do
SYST_OPTION="${SYST_OPTION} \" -s ${syst}\" "
done
SYST_OPTION="${SYST_OPTION} \" -p\" "
SYST_OPTION="${SYST_OPTION} \" -P\" "

#  Scale factor systematics `
for syst in `seq 1 12 `
do
SYST_OPTION="${SYST_OPTION} \" -S ${syst}\" "
done

fi


if [ `echo ${channel[$i]} | grep  Fake ` ]
then
   # in addition, this is also Fake estimation.
    EX_OPTION="${EX_OPTION} -F "

    # Add FakeEstimationRelated uncertainty.
    for syst in `seq 101 104 `
    do
	SYST_OPTION="${SYST_OPTION} \" -S ${syst}\" "
    done

fi

fi


# 
EX_OPTION="${EX_OPTION} "`getSpecialOption ${channel[$i]} `


cat TEMPLATE_scriptRunByCondor.sh \
  | sed "s|PROGRAMTORUN|$executable|g" \
  | sed "s|INPUTOPTIONS|${inputs[$jJob]}|g" \
  | sed "s|OUTPUTFILENAME|${outputoption}|g" \
  | sed "s|__MYCMSSWDIR_|${CMSSW_DIR}| g" \
  | sed "s|__JOBDIR___|${JOBDIR}|" \
  | sed "s|NICKNAMEOFJOBS|${nickname}/${channel[$i]}|" \
  | sed "s|SPECIALOPTIONS|${EX_OPTION}|" \
  | sed "s|SYSTEMATICSOPTION|${SYST_OPTION}|" \
  > ${scriptname}

cat template__condorSetting \
 | sed "s|__SCRIPT_TO_RUN__|${scriptname}|g" \
 > ${condorscript}

jJob=${jJob}+1

done


done

cat <<EOF
#### execute this ###
for file in __condorJobScript*py
do
condor_submit \$file
done
###############
EOF
