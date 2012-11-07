#!/bin/bash

##########################################################################################
### Author........Nil Valls <nvallsve@nd.edu>                                          ###
### Created.......2012/10/28                                                           ###
### Modified......2012/11/01                                                           ###
### Description...This script will generate the multicrab output from a text file      ###
###    containing the pasted text from the ttH twiki tables containing the dataset     ###
###    information. The sample number must be in the first column, and the primary     ###
###    dataset in the second. The script will check DBS and figure out what dbs        ###
###    instance to find such PDS.                                                      ###
##########################################################################################

### User-defined parameters (change these as needed)
project="BEAN_53xOn52x_V01_CV01"
workingDir="/afs/crc.nd.edu/user/n/nvallsve/data/Production/"
pset="makeBEAN_cfg.py"
storageElement="T3_US_NotreDame"
numJobs=450
json="Collisions12/8TeV/Prompt/Cert_190456-205618_8TeV_PromptReco_Collisions12_JSON.txt"


##########################################################################################
### After this point, any modifications risk messing up the script. Hack responsibly!  ###
##########################################################################################
RED="\033[0;31m"
GRAY="\033[0;30m"
BLUE="\033[0;34m"
ORANGE="\033[0;33m"
PURPLE="\033[0;35m"
GREEN="\033[0;32m"
WHITE="\033[1;37m"
NOCOLOR="\e[0m"

#function echoWar(){ echo -e "$PURPLE$1$NOCOLOR"; }
#function echoInf(){ echo -e "$BLUE$1$NOCOLOR"; }
#function echoSuc(){ echo -e "$GREEN$1$NOCOLOR"; }
#function echoWar(){ echo -e "$PURPLE[ WARNING ] $1$NOCOLOR" >&2; exit 1; }
function echoWar(){ echo -e "$PURPLE[ WARNING ] $1$NOCOLOR" >&2; }
function echoErr(){ echo -e "$RED[  ERROR  ] $1$NOCOLOR" >&2; rm -rf "$tempfile" "$output"; exit 1; }

input="$1"
output="$2"
if [ -z "$input" ]; then echoErr "must provide input file as the first argument."; fi
if [ -z "$output" ]; then echoErr "must provide output file as the second argument."; fi
touch "$output" > /dev/null
if [ $? -ne 0 ]; then echoErr "output file '$output' unwritable."; fi

## Check that we have a valid proxy
proxyExpiration=`voms-proxy-info | grep timeleft | sed -e 's/.*:\ //g' -e 's/://g'`
if [[ $proxyExpiration -eq 0 ]]; then echo -n "[ ATTENTION ] You don't seem to have a valid proxy. Would you like to get one (Y) or skip the DBS instance autoquery (N)? ";
	read getProxy;
	while [[ "$getProxy" != "Y" ]] && [[ "$getProxy" != "y" ]] && [[ "$getProxy" != "N" ]] && [[ "$getProxy" != "n" ]]; do
		echo -n "Invalid answer, please say 'Y' or 'N': "; read getProxy;
	done
	if [[ "$getProxy" == "Y" ]] || [[ "$getProxy" == "y" ]]; then
		skipDBS=0;
		proxy=`voms-proxy-init -voms cms`
		if [ $? -ne 0 ]; then echoErr "Terminating."; fi
		echo "":
	else
		skipDBS=1;
	fi
fi


mkdir -p "/tmp/$USER/"
tempfile="/tmp/$USER/.samples_from_twiki"

## Parse input file into temp file
#grep 'AOD\|USER' $input | sed 's/AOD\s.*/AOD/' | sed 's/AODSIM\s.*/AODSIM/' | sed 's/USER\s.*/USER/' > "$tempfile"
grep 'AOD\|USER' $input > "$tempfile"

## Common parts for multicrab
header="[MULTICRAB]
cfg                             = mcrab.cfg

[COMMON]
CRAB.jobtype                    = cmssw
CRAB.scheduler                  = condor

CMSSW.pset                      = ${pset}
CMSSW.output_file               = ttH_pat2bean_53x.root

USER.ui_working_dir             = $workingDir$(date +%Y%m%d)_${project}/
USER.return_data                = 0
USER.copy_data                  = 1
USER.check_user_remote_dir      = 0
USER.dontCheckSpaceLeft         = 1
USER.storage_element            = ${storageElement}
USER.publish_data               = 1
USER.dbs_url_for_publication    = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_01_writer/servlet/DBSServlet"


## This function echoes the first dbs instance where the input dataset is found (global, ph1, ph2)
function getDBSinstance(){
	if [[ $skipDBS -eq 1 ]]; then echo "<fill me>"; return; fi
	input_ds="$1";
	if [ -z "$input_ds" ]; then echoErr "getDBSinstance needs a dataset as argument"; fi

	urls=("http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet" \
	"https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_01_writer/servlet/DBSServlet" \
	"https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet")

	for url in ${urls[@]}; do
			dbsout=`python $DBSCMD_HOME/dbsCommandLine.py -c search --url="$url" --query="find dataset where dataset=$input_ds" | grep "$input_ds"`
			#e="$?"; if [ $e -ne 0 ]; then echoErr "Problem with dbsCommandLine.py '$e'"; fi
			if [[ "$dbsout" == "$input_ds" ]]; then echo "$url"; return 0; fi
	done
	echo "";
	echoErr "dataset '$input_ds' not found in any DBS instance. Please check.";
}

function getRecoType(){
	ds="$1"
	dssub=`echo "$ds" | grep Run20 | sed -e 's/.*Run20....//' -e 's/\/AOD.*//'`

	if [ -z "$dssub" ]; then echoErr "Trying to obtain reco type for non collision dataset: $1"; fi

	promtreco=`echo "$dssub" | grep -i "PromptReco"`
	if [ ! -z "$promptreco" ]; then echo "PR"; return 0; fi

	recover=`echo "$dssub" | grep -i "recover"`
	if [ ! -z "$recover" ]; then echo "RRr"; return 0; fi

	echo "RR"; return 0;

}


function getJSON(){
	ds="echo $1 | awk '{print $2}'"
	dssub=`echo "$ds" | grep Run20 | sed -e 's/.*Run20....//' -e 's/\/AOD.*//'`
	if [ -z "$dssub" ]; then echoErr "Trying to obtain reco type for non collision dataset: $1"; fi

	jsonPath="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions"
	era=`echo "$ds" | grep Run20 | sed -e 's/.*Run20\(..\).*/\1/'`
	   if [ "$era" == "11" ]; then jsonPath+="11/7TeV/";
	elif [ "$era" == "12" ]; then jsonPath+="12/8TeV/";
	else echoErr "Could not determine era while attempting to get JSON file."
	   fi

	recoType=`getRecoType "$ds"`
	if [ "$recoType" == "PR" ]; then jsonPath+="Prompt/";
	else jsonPath+="Reprocessing/"; fi

	filename=`echo "$line" | sed -e 's/.*Cert/Cert/' -e 's/\(\.txt\).*/\1/'`
	jsonInFilename=`echo "$filename" | grep -i json`
	if [ -z "$jsonInFilename" ]; then filename=""; fi

	if [ ! -f "$jsonPath$filename" ]; then
		echo "" >&2;
		echoWar "JSON file does not exist: $jsonPath$filename"
		echoWar "lumi_mask for this task will contain its directory only!"
	else
		jsonPath+="$filename"
	fi

	echo "$jsonPath"
}


## Function called for each task
function getBlock(){

	## Set up values
	DS="$1"
	NUM="$2"
	JSON="$3"
	DSU=`echo $DS | sed "s/^\///" | sed "s/\//_/g" | sed "s/_AODSIM//g" | sed "s/_AOD//g" | sed "s/_USER//g"`

	## Determine era
	if [[ "$DS" == *Run2011* ]] || [[ "$DS" == *7TeV* ]]; then
		era="2011";
	elif [[ "$DS" == *Run2012* ]] || [[ "$DS" == *8TeV* ]]; then
		era="2012";
	fi

	## Determine subera
	  if [[ "$DS" == *Run${era}A* ]]; then 
		subEra="A";
	elif [[ "$DS" == *Run${era}B* ]]; then 
		subEra="B";
	elif [[ "$DS" == *Run${era}C* ]]; then 
		subEra="C";
	fi

## Set up task
common="[$DSU]
CMSSW.datasetpath               = $DS
CMSSW.dbs_url                   = $(getDBSinstance $DS)
USER.publish_data_name          = ${DSU}_${project}
CMSSW.number_of_jobs            = $numJobs"

collisions="$common
CMSSW.total_number_of_lumis     = -1
#CMSSW.runselection              = 190456-196531
CMSSW.lumi_mask                 = $JSON
CMSSW.pycfg_params              = jobParams=${era}_${subEra}_data-$(getRecoType $DS)_$NUM"

background="$common
CMSSW.total_number_of_events    = -1
CMSSW.pycfg_params              = jobParams=${era}_X_MC-bg_$NUM"

if [[ "$DS" == *FastSim* ]]; then SIM="FastSim"; else SIM="FullSim"; fi
signal="$common
CMSSW.total_number_of_events    = -1
CMSSW.pycfg_params              = jobParams=${era}_X_MC-sig${SIM}_${NUM}"

	if [[ "$DS" == *Run20* ]]; then 
		echo "$collisions"
	elif [[ "$DS" == /TTH* ]]; then 
		echo "$signal"
	elif [[ "$DS" == */AODSIM ]]; then 
		echo "$background"
	else
		echoErr "could not figure out sample type for '$DS'";
	fi

	echo -ne "\n\n"
}


## Print the header
echo -ne "$header\n\n\n" > "$output"

## Print each task
echo "Found $(cat $tempfile | wc -l) datasets:"
while read line; do
	num=`echo $line | awk '{print $1}'`
	ds=`echo $line | awk '{print $2}'`
	echo -ne "${PURPLE}Preparing '${NOCOLOR}${ORANGE}$ds${PURPLE}'...${NOCOLOR}"
	json="$(getJSON $line)"
	getBlock "$ds" "$num" "$json" >> "$output"
	echo -e "${GREEN} done!${NOCOLOR}";
done < "$tempfile"

echo -e "\n${GREEN}\tAll done!\n${NOCOLOR}";


## Clean up
rm -rf "$tempfile"

#EOF
