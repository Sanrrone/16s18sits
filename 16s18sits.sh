#!/bin/bash
set -e
shopt -s direxpand
shopt -s expand_aliases
#if [ -f ~/.bash_profile ]; then source ~/.bash_profile; fi
#if [ -f ~/.bashrc ]; then source ~/.bashrc; fi
#if [ -f ~/.bash_alias ]; then source ~/.bash_alias; fi

export ESCLAVOHOME=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
source $ESCLAVOHOME/modules/checkVariables.sh


#usage: bash 16s18sits.sh --force -p /home/sandro/Programas/ESCLAVO/projects -f /home/sandro/Programas/ESCLAVO/0-raw -pt .fastq.gz
#software requirements
# FASTQC (better download and create alias), MULTIQC
# sudo apt install pandoc
#R packages: DADA2, rmarkdown, rmdformats, DT, dplyr, ggplot2, plotly, stringr, Biostrings
if [ "$1" == "" ];then
    printHelp
fi
actualfolder=$(pwd)
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -p|--projectfolder)
    cd $2
    PROJECTFOLDER=$(pwd)
    cd $actualfolder
    shift # past argument
    shift # past value
    ;;
    -f|--fastqfolder)
    cd $2
    FASTQFOLDER=$(pwd)
    cd $actualfolder
    shift # past argument
    shift # past value
    ;;
    -fp|--fastqpattern)
    PATTERN="$2"
    shift # past argument
    shift # past value
    ;;
    -to|--tolerance)
    TOLERANCE="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--module)
    MODULE="$2"
    shift # past argument
    shift # past value
    ;;
    -ler|--learnErrorReads)
    READS4LEARN="$2"
    shift # past argument
    shift # past value
    ;;   
    -F|-force|--force)
    FORCE=true
    shift # past argument
    ;;
    --debug)
    DEBUG=true
    set -ex
    shift # past argument
    ;;
    -h|--help)
    printHelp
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -n $1 ]]; then
    echo "the following argument lack of parameter: $1"
    exit
fi

checkVariables


cd $PROJECTFOLDER
export PNAME=$(pwd | awk -F"/" '{print $NF"_eConf.tsv"}')
export PCONF=$(pwd | awk -v pname=$PNAME '{print $0"/"pname}')

if [ ! -f "$PNAME" ] || [ $FORCE ]; then
	echo "	value" > $PNAME
    echo "pname	$(pwd | awk -F"/" '{print $NF}')" >> $PNAME
	echo "pfolder	$(pwd)" >> $PNAME
	echo "ffolder	$FASTQFOLDER" >> $PNAME
	echo "fqpattern	$PATTERN" >> $PNAME
	echo "analysis	16s18sits" >> $PNAME
	echo "aversion	1.0" >> $PNAME
	echo "created	$(date --iso-8601)" >> $PNAME
	echo "status	open" >> $PNAME
	echo "pPercent	0" >> $PNAME
	echo "lastStep	statusb" >> $PNAME
fi

for mod in $MODULE
do
    case $mod in
	   "statusb")
		  source $ESCLAVOHOME/modules/statusb.sh
		  statusb $FORCE
	   ;;
	   "humanDecont")
		  source $ESCLAVOHOME/modules/humanDecont.sh
		  humanDecont 1-qc
	   ;;
	   "qc")
		  source $ESCLAVOHOME/modules/qc.sh
		  qc
	   ;;
	   "statusa")
		  source $ESCLAVOHOME/modules/statusa.sh
		  statusa $FORCE
	   ;;
	   "assignTaxonomy")
		  source $ESCLAVOHOME/modules/assignTaxonomy.sh
		  assignTaxonomy 3-taxInsight
	   ;;
	   "report")
		  source $ESCLAVOHOME/modules/report.sh
		  makeReport 0-raw 2-decont 3-taxInsight
	   ;;
	   *)
		  echo "Module $mod not recognized"
	   ;;
    esac
    cd $PROJECTFOLDER
done

echo "ESCLAVO: Pipeline done :)"

    
#	    ___	_,.--.,_			 ___	_,.--.,_			  ___	_,.--.,_			___	_,.--.,_    
#	 .-~   ~--"~-.   ._ "-.		.-~   ~--"~-.   ._ "-.		.-~   ~--"~-.   ._ "-.		.-~   ~--"~-.   ._ "-.    
#	/	 ./_    Y    "-. \	   /	 ./_    Y    "-. \	   /	 ./_    Y    "-. \	   /	 ./_    Y    "-. \    
#    Y	  :~	!	    Y	   Y	  :~	!	    Y	 Y	  :~	!	    Y	 Y	  :~	!	    Y    
#    lq p    |	/	    .|	   lq p    |	/	    .|	 lq p    |	/	    .|	 lq p    |	/	    .|    
# _   \. .-, l    /		|j   _   \. .-, l    /		|j   _   \. .-, l    /		|j   _   \. .-, l    /		|j    
#()\___) |/   \_/";		!   ()\___) |/   \_/";		!   ()\___) |/   \_/";		!   ()\___) |/   \_/";		!    
# \._____.-~\  .  ~\.	 ./	\._____.-~\  .  ~\.	 ./	\._____.-~\  .  ~\.	 ./	\._____.-~\  .  ~\.	 ./    
#		  Y_ Y_. "vr"~  T			  Y_ Y_. "vr"~  T			  Y_ Y_. "vr"~  T			  Y_ Y_. "vr"~  T    
#		  (  (    |L    j			  (  (    |L    j			  (  (    |L    j			  (  (    |L    j    
#		  [nn[nn..][nn..]			  [nn[nn..][nn..]			  [nn[nn..][nn..]			  [nn[nn..][nn..]    
#	   ~~~~~~~~~~~~~~~~~~~~~~~	    ~~~~~~~~~~~~~~~~~~~~~~~	    ~~~~~~~~~~~~~~~~~~~~~~~	    ~~~~~~~~~~~~~~~~~~~~~~~   
