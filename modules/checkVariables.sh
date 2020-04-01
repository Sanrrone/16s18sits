function checkVariables {
	set -e

if [ "$PROJECTFOLDER" == "" ];then
	export PROJECTFOLDER="."
fi

if [ "$TOLERANCE" == "" ];then
    TOLERANCE=0.1
fi

if [ "$(echo $TOLERANCE | awk '{if($1>1 || $1<0)print "fail"}')" == "fail"  ];then
    echo "TOLERANCE must be between 0 and 1"
    exit
fi

if [ "$MODULE" == "" ] || [ "$MODULE" == "all" ] || [ "$MODULE" == "ALL" ]; then
	#MODULE="statusb qc humanDecont statusa assignTaxonomy report"
	MODULE="statusb qc humanDecont statusa assignTaxonomy"
fi

if [ "$PAIRED" == "" ]; then
	PAIRED="TRUE"
fi

}

function printHelp {
	echo "usage: bash 16s18sits.sh -p [project folder] -f [fastq foder] -fp [fastq pattern]"
    echo "example: bash 16s18sits.sh -p mynewproject -f 0-raw -fp .fastq.gz"
    echo -e "\noptions available:\n"
    echo "-p Project folder, if no exist, the Pipeline will assume is the actual folder"
    echo "-f fastq folder, mandatory parameter"
    echo "-fp fastq pattern, it could be '.fastq' or '.fastq.gz'"
    exit
}
