function checkVariables {
	set -e

if [ "$PROJECTFOLDER" == "" ];then
	export PROJECTFOLDER="."
fi

if [ "$MODULE" == "" ] || [ "$MODULE" == "all" ] || [ "$MODULE" == "ALL" ]; then
	MODULE="statusb qc humanDecont statusa assignTaxonomy report"
fi

}

function printHelp {
	echo "usage: bash 16s18sits.sh -p [project folder] -f [fastq foder] -pt [fastq pattern]"
    echo "example: bash 16s18sits.sh -p mynewproject -f 0-raw -pt .fastq.gz"
    echo -e "\noptions available:\n"
    echo "-p Project folder, if no exist, the Pipeline will assume is the actual folder"
    echo "-f fastq folder, mandatory parameter"
    echo "-pt fastq pattern, it could be '.fastq' or '.fastq.gz'"
    exit
}