function checkVariables {
	set -e

if [ "$PROJECTFOLDER" == "" ];then
	export PROJECTFOLDER="."
fi

if [ "$MODULE" == "" ] || [ "$MODULE" == "all" ] || [ "$MODULE" == "ALL" ]; then
	MODULE="statusb qc humanDecont statusa assignTaxonomy report"
fi

}