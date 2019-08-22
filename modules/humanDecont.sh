function humanDecont {
	set -e
	echo "ESCLAVO: humanDecont begin"

	centrifugeBIN=$(which centrifuge)
	if [ "$centrifugeBIN" == "" ];then
		echo "ESCLAVO: no centrifuge installed, aborting this step"
		exit
	fi

	if [ ! -d 3-decont ];then
		mkdir 3-decont
	fi

	cd 3-decont

	if [ ! -f $ESCLAVOHOME/DB/humanDB.[123].cf];then
		echo "No humanDB found, donwloading..."
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1QZUKHj040nFSU69ojBTU2uOhb2MYNIOa' -O humanDB.1.cf
		mv humanDB.1.cf $ESCLAVOHOME/DB/.
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1fA-0CyIULPSWKryru0q14bqyWr7z1MnE' -O humanDB.2.cf
		mv humanDB.2.cf $ESCLAVOHOME/DB/.
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1OksPc2EOJHYtOXPtFL2QE9AMLsT4dCTt' -O humanDB.3.cf
		mv humanDB.3.cf $ESCLAVOHOME/DB/.
		echo "Done"
	fi


	centrifuge -x $ESCLAVOHOME/DB/humanDB -q -m1 readsf1 -m2 readsf2 --un decontaminated -p $(nproc) --out-fmt tab


	echo "ESCLAVO: humanDecont end"
}