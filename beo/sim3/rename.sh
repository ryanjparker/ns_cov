#!/bin/bash

for ofile in `ls sim2_*.R`
do
	nfile=`echo $ofile | sed 's/sim2/sim3/g'`
	nc=`cat "${ofile}" | sed 's/which_exp <- 2/which_exp <- 3/g'`
	echo "$ofile to $nfile"
	echo "$nc" >$nfile
done
