#!/bin/bash

for ofile in `ls sim1_*.R`
do
	nfile=`echo $ofile | sed 's/sim1/sim2/g'`
	nc=`cat "${ofile}" | sed 's/which_exp <- 1/which_exp <- 2/g'`
	echo "$ofile to $nfile"
	echo "$nc" >$nfile
done
