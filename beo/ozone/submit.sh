#!/bin/bash

rm ozone_*.R
rm *.out

cat sample.R | sed "s/WT/0/g" | sed "s/WL/0/g" | sed "s/WR/10/g" > ozone_0_0_1.R
CMD="bwsubmit_bash r ozone_0_0_1.R"
OUT=`$CMD`

for j in {4,8,16,32}; do
	for i in {4..-4}; do
		for t in {1..2}; do
			RFILE="ozone_${t}_${i}_${j}.R"
			cat sample.R | sed "s/WT/${t}/g" | sed "s/WL/${i}/g" | sed "s/WR/${j}/g" > $RFILE
			CMD="bwsubmit_bash r ${RFILE}"
			echo $CMD
			OUT=`$CMD`
		done
	done
done
