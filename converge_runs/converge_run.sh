#!/bin/bash
#./2d_CSD/csd -Nx 50 -Ny 50 
rm timing.txt
# for i in "64" "128"
for i in "32" "64" "128"
do
	opts="-Nx $i -Ny $i -Nt 0.005"
	echo $opts
	./csd $opts
	mv data_csd.txt "data_csd_${i}.txt"
done
mv data_csd_* dt005/
mv timing.txt dt005/
