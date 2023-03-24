#!/bin/bash

for i in tPh0*;                                                                    
do 
	tPh=$(echo "$i" | cut -d "_" -f1)
	tPh=${tPh:3}
	tso=$(echo "$i" | cut -d "_" -f2)
	tso=${tso:3}
	tsoPh=$(echo "$i" | cut -d "_" -f3)
	tsoPh=${tsoPh:5}
	P=$(echo "$i" | cut -d "_" -f4)
	P=${P:1}
	if [ "$tPh" = 0 ];
	then 
		tPh=0.0;
	fi
	if [ "$tso" = 0 ];
	then 
		tso=0.0;
	fi 
	if [ "$tsoPh" = 0 ];
	then	
		tsoPh=0.0;
	fi
		
	cp -r plot.jl $i;
	cd $i;
	sed -i "s/tPh = 0.0/tPh = ${tPh}/g" plot.jl
	sed -i "s/SOC = 0.0/SOC = ${tso}/g" plot.jl
	sed -i "s/SOCPh = 0.0/SOCPh = ${tsoPh}/g" plot.jl
	sed -i "s/P = 0.0/P = ${P}/g" plot.jl
	cd ../
done
