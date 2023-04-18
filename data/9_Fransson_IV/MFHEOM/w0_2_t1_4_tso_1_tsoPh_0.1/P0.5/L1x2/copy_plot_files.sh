#!/bin/bash

for i in {48,55,65}; 
do 	
	cp -r plot_electrons_up.jl V$i;
	cp -r plot_electrons_dn.jl V$i;
	cp -r plot_phonons_phspace.jl V$i;
	cd V$i;
	sed -i '' "s/vcur = 48/vcur = ${i}/g" plot_electrons_up.jl;	
	sed -i '' "s/vcur = 48/vcur = ${i}/g" plot_electrons_dn.jl;
	sed -i '' "s/vcur = 48/vcur = ${i}/g" plot_phonons_phspace.jl;
	cd ../;
done
