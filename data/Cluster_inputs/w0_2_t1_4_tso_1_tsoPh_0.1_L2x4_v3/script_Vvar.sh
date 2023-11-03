#!/bin/bash

for V in {60,65}; #,75,85};
	do mkdir V${V}/;
	cp -r generate_traj0.jl V${V}/;
	cp -r input0P.jl V${V}/;
	cp -r input0M.jl V${V}/;
	cp -r job_submit_scratch0_P.sh V${V}/;
	cp -r job_submit_scratch0_M.sh V${V}/;
	cp -r job_submit_gen0.sh V${V}/; 
	cd V${V}/;
	sed -i "s/V = 0/V = ${V}/g" input0P.jl;
	sed -i "s/V = 0/V = ${V}/g" input0M.jl;
	for j in {1..225};
		do cp -r input0P.jl input${j}P.jl;
		cp -r input0M.jl input${j}M.jl;
		cp -r generate_traj0.jl generate_traj${j}.jl; 
		cp -r job_submit_gen0.sh job_submit_gen${j}.sh
		sed -i "s/version = 0/version = ${j}/g" generate_traj${j}.jl;
		sed -i "s/version = 0/version = ${j}/g" input${j}P.jl;
		sed -i "s/version = 0/version = ${j}/g" input${j}M.jl;
		cp -r job_submit_scratch0_P.sh job_submit_scratch${j}_P.sh;
		cp -r job_submit_scratch0_M.sh job_submit_scratch${j}_M.sh;
		sed -i "s/data0/data${j}/g" job_submit_scratch${j}_P.sh;
		sed -i "s/data0/data${j}/g" job_submit_scratch${j}_M.sh;
		sed -i "s/init_X_traj0/init_X_traj${j}/g" job_submit_scratch${j}_P.sh;
                sed -i "s/init_P_traj0/init_P_traj${j}/g" job_submit_scratch${j}_P.sh;
		sed -i "s/init_X_traj0/init_X_traj${j}/g" job_submit_scratch${j}_M.sh;
		sed -i "s/init_P_traj0/init_P_traj${j}/g" job_submit_scratch${j}_M.sh; 
		sed -i "s/input0/input${j}/g" job_submit_scratch${j}_P.sh;
		sed -i "s/input0/input${j}/g" job_submit_scratch${j}_M.sh;
		sed -i "s/traj0/traj${j}/g" job_submit_gen${j}.sh;
	done
#	rm input0P.jl;
#	rm input0M.jl;
#	rm generate_traj0.jl;
#	rm job_submit_gen0.sh;
#	rm job_submit_scratch0_P.sh;
#	rm job_submit_scratch0_M.sh;
	cd ../;
done
