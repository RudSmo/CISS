# Running on the cluster
## compiling
for running on the cluster, CLUSTER_MODE should/must be enabled in src/sim_MCEv1_1e.c (sigaction unknown on cluster)
also check MAX_NUM_THREADS and OVERRIDE_OMP_NUM_THREADS.
Makefile must not use the cuda flags (espically th linker flag)

rsync -vax src stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/

ssh itp

ssh ds10  (ds9 doesn't have module devtoolset)

cd /net/theorie/rocks/stefanarne.graeber


module load devtoolset
module load intel/mkl/intel64/lp64/2018.0.128
module load rh-python

or short:
source load_compile_modules.sh

make


## running simulation

from itp:
ssh rocks

edit jobscript

qsub jobscript

OR to run it with python:
module load rh-python

## getting data back

rsync -vax stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/data_export/ data_export/

