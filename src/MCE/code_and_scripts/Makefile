MKDIR_P = mkdir -p
BUILD_DIR = build
FINAL_LIB_NAME=simulate.so
MCE_BINARY_START_NAME=MCE_start_binary.out
MKL_LINK_FLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
LINK_FLAGS=${MKL_LINK_FLAGS} ${CUDA_LINK_FLAGS} -lgomp -lpthread -lm -ldl

# for stefans laptop: (cuda can be removed if USE_CUDA == 0 in src/sim_MCE_defs.h
PYTHON_COMPILE_FLAG=-I/usr/include/python3.8
#CUDA_COMPILE_FLAGS=-I/usr/local/cuda/include
#CUDA_LINK_FLAGS=-L/usr/local/cuda/lib64 -lcusolver -lcudart -lcublas -lnvToolsExt
#COMPILE_FLAGS=-g -Wno-unknown-pragmas -DCLUSTER_MODE=0
#COMPILE_FLAGS=-Ofast -march=native -Wno-unknown-pragmas -DCLUSTER_MODE=0
COMPILE_FLAGS=-fopenmp -Ofast -march=native -DCLUSTER_MODE=0

# for rocks cluster:
#PYTHON_COMPILE_FLAG=-I/usr/include/python3.6m
CUDA_COMPILE_FLAGS=
CUDA_LINK_FLAGS=
#COMPILE_FLAGS=-fopenmp -Ofast -DCLUSTER_MODE=1

# end of alternative compile flags

py_lib: ${BUILD_DIR}/simulate.o ${BUILD_DIR}/sim_MTE.o ${BUILD_DIR}/sim_MCE.o ${BUILD_DIR}/sim_MCEv1_1e.o ${BUILD_DIR}/sim_MCEv2_1e.o ${BUILD_DIR}/sim_MCE_util.o ${BUILD_DIR}/zpotrf_omp.o
	gcc -shared -O1 $^ -o ${FINAL_LIB_NAME} ${LINK_FLAGS}

mce_profiler: src/mce_profiler.c ${BUILD_DIR}/sim_MCEv1_1e.o directories ${BUILD_DIR}/zpotrf_omp.o
	gcc ${COMPILE_FLAGS} -c src/mce_profiler.c -o ${BUILD_DIR}/mce_profiler.o
	gcc -O1 ${BUILD_DIR}/mce_profiler.o ${BUILD_DIR}/sim_MCEv1_1e.o ${BUILD_DIR}/zpotrf_omp.o ${LINK_FLAGS}

# binary start for server usage. probably no CUDA support
mce_binary_start: src/MCE_binary_start.c ${BUILD_DIR}/sim_MCEv1_1e.o directories ${BUILD_DIR}/zpotrf_omp.o
	gcc ${COMPILE_FLAGS} -c src/MCE_binary_start.c -o ${BUILD_DIR}/MCE_binary_start.o
	gcc -o ${MCE_BINARY_START_NAME} -O1 ${BUILD_DIR}/MCE_binary_start.o ${BUILD_DIR}/sim_MCEv1_1e.o ${BUILD_DIR}/zpotrf_omp.o ${LINK_FLAGS}


.PHONY: directories
.PHONY: clean
all: py_lib

directories: ${BUILD_DIR}

${BUILD_DIR}:
	${MKDIR_P} ${BUILD_DIR}

src/simulate.c : src/simulate.pyx
	cythonize src/simulate.pyx

${BUILD_DIR}/simulate.o : src/simulate.c directories
	gcc -Wall ${COMPILE_FLAGS} -m64 -pthread -fwrapv -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC ${PYTHON_COMPILE_FLAG} -c src/simulate.c -o ${BUILD_DIR}/simulate.o -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

${BUILD_DIR}/%.o: src/%.c directories
	gcc -Wall ${COMPILE_FLAGS} -m64 -I"${MKLROOT}/include" -fPIC -c $< -o $@


clean:
	rm -rf build/
	rm -f ${FINAL_LIB_NAME}
	rm -f src/simulate.c
	rm -rf __pycache__
	rm -rf src/__pycache__
	rm -f a.out
