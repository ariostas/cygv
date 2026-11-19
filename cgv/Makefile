CC ?= gcc
CFLAGS ?= -O3 -march=native -Wall -Wextra -Wno-unused-parameter -g
NVCC ?= $(HOME)/miniforge3/envs/llamacpp-cuda/bin/nvcc
CUDA_ARCH ?= sm_120

# cgv: CPU build, lane count chosen at run time (main.c)
cgv: main.c gv.c cgv.h
	for n in 2 3 4; do $(CC) $(CFLAGS) -DNL=$$n -DCGV_ENTRY=cgv_entry_nl$$n -c gv.c -o gv_nl$$n.o || exit 1; done
	$(CC) $(CFLAGS) -c main.c -o main.o
	$(CC) $(CFLAGS) -o $@ main.o gv_nl2.o gv_nl3.o gv_nl4.o -lpthread -lm

# cgv_gpu: same, with the extraction on a CUDA GPU (-g [device])
cgv_gpu: main.c gv.c cgv.h gpu.cu gpu.h
	for n in 2 3 4; do \
	  $(CC) $(CFLAGS) -DUSE_GPU -DNL=$$n -DCGV_ENTRY=cgv_entry_nl$$n -Dgpu_extract=gpu_extract_nl$$n -Dgpu_free_gb=gpu_free_gb_nl$$n -c gv.c -o gvg_nl$$n.o || exit 1; \
	  $(NVCC) -O3 -g -arch=$(CUDA_ARCH) -DNL=$$n -Dgpu_extract=gpu_extract_nl$$n -Dgpu_free_gb=gpu_free_gb_nl$$n -c gpu.cu -o gpu_nl$$n.o || exit 1; \
	done
	$(CC) $(CFLAGS) -c main.c -o main.o
	$(NVCC) -arch=$(CUDA_ARCH) -o $@ main.o gvg_nl2.o gvg_nl3.o gvg_nl4.o gpu_nl2.o gpu_nl3.o gpu_nl4.o -lpthread -lm

# cgv_maj: the same pipeline over upward-rounded nonnegative reals, for proven |GV| bounds
cgv_maj: gv.c cgv.h
	$(CC) $(CFLAGS) -frounding-math -DNL=1 -DCGV_MAJ -o $@ gv.c -lpthread -lm

# cgv_hip: the same GPU code built for AMD with HIP (32-lane waves: RDNA, e.g. gfx1151)
HIPCC ?= hipcc
HIP_ARCH ?= gfx1151
cgv_hip: main.c gv.c cgv.h gpu.cu gpu.h gpu_compat.h
	for n in 2 3 4; do \
	  $(CC) $(CFLAGS) -DUSE_GPU -DNL=$$n -DCGV_ENTRY=cgv_entry_nl$$n -Dgpu_extract=gpu_extract_nl$$n -Dgpu_free_gb=gpu_free_gb_nl$$n -c gv.c -o gvh_nl$$n.o || exit 1; \
	  $(HIPCC) -O3 -x hip --offload-arch=$(HIP_ARCH) -DNL=$$n -Dgpu_extract=gpu_extract_nl$$n -Dgpu_free_gb=gpu_free_gb_nl$$n -c gpu.cu -o hip_nl$$n.o || exit 1; \
	done
	$(CC) $(CFLAGS) -c main.c -o main.o
	$(HIPCC) --offload-arch=$(HIP_ARCH) -o $@ main.o gvh_nl2.o gvh_nl3.o gvh_nl4.o hip_nl2.o hip_nl3.o hip_nl4.o -lpthread -lm

clean:
	rm -f cgv cgv_gpu cgv_hip cgv_maj *.o
