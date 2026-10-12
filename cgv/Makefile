CC ?= gcc
CFLAGS ?= -O3 -march=native -Wall -Wextra -Wno-unused-parameter -g
# always (also with CFLAGS given): no fused multiply-adds, the same LLL basis everywhere
override CFLAGS += -ffp-contract=off
NVCC ?= $(HOME)/miniforge3/envs/llamacpp-cuda/bin/nvcc
CUDA_ARCH ?= sm_120
HIPCC ?= hipcc
HIP_ARCH ?= gfx1151
NVFLAGS = -O3 -g -arch=$(CUDA_ARCH)
HIPFLAGS = -O3 -x hip --offload-arch=$(HIP_ARCH)

.DEFAULT_GOAL := cgv

# flag stamps: objects are rebuilt when the compiler, CFLAGS or the GPU arch change
stamp = $(shell echo '$(2)' | cmp -s - $(1) 2>/dev/null || echo '$(2)' > $(1))
$(call stamp,.cflags,$(CC) $(CFLAGS))
$(call stamp,.nvflags,$(NVCC) $(NVFLAGS))
$(call stamp,.hipflags,$(HIPCC) $(HIPFLAGS))

# gv.c is built once per lane count n (2, 3, 4; chosen at run time in main.c) and key width:
# gv_nl<n> 128-bit keys, gvw_nl<n> 256, gvx_nl<n> 512; gvg*_nl<n> the same with the GPU extraction (-g);
# gpu.cu likewise (gpu*_nl<n> CUDA, hip*_nl<n> HIP; w/x: 256/512-bit global keys)
NLS = 2 3 4
W_ =
W_w = -DCGV_WIDE=256
W_x = -DCGV_WIDE=512
GPUDEF = -Dgpu_extract=gpu_extract_nl$*$(1) -Dgpu_free_gb=gpu_free_gb_nl$*$(1)
define WIDTH
gv$(1)_nl%.o: gv.c cgv.h .cflags
	$$(CC) $$(CFLAGS) -DNL=$$* $$(W_$(1)) -DCGV_ENTRY=cgv_entry_nl$$*$(1) -c $$< -o $$@
gvg$(1)_nl%.o: gv.c cgv.h gpu.h .cflags
	$$(CC) $$(CFLAGS) -DUSE_GPU $$(call GPUDEF,$(1)) -DNL=$$* $$(W_$(1)) -DCGV_ENTRY=cgv_entry_nl$$*$(1) -c $$< -o $$@
gpu$(1)_nl%.o: gpu.cu gpu.h gpu_compat.h .nvflags
	$$(NVCC) $$(NVFLAGS) -DNL=$$* $$(W_$(1)) $$(call GPUDEF,$(1)) -c $$< -o $$@
hip$(1)_nl%.o: gpu.cu gpu.h gpu_compat.h .hipflags
	$$(HIPCC) $$(HIPFLAGS) -DNL=$$* $$(W_$(1)) $$(call GPUDEF,$(1)) -c $$< -o $$@
endef
$(eval $(call WIDTH,))
$(eval $(call WIDTH,w))
$(eval $(call WIDTH,x))
ALLW = $(foreach w,_ w_ x_,$(NLS:%=$(1)$(w:_=)_nl%.o))
main.o: main.c cgv.h .cflags
	$(CC) $(CFLAGS) -c $< -o $@

# cgv: CPU build
cgv: main.o $(call ALLW,gv)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread -lm

# cgv_gpu: same, with the extraction on a CUDA GPU (-g [device])
cgv_gpu: main.o $(call ALLW,gvg) $(call ALLW,gpu)
	$(NVCC) -arch=$(CUDA_ARCH) -o $@ $^ -lpthread -lm

# cgv_hip: the same GPU code built for AMD with HIP (32-lane waves: RDNA, e.g. gfx1151)
cgv_hip: main.o $(call ALLW,gvg) $(call ALLW,hip)
	$(HIPCC) --offload-arch=$(HIP_ARCH) -o $@ $^ -lpthread -lm

# cgv_maj: the same pipeline over upward-rounded nonnegative reals, for proven |GV| bounds
cgv_maj: gv.c cgv.h .cflags
	$(CC) $(CFLAGS) -frounding-math -DNL=1 -DCGV_MAJ -o $@ gv.c -lpthread -lm

clean:
	rm -f cgv cgv_gpu cgv_hip cgv_maj *.o .cflags .nvflags .hipflags

.PHONY: clean
