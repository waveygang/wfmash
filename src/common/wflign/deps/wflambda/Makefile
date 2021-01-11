###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BIN=bin
FOLDER_BUILD=build
FOLDER_ROOT=./

UNAME=$(shell uname)

CC=g++
CPP=g++

LD_FLAGS=-lm
CC_FLAGS=-Wall -g
ifeq ($(UNAME), Linux)
  LD_FLAGS+=-lrt 
endif

AR=ar
AR_FLAGS=-rsc

###############################################################################
# Compile rules
###############################################################################

LIB_WFA=$(FOLDER_BUILD)/libwflambda.a

all: CC_FLAGS+=-O3 -Wno-pointer-arith -Wno-literal-suffix
all: MODE=all
all: setup
all: lib_wfa

$(FOLDER_BUILD): setup

debug: setup
debug: MODE=all
#debug: $(SUBDIRS) lib_wfa tools

OBJS=$(FOLDER_BUILD)/commons.o \
	$(FOLDER_BUILD)/vector.o \
	$(FOLDER_BUILD)/string_padded.o \
	$(FOLDER_BUILD)/affine_wavefront_backtrace.o \
	$(FOLDER_BUILD)/affine_wavefront_reduction.o \
	$(FOLDER_BUILD)/affine_wavefront_penalties.o \
	$(FOLDER_BUILD)/affine_penalties.o \
	$(FOLDER_BUILD)/affine_wavefront_utils.o \
	$(FOLDER_BUILD)/affine_table.o \
	$(FOLDER_BUILD)/wavefront_stats.o \
	$(FOLDER_BUILD)/swg.o \
	$(FOLDER_BUILD)/affine_wavefront_display.o \
	$(FOLDER_BUILD)/affine_wavefront_align.o \
	$(FOLDER_BUILD)/affine_wavefront.o \
	$(FOLDER_BUILD)/affine_wavefront_extend.o \
	$(FOLDER_BUILD)/nw.o \
	$(FOLDER_BUILD)/benchmark_edit.o \
	$(FOLDER_BUILD)/benchmark_utils.o \
	$(FOLDER_BUILD)/benchmark_gap_lineal.o \
	$(FOLDER_BUILD)/benchmark_gap_affine.o \
	$(FOLDER_BUILD)/profiler_timer.o \
	$(FOLDER_BUILD)/mm_allocator.o \
	$(FOLDER_BUILD)/profiler_counter.o \
	$(FOLDER_BUILD)/edit_dp.o \
	$(FOLDER_BUILD)/edit_cigar.o \
	$(FOLDER_BUILD)/edit_table.o \

$(FOLDER_BUILD)/dna_text.o: ./utils/dna_text.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/commons.o: ./utils/commons.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/vector.o: ./utils/vector.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/string_padded.o: ./utils/string_padded.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

CC_XFLAGS=-march=native ##-fopt-info-vec-optimized

$(FOLDER_BUILD)/affine_wavefront_backtrace.o: ./gap_affine/affine_wavefront_backtrace.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_reduction.o: ./gap_affine/affine_wavefront_reduction.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_penalties.o: ./gap_affine/affine_wavefront_penalties.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_penalties.o: ./gap_affine/affine_penalties.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_utils.o: ./gap_affine/affine_wavefront_utils.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_table.o: ./gap_affine/affine_table.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/wavefront_stats.o: ./gap_affine/wavefront_stats.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/swg.o: ./gap_affine/swg.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_display.o: ./gap_affine/affine_wavefront_display.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_align.o: ./gap_affine/affine_wavefront_align.cpp
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront.o: ./gap_affine/affine_wavefront.cpp
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/affine_wavefront_extend.o: ./gap_affine/affine_wavefront_extend.cpp
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/nw.o: ./gap_lineal/nw.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/benchmark_edit.o: ./benchmark/benchmark_edit.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/benchmark_utils.o: ./benchmark/benchmark_utils.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/benchmark_gap_lineal.o: ./benchmark/benchmark_gap_lineal.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/benchmark_gap_affine.o: ./benchmark/benchmark_gap_affine.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/profiler_timer.o: ./system/profiler_timer.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/mm_allocator.o: ./system/mm_allocator.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/profiler_counter.o: ./system/profiler_counter.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/edit_dp.o: ./edit/edit_dp.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/edit_cigar.o: ./edit/edit_cigar.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@

$(FOLDER_BUILD)/edit_table.o: ./edit/edit_table.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@



setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD)

lib_wfa: $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD)

###############################################################################
# Subdir rule
###############################################################################
#export
#$(SUBDIRS):
#	$(MAKE) --directory=$@ all

#tools: $(SUBDIRS)
#	$(MAKE) --directory=$@ $(MODE)

#.PHONY: $(SUBDIRS) tools

