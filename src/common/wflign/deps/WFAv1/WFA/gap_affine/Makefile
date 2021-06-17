###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build

###############################################################################
# Modules
###############################################################################
MODULES=affine_penalties \
        affine_wavefront \
        affine_table \
        affine_wavefront_align \
        affine_wavefront_backtrace \
        affine_wavefront_display \
        affine_wavefront_extend \
        affine_wavefront_penalties \
        affine_wavefront_reduction \
        affine_wavefront_utils \
        swg \
        wavefront_stats

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

CC_XFLAGS=-march=native ##-fopt-info-vec-optimized

###############################################################################
# Rules
###############################################################################
        
all: $(OBJS)

$(FOLDER_BUILD)/affine_wavefront.o : affine_wavefront.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD)/affine_wavefront_align.o : affine_wavefront_align.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
$(FOLDER_BUILD)/affine_wavefront_extend.o : affine_wavefront_extend.c
	$(CC) $(CC_FLAGS) $(CC_XFLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
	

