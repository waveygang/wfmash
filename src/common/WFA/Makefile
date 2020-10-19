###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BIN=bin
FOLDER_BUILD=build

UNAME=$(shell uname)

CC=gcc
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
SUBDIRS=benchmark \
        edit \
        gap_affine \
        gap_lineal \
        system \
        utils
       
LIB_WFA=$(FOLDER_BUILD)/libwfa.a

all: CC_FLAGS+=-O3
all: MODE=all
all: setup
all: $(SUBDIRS) tools $(LIB_WFA)

$(FOLDER_BUILD)/*.o: $(SUBDIRS)

debug: setup
debug: MODE=all
debug: $(SUBDIRS) tools $(LIB_WFA)

$(LIB_WFA): $(FOLDER_BUILD)/*.o
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null

setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD)

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
tools:
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS) tools
.FORCE:

