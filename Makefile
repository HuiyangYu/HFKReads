BIN_DIR := $(shell pwd)/bin
LIB_DIR := $(shell pwd)/lib
SRC_DIR := $(shell pwd)/src
INCLUDE_DIR := $(shell pwd)/include

OS := $(shell uname)

ifeq ($(OS),Darwin)
    ZLIB_LIBS := -lz
else
    ifeq ($(shell pkg-config --exists zlib && echo 1), 1)
        ZLIB_LIBS := $(shell pkg-config --libs zlib)
    else
        ZLIB_LIBS := -L$(LIB_DIR) -lz
    endif
endif

.PHONY: all clean

all: $(BIN_DIR)/hfkreads

$(BIN_DIR)/hfkreads: $(SRC_DIR)/HFKReads.cpp | $(BIN_DIR)
	g++ --std=c++11 -g -O3 $< $(ZLIB_LIBS) -pthread -I$(INCLUDE_DIR) -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)

