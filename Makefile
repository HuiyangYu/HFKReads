BIN_DIR := $(shell pwd)/bin
SRC_DIR := $(shell pwd)/src
INCLUDE_DIR := $(shell pwd)/include

.PHONY: all clean

all: $(BIN_DIR)/hfkreads

$(BIN_DIR)/hfkreads: $(SRC_DIR)/HFKReads.cpp | $(BIN_DIR)
	g++ --std=c++11 -g -O3 $< -lz -pthread -I$(INCLUDE_DIR) -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)


