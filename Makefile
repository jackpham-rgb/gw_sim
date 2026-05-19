CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -Isrc -pthread
LDFLAGS  ?= -pthread

# raylib flags (only needed for the 3d target)
RAYLIB_LIBS ?= -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

CORE_OBJ_2D := src/binary_system.o
CORE_OBJ_3D := src/nbody.o

BINS_2D := gw_sim gw_visualize
BIN_3D  := gw_3d

.PHONY: all 2d 3d run visualize visualize3d clean help

all: 2d
	@echo ""
	@echo "Built 2D targets. To build the 3D visualizer:"
	@echo "    sudo apt install libraylib-dev"
	@echo "    make 3d"

2d: $(BINS_2D)

3d: $(BIN_3D)

# ---- 2D / data export targets ------------------------------------------
gw_sim: $(CORE_OBJ_2D) src/main.o
	@mkdir -p output
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

gw_visualize: $(CORE_OBJ_2D) src/main_visualize.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ---- 3D raylib target --------------------------------------------------
gw_3d: $(CORE_OBJ_3D) src/main_3d.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(RAYLIB_LIBS)

# ---- Generic compile rule ----------------------------------------------
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- Convenience runners -----------------------------------------------
run: gw_sim
	@mkdir -p output
	./gw_sim

visualize: gw_visualize
	./gw_visualize

visualize3d: gw_3d
	./gw_3d

clean:
	rm -f src/*.o $(BINS_2D) $(BIN_3D) test_nbody
	rm -rf output

help:
	@echo "Targets:"
	@echo "  make            -- build the two 2D binaries"
	@echo "  make 2d         -- same as 'make'"
	@echo "  make 3d         -- build the raylib 3D visualizer (needs libraylib-dev)"
	@echo "  make run        -- run gw_sim, write CSV + WAV"
	@echo "  make visualize  -- run gw_visualize, ASCII orbit + chirp"
	@echo "  make visualize3d-- run gw_3d, real 3D window with scenario switcher"
	@echo "  make clean      -- remove build artifacts and output/"
