ifeq ($(OS),Windows_NT)
    uname_S := Windows
else
    uname_S := $(shell uname -s)
endif

ifeq ($(uname_S), Windows)
    EXECEND = .exe
endif
ifeq ($(uname_S), Linux)
    EXECEND = .out
endif
#ifeq ($(uname_S), .....) #See https://stackoverflow.com/a/27776822/938111
#    target = .....
#endif

IDIR    := ./depends/eigen/ ./depends/ezOptionParser/
TEST_DIR:= ../testCase/
CC      := g++
# Substitiutions used below. See https://www.gnu.org/software/make/manual/html_node/Substitution-Refs.html
CCFLAGS :=-std=c++11 $(IDIR:%=-I%) -O3 -fopenmp
LDFLAGS :=

NAMES   := 3DPOD_U 3DPOD_s 3DPOD
TARGETS := $(NAMES:%=%$(EXECEND))

.PHONY: all clean install_to_test

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

install_to_test: all
	install $(TARGETS) $(TEST_DIR)

# Pattern rule used below. See https://www.gnu.org/software/make/manual/html_node/Pattern-Intro.html#Pattern-Intro
%$(EXECEND): %.cpp
	$(CC) -o $@ $^ $(CCFLAGS)
    
