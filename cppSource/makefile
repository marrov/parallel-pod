ifeq ($(OS),Windows_NT)
    uname_S := Windows
else
    uname_S := $(shell uname -s)
endif

ifeq ($(uname_S), Windows)
    TARGETS = 3DPOD.exe
endif
ifeq ($(uname_S), Linux)
    TARGETS = 3DPOD.out
endif
#ifeq ($(uname_S), .....) #See https://stackoverflow.com/a/27776822/938111
#    target = .....
#endif

IDIR    := ../testCase/eigen/ ./depends/ezOptionParser/
TEST_DIR:= ../testCase/
CC      := g++
# Substitiutions used below. See https://www.gnu.org/software/make/manual/html_node/Substitution-Refs.html
CCFLAGS :=-std=c++11 $(IDIR:%=-I%) -O3
LDFLAGS :=

#TARGETS := 3DPOD.out

.PHONY: all clean install_to_test

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

install_to_test: all
	install $(TARGETS) $(TEST_DIR)

$(TARGETS): 3DPOD.cpp
	$(CC) -o $@ $^ $(CCFLAGS)

