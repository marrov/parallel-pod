IDIR    := ../testCase/eigen/
TEST_DIR:= ../testCase/
CC      := g++
CCFLAGS :=-std=c++11 -I$(IDIR)
LDFLAGS :=

TARGETS := 3DPOD.out

.PHONY: all clean install_to_test

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

install_to_test: all
	install ./3DPOD.out $(TEST_DIR)
 

$(TARGETS): 3DPOD.cpp
	$(CC) -o $@ $^ $(CCFLAGS)

