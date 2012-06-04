.PHONY: all examples html memtest test clean

APP = ezOptionParser
DEVDIR = ezoptionparser-code

all: examples

examples:
	cd examples && $(MAKE) all

clean:
	cd examples && $(MAKE) clean

html:
	cd examples && $(MAKE) html
	pygmentize -O full,linenos=1,style=manni -o html/ezOptionParser.html ezOptionParser.hpp

memtest:
	cd examples && $(MAKE) memtest

test:
	cd examples && $(MAKE) test

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.1.0"
  else
		cd ..; rm -fr $(APP)-$(VER); mkdir $(APP)-$(VER); rsync --cvs-exclude -a $(DEVDIR)/ $(APP)-$(VER); tar zcvf $(APP)-$(VER).tar.gz $(APP)-$(VER); rm -fr $(APP)-$(VER);
  endif
