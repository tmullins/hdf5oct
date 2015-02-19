CXX=h5c++
src= h5read.cc
headers=h5read.h
octs=$(src:.cc=.oct)
objs=$(src:.cc=.o)
testobjs=test/write_test.o
MKOCTFILE=CXX=$(CXX) mkoctfile -g

VERSION=0.2.0
PACKAGEFILE=hdf5oct-$(VERSION).tar.gz

.PHONY: test clean install uninstall package

all: $(octs) package

%.oct: $(objs)
	$(MKOCTFILE) -o $@ $(objs)

%.o: %.cc $(headers)
	$(MKOCTFILE) -c $<

clean:
	rm -f *.o *.oct package/inst/* $(testobjs) test/write_test test/test*.h5 $(PACKAGEFILE)

install: $(PACKAGEFILE)
	@echo "-- Install Octave Package ------------"
	octave --silent --no-gui --eval "pkg install $(PACKAGEFILE)"

uninstall:
	@echo "-- Uninstall Octave Package ----------"
	octave --silent --no-gui --eval "pkg uninstall hdf5oct"

package: $(PACKAGEFILE)

$(PACKAGEFILE): $(octs)
	@echo "-- Create Octave Package Archive ------------"
	mkdir -p package/inst
	cp *.oct package/inst
	tar -czf $(PACKAGEFILE) package/

# TESTING ###########

# a minimal program to generate some testdata
test/write_test: test/write_test.cc
	$(CXX) -o $@ $< -lhdf5

# a target to test the octave functions
test: test/write_test
	@echo "-- Perform Tests --------------"
	rm -f test/test*.h5
	cd test && ./write_test
	cd test && octave --silent --no-gui h5test.m
