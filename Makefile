CXXFLAGS := -O2 -std=c++20 $(shell fastjet-config --cxxflags) $(shell pythia8-config --cxxflags) $(shell root-config --cflags)
LDLIBS   := $(shell pythia8-config --libs) $(shell root-config --libs) $(shell fastjet-config --libs)

bin/runPythia: src/runPythia.cc
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -rf bin
