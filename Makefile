# Makefile syntax:
# $< = first dependency of target
# $@ = name of target

release: clusterRelease
dependencies: sdsl ropebwt_lib metacluster_lib dbwt_lib

clusterSrc = UnionFind.cpp tools.cpp Parser.cpp bwt.cpp fasta_tools.cpp K_means_metacluster.cpp K_means.cpp main.cpp
clusterObjects=$(patsubst %.cpp,build/bwtCluster/%.o,$(clusterSrc))
parserObjects=build/bwtCluster/Parser.o build/bwtCluster/tools.o

dbwtSrc = dbwt.c dbwt_queue.c dbwt_utils.c sais.c
cppflags = -O3 -pthread -fopenmp -march=native -std=c++11 -g -Wall -Wextra -Wno-sign-compare -MMD
link = -L ./lib -ldivsufsort64 -lsdsl -ldbwt -lropebwt2 -lz -lmetacluster
includes = -I include -I src/bwtCluster -I src/dataStructures -I src/tools -I ./dbwt -I ropebwt2 -I src/MetaCluster/include

# Search paths for source files
VPATH=src/bwtCluster:src/tools:src/dataStructures

# Include header dependencies
-include $(clusterObjects:%.o=%.d)

clusterRelease: $(clusterObjects)
	$(CXX) $(cppflags) $(clusterObjects) $(link) -o bin/BWT_cluster

	
# Compilation rule for object files in directory build/bwtCluster
build/bwtCluster/%.o : %.cpp build
	$(CXX) -c $< $(cppflags) $(includes) -o $@

# Make target to create the build directory if does not already exist
build:
	@mkdir $@ -p
	@mkdir $@/bwtCluster -p
	
sdsl:
	cd sdsl-lite; sh install.sh; cd ..;
	mkdir -p include
	cp -r sdsl-lite/include/sdsl include/;
	cp sdsl-lite/build/external/libdivsufsort/include/divsufsort.h include
	cp sdsl-lite/build/external/libdivsufsort/include/divsufsort64.h include
	cp sdsl-lite/build/lib/libsdsl.a lib;
	cp sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort.a lib
	cp sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a lib

dbwt_lib:
	cd dbwt; gcc -std=c99 -O3 -g $(dbwtSrc) -c; ar rcs ../lib/libdbwt.a *.o; rm *.o;

ropebwt:
	cd ropebwt2; make;
	cp ropebwt2/ropebwt2 bin

ropebwt_lib:
	cd ropebwt2; gcc -g -O2 *.c -c -lz -lpthread; 
	ar rcs lib/libropebwt2.a ropebwt2/*.o
	rm ropebwt2/*.o

metacluster: $(parserObjects)
	$(CXX) -std=c++11 -fopenmp -I src/MetaCluster/include -I src/tools src/MetaCluster/MetaCluster.cpp src/MetaCluster/USet.cpp src/MetaCluster/ULL4.cpp src/MetaCluster/TomAlgorithm.cpp src/MetaCluster/MetaCluster_main.cpp $(parserObjects) -o bin/metacluster_main

metacluster_lib:
	$(CXX) -std=c++11 -fopenmp -I src/MetaCluster/include -I src/tools src/MetaCluster/MetaCluster.cpp src/MetaCluster/USet.cpp src/MetaCluster/ULL4.cpp src/MetaCluster/TomAlgorithm.cpp -c
	ar rcs lib/libmetacluster.a *.o
	rm *.o

