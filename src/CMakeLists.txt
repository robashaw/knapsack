blt_add_library( NAME     libknapsack
		 SOURCES  constants.f90 hashtable.f90 random.f90 ioutil.f90 radiative.f90 knapsack.f90 sample.f90 nonradiative.f90)

blt_add_executable( NAME knapsack
		    SOURCES main.f90
		    DEPENDS_ON libknapsack )
		    
install(TARGETS knapsack DESTINATION bin)

