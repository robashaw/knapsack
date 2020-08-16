blt_add_executable( NAME     knapsack
		    SOURCES  constants.f90 hashtable.f90 random.f90 ioutil.f90 radiative.f90 knapsack.f90 sample.f90 nonradiative.f90 main.f90)

install(TARGETS knapsack DESTINATION bin)

