add_executable(knapsack constants.f90 ioutil.f90 hashtable.f90 radiative.f90 knapsack.f90 sample.f90 nonradiative.f90 random.f90 progress.f90 main.f90)
set_target_properties(knapsack PROPERTIES LINKER_LANGUAGE Fortran)
install(TARGETS knapsack DESTINATION bin)

