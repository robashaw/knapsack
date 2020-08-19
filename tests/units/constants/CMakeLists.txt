blt_add_executable( NAME       const_test
                    SOURCES    driver.f90
                    DEPENDS_ON libknapsack fruit )

blt_add_test( NAME  const_test
	      COMMAND const_test )
