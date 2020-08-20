blt_add_executable( NAME       random_test
                    SOURCES    driver.f90
                    DEPENDS_ON libknapsack fruit )

blt_add_test( NAME  random_test
	      COMMAND random_test )
