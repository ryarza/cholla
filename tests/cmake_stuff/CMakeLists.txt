project(cholla)

add_executable(cholla ../../src/main.cpp)


enable_testing()

add_test(poisson64 ../../cholla.sor poissonParameterFiles/poisson64.txt)
set_tests_properties(poisson64 PROPERTIES WILL_FAIL FALSE)

add_test(poisson128 ../../cholla.sor poissonParameterFiles/poisson128.txt)
set_tests_properties(poisson128 PROPERTIES WILL_FAIL FALSE)

add_test(poisson256 ../../cholla.sor poissonParameterFiles/poisson256.txt)
set_tests_properties(poisson256 PROPERTIES WILL_FAIL FALSE)

add_test(poisson512 ../../cholla.sor poissonParameterFiles/poisson512.txt)
set_tests_properties(poisson512 PROPERTIES WILL_FAIL FALSE)
