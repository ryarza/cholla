# CMake generated Testfile for 
# Source directory: /home/rcastroy/src/cholla/tests/cmake_stuff
# Build directory: /home/rcastroy/src/cholla/tests/cmake_stuff
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(poisson64 "../../cholla.sor" "poissonParameterFiles/poisson64.txt")
set_tests_properties(poisson64 PROPERTIES  WILL_FAIL "FALSE" _BACKTRACE_TRIPLES "/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;8;add_test;/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;0;")
add_test(poisson128 "../../cholla.sor" "poissonParameterFiles/poisson128.txt")
set_tests_properties(poisson128 PROPERTIES  WILL_FAIL "FALSE" _BACKTRACE_TRIPLES "/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;11;add_test;/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;0;")
add_test(poisson256 "../../cholla.sor" "poissonParameterFiles/poisson256.txt")
set_tests_properties(poisson256 PROPERTIES  WILL_FAIL "FALSE" _BACKTRACE_TRIPLES "/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;14;add_test;/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;0;")
add_test(poisson512 "../../cholla.sor" "poissonParameterFiles/poisson512.txt")
set_tests_properties(poisson512 PROPERTIES  WILL_FAIL "FALSE" _BACKTRACE_TRIPLES "/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;17;add_test;/home/rcastroy/src/cholla/tests/cmake_stuff/CMakeLists.txt;0;")
