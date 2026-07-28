# CMake generated Testfile for 
# Source directory: /home/viventus/WHQ/src/WHQ5402/tests/ctests
# Build directory: /home/viventus/WHQ/src/WHQ5402/src/tests/ctests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fortran_ctest_should_fail "/home/viventus/WHQ/src/WHQ5402/src/tests/ctests/fortran_ctest_should_fail")
set_tests_properties(fortran_ctest_should_fail PROPERTIES  WILL_FAIL "TRUE" _BACKTRACE_TRIPLES "/home/viventus/WHQ/src/WHQ5402/tests/ctests/CMakeLists.txt;37;add_test;/home/viventus/WHQ/src/WHQ5402/tests/ctests/CMakeLists.txt;0;")
add_test(fortran_ctest_should_not_fail "/home/viventus/WHQ/src/WHQ5402/src/tests/ctests/fortran_ctest_should_not_fail")
set_tests_properties(fortran_ctest_should_not_fail PROPERTIES  _BACKTRACE_TRIPLES "/home/viventus/WHQ/src/WHQ5402/tests/ctests/CMakeLists.txt;44;add_test;/home/viventus/WHQ/src/WHQ5402/tests/ctests/CMakeLists.txt;0;")
