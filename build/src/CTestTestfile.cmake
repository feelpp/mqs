# CMake generated Testfile for 
# Source directory: /home/anki/csmi/mqs/src
# Build directory: /home/anki/csmi/mqs/build/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(feelpp_mqs_feelpp2d-default-np-6 "/usr/bin/mpiexec" "-n" "6" "/home/anki/csmi/mqs/build/src/feelpp_mqs_feelpp2d" "--config-file" "feelpp2d.cfg")
set_tests_properties(feelpp_mqs_feelpp2d-default-np-6 PROPERTIES  ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0;LSAN_OPTIONS=suppressions=/home/anki/csmi/mqs/feelpp/tools/lsan/suppressions.txt" _BACKTRACE_TRIPLES "/usr/share/feelpp/feel/cmake/modules/feelpp.macros.cmake;199;add_test;/home/anki/csmi/mqs/src/CMakeLists.txt;6;feelpp_add_application;/home/anki/csmi/mqs/src/CMakeLists.txt;0;")
add_test(feelpp_mqs_feelpp2d-default-np-1 "/home/anki/csmi/mqs/build/src/feelpp_mqs_feelpp2d" "--config-file" "feelpp2d.cfg")
set_tests_properties(feelpp_mqs_feelpp2d-default-np-1 PROPERTIES  ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0;LSAN_OPTIONS=suppressions=/home/anki/csmi/mqs/feelpp/tools/lsan/suppressions.txt" _BACKTRACE_TRIPLES "/usr/share/feelpp/feel/cmake/modules/feelpp.macros.cmake;204;add_test;/home/anki/csmi/mqs/src/CMakeLists.txt;6;feelpp_add_application;/home/anki/csmi/mqs/src/CMakeLists.txt;0;")
