B. Johnson JCSDA 10/2023

Synopsis:
Application, Unit, and Regression tests, largely culled from Paul van Delst's and Dave Groff's CRTM tests, modified to work with CRTM v3.0.1 in a CMake environment.
Not a complete or comprehensive suite of tests, please add a test each time you add a new code element or substantially change code.


Layout:

  cmake/                 (cmake Modules directory for storing the required cmake modules for the testing framework, do not delete!)
  CMakeLists.txt         (primary CMakeLists.txt file, tests, parameters, etc.  are defined here)
  mains/                 (top level of source files)
    application/            (Application tests go here, larger, end-to-end, more punishing tests [large scope ])
    regression/             (Regression tests go here to check various functionalities of CRTM   [medium scope])
    unit/                   (Unit tests go here to test small units of code, convergence, etc.   [small scope ])
  readme_crtm_tests.txt  (this file)

Instructions:

  Tests do not have to have a success metric, but it's nice to have in order to test for failure.
  STOP 0 and STOP 1 signal ctest for success or failure, respectively.  Any "failure" messages or "error" messages should have a STOP 1 following it.
	Any existing STOP should be converted to a STOP 1 as well.   


Prerequisites:

  Make and _install_ the CRTM library following the directions provided by the CRTM package.   
  These tests use both the generated modules and the libcrtm.a, so it needs to see both of these.
  The cmake step during building the CRTM library will also create the necessary makefile for building tests.  

Running Tests:
  cd ./build   (directory where CRTM was built using cmake)
  cmake ..; make -j12    (if not already )
  ctest -j4    (run 4 tests in parallel)

Slightly more verbose Mode (-VV):

  ls -l  (find your test name)
  ctest -VV -R testname (run an individual test)

Cleanup:

  In the build/ directory, you may need to "make clean" to remove compiled executables, modules, etc.  This will not remove the binary data (coefficient files).


Troubleshooting/Support:
  Please feel free to contact us at:
    https://forums.jcsda.org/
    or 
    crtm-support@groups.google.com
  
