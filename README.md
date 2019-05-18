Buind Status.
-------------
[![Build Status](https://travis-ci.com/sht-cs133/hw5-6-mathematical-optimization-YuXinFan.svg?token=appyqywAiysphxXppp9y&branch=master)](https://travis-ci.com/sht-cs133/hw5-6-mathematical-optimization-YuXinFan)

## Notices

1. This time we provided 30 test cases under `test/` directory. During implementing
   you might not want to compile all these test cases. To filter these cases by
   a name template, you can pass a *regular expression* variable `ENABLE_TESTS` to `cmake`.
   
   E.g. to only compile the test cases in `test/test_po_gd_an.cpp`,
   `test/test_po_gd_nc.cpp` and `test/test_po_gd_nf.cpp`, we set the name template
   as `po_gd[a-z_]*`:
   ```bash
   mkdir build && cd build
   cmake -DENABLE_TESTS="po_gd[a-z_]*" ..
   make -j `nproc --all` all && make test
   ```
   
   Please note that the name to be matched is extracted from *file name* of these test cases,
   e.g. the name of test case `test/test_po_gd_an.cpp` is `po_gd_an`. The provided
   `ENABLE_TESTS` regexp is then matched with these names, and only matched tests
   will be compiled and activated.
   
   One more notice: all test cases are enabled in Travis-CI by default. You can
   find the syntax of CMake regexp [here](https://cmake.org/cmake/help/latest/command/string.html#regex-specification).
2. You are encouraged to add new test cases for your own implementation details,
   e.g. helper functions and numerical differentiations. Don't forget to add them
   to `CMakeList.txt`.


