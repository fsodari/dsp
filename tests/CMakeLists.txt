#[[ ================================================================================
    MIT License

    Copyright (c) 2019 Frank Sodari

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
==================================================================================]]

find_package(CppUTest REQUIRED)

add_executable(run_all_tests EXCLUDE_FROM_ALL run_all_tests.cpp)

# Link the tests here
target_link_libraries(run_all_tests matrix_test)

target_link_libraries(run_all_tests ${CppUTest_LIBRARIES})

add_custom_target(check DEPENDS run_all_tests)
add_custom_command(TARGET check POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND run_all_tests)
