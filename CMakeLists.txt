include(CheckCCompilerFlag)

cmake_minimum_required(VERSION 2.6)

project(TABU)

if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${TABU_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TABU_SOURCE_DIR})

macro(c_compiler_flag_test flag_name flag fallback_flag required)
  check_c_compiler_flag(${flag} ${flag_name})
  if("${flag_name}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  elseif(NOT "${fallback_flag}" STREQUAL "")
    c_compiler_flag_test(${flag_name}_FALLBACK ${fallback_flag} "" ${required})
  elseif("${required}")
    message(FATAL_ERROR "${flag} support is required by ${CMAKE_C_COMPILER}")
  endif()
endmacro(c_compiler_flag_test)

c_compiler_flag_test(C_STD_GNU99 "-std=gnu99" "" 1)
c_compiler_flag_test(C_PEDANTIC "-pedantic" "" 0)
c_compiler_flag_test(C_MARCH "-march=native" "" 0)
c_compiler_flag_test(C_CLANG_ASSEMBLER "-Wa,-q" "" 0)

if (CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
  c_compiler_flag_test(C_DEBUG "-g" "" 1)
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "SAFE")
  c_compiler_flag_test(C_O2 "-O2" "" 0)
else()
  c_compiler_flag_test(C_FUNROLL_LOOPS "-funroll-loops" "" 0)
  c_compiler_flag_test(C_OFAST "-Ofast" "-O3" 1)
  c_compiler_flag_test(C_FAST_MATH "-ffast-math" "" 0)
endif()

find_library(RNA libRNA.a PATHS)

add_executable(get_barrier get_barrier.c fold_vars.c foldMod.c utils.c energy_par.c read_epars.c params.c)

target_link_libraries(get_barrier ${RNA})
target_link_libraries(get_barrier m)
install(TARGETS get_barrier DESTINATION bin)
