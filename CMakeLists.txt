cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

# check dependencies
find_package(MPI REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()

# ParHIP
option(PARHIP "build ParHIP" OFF)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)

set(LIBKAFFPA_SOURCE_FILES
  lib/data_structure/graph_hierarchy.cpp
  lib/algorithms/strongly_connected_components.cpp
  lib/algorithms/topological_sort.cpp
  lib/algorithms/push_relabel.cpp
  lib/io/graph_io.cpp
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp
  lib/tools/graph_extractor.cpp
  lib/tools/misc.cpp
  lib/tools/partition_snapshooter.cpp
  lib/partition/graph_partitioner.cpp
  lib/partition/w_cycles/wcycle_partitioner.cpp
  lib/partition/coarsening/coarsening.cpp
  lib/partition/coarsening/contraction.cpp
  lib/partition/coarsening/edge_rating/edge_ratings.cpp
  lib/partition/coarsening/matching/matching.cpp
  lib/partition/coarsening/matching/random_matching.cpp
  lib/partition/coarsening/matching/gpa/path.cpp
  lib/partition/coarsening/matching/gpa/gpa_matching.cpp
  lib/partition/coarsening/matching/gpa/path_set.cpp
  lib/partition/coarsening/clustering/node_ordering.cpp
  lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp
  lib/partition/initial_partitioning/initial_partitioning.cpp
  lib/partition/initial_partitioning/initial_partitioner.cpp
  lib/partition/initial_partitioning/initial_partition_bipartition.cpp
  lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp
  lib/partition/initial_partitioning/bipartition.cpp
  lib/partition/initial_partitioning/initial_node_separator.cpp
  lib/partition/uncoarsening/uncoarsening.cpp
  lib/partition/uncoarsening/separator/area_bfs.cpp
  lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp
  lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp
  lib/partition/uncoarsening/refinement/mixed_refinement.cpp
  lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp
  lib/partition/uncoarsening/refinement/refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp
  lib/partition/uncoarsening/refinement/node_separators/greedy_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/fm_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/localized_fm_ns_local_search.cpp
  lib/algorithms/cycle_search.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp
  lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp
  extern/argtable3-3.0.3/argtable3.c)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})

set(LIBKAFFPA_PARALLEL_SOURCE_FILES
  lib/parallel_mh/parallel_mh_async.cpp
  lib/parallel_mh/population.cpp
  lib/parallel_mh/galinier_combine/gal_combine.cpp
  lib/parallel_mh/galinier_combine/construct_partition.cpp
  lib/parallel_mh/exchange/exchanger.cpp
  lib/tools/graph_communication.cpp
  lib/tools/mpi_tools.cpp)
add_library(libkaffpa_parallel OBJECT ${LIBKAFFPA_PARALLEL_SOURCE_FILES})
target_include_directories(libkaffpa_parallel PUBLIC ${MPI_CXX_INCLUDE_PATH})

set(LIBMAPPING_SOURCE_FILES
  lib/mapping/local_search_mapping.cpp
  lib/mapping/full_search_space.cpp
  lib/mapping/full_search_space_pruned.cpp
  lib/mapping/communication_graph_search_space.cpp
  lib/mapping/fast_construct_mapping.cpp
  lib/mapping/construct_distance_matrix.cpp
  lib/mapping/mapping_algorithms.cpp
  lib/mapping/construct_mapping.cpp)
add_library(libmapping OBJECT ${LIBMAPPING_SOURCE_FILES})

set(LIBSPAC_SOURCE_FILES OBJECT lib/spac/spac.cpp)
add_library(libspac ${LIBSPAC_SOURCE_FILES})



add_executable(strongmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(strongmap PRIVATE "-DMODE_STRONGMAP")
target_link_libraries(strongmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS strongmap DESTINATION bin)

add_executable(ecomap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(ecomap PRIVATE "-DMODE_ECOMAP")
target_link_libraries(ecomap ${OpenMP_CXX_LIBRARIES})
install(TARGETS ecomap DESTINATION bin)

add_executable(fastmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(fastmap PRIVATE "-DMODE_FASTMAP")
target_link_libraries(fastmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS fastmap DESTINATION bin)

add_executable(fastestmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(fastestmap PRIVATE "-DMODE_FASTESTMAP")
target_link_libraries(fastestmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS fastestmap DESTINATION bin)



add_executable(ssocialmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(ssocialmap PRIVATE "-DMODE_SSOCIALMAP")
target_link_libraries(ssocialmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS ssocialmap DESTINATION bin)

add_executable(esocialmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(esocialmap PRIVATE "-DMODE_ESOCIALMAP")
target_link_libraries(esocialmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS esocialmap DESTINATION bin)

add_executable(fsocialmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(fsocialmap PRIVATE "-DMODE_FSOCIALMAP")
target_link_libraries(fsocialmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS fsocialmap DESTINATION bin)

add_executable(ffsocialmap app/integratedmapping.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(ffsocialmap PRIVATE "-DMODE_FFSOCIALMAP")
target_link_libraries(ffsocialmap ${OpenMP_CXX_LIBRARIES})
install(TARGETS ffsocialmap DESTINATION bin)



