

set(CMAKE_BUILD_TYPE Debug)


set(MPI_SCHEDULER_SOURCES main.cpp
  Checkpoint.cpp
  Command.cpp 
  CommandsRunner.cpp 
  Common.cpp 
  DynamicLibrary.cpp 
  Instance.cpp
  ParallelImplementation.cpp
  RunStatistics.cpp
  SVGDrawer.cpp
  fork_implem/ForkImplem.cpp
)

set (MPI_SCHEDULER_LIBRARIES dl)

if (WITH_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  set(MPI_SCHEDULER_SOURCES ${MPI_SCHEDULER_SOURCES}
    split_implem/SplitImplem.cpp
    onecore_implem/OneCoreImplem.cpp
  )
  set(MPI_SCHEDULER_LIBRARIES ${MPI_LIBRARIES} ${MPI_SCHEDULER_LIBRARIES})
endif()


if (MPISCHEDULER_BUILD_AS_LIBRARY)
  add_library(mpi-scheduler ${MPI_SCHEDULER_SOURCES})
else()
  add_executable(mpi-scheduler ${MPI_SCHEDULER_SOURCES}) 
endif()

if (MPISCHEDULER_BUILD_AS_LIBRARY AND NOT APPLE)
  target_link_libraries(mpi-scheduler gcc_s gcc)
endif()

target_link_libraries(mpi-scheduler ${MPI_SCHEDULER_LIBRARIES})

if (WITH_MPI)
  if(MPI_COMPILE_FLAGS)
    set_target_properties(mpi-scheduler PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(mpi-scheduler PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()


