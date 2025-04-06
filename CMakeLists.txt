find_package(FFTW REQUIRED)
find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)

include_directories(${FFTW_INCLUDES})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

# ------------ Compile library ------------
add_definitions(-DCOMPUTE_STATISTICS)


set(SRC 
    ./st_dsp.c
    ./st_io.c
    )

# Force c++
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )

add_library(st SHARED ${SRC})


# ------------ Compile tone_track ------------

set(SRC 
    ./tone_track.c
    )

# Force c++
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )

add_executable(tone_track ${SRC})
target_link_libraries(tone_track st)
target_link_libraries(tone_track ${FFTW_LIBRARIES})
target_link_libraries(tone_track ${ARMADILLO_LIBRARIES})
target_link_libraries(tone_track ${GSL_LIBRARIES})
