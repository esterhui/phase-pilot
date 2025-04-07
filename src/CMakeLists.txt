find_package(FFTW3 REQUIRED)

include_directories(${FFTW_INCLUDES})

# ------------ Compile library ------------
#add_definitions(-DCOMPUTE_STATISTICS)
#add_definitions(-DFITTER)


set(SRC
    st_dsp.c
    st_io.c
    )

# Force c++
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )

add_library(st SHARED ${SRC})
target_link_libraries(st FFTW3::fftw3)


# ------------ Compile tone_track ------------

set(SRC tone_track.c )
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )
add_executable(tone_track ${SRC})
target_link_libraries(tone_track st )

set(SRC prn_track.c )
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )
add_executable(prn_track ${SRC})
target_link_libraries(prn_track st )

set(SRC fft_acquire.c )
SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE CXX )
add_executable(fft_acquire ${SRC})
target_link_libraries(fft_acquire st )
