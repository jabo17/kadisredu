file(GLOB_RECURSE KADISREDU_SOURCE_FILES CONFIGURE_DEPENDS
        *.cpp *.h)

message(${KADISREDU_SOURCE_FILES})

add_library(kadisredu ${KADISREDU_SOURCE_FILES})
target_include_directories(kadisredu PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")

# we need to include kamis first, because header names are not fully qualified
target_link_libraries(kadisredu PRIVATE kamis::wmis)
target_include_directories(kadisredu PRIVATE ${KAMIS_WMIS_KAHIP_INCLUDES})
target_link_libraries(kadisredu PUBLIC KaMinPar::dKaMinPar kamping KaGen::KaGen message-queue::message-queue)
target_link_libraries(kadisredu PUBLIC absl::flat_hash_map)
target_link_libraries(kadisredu PUBLIC fmt::fmt)
target_link_libraries(kadisredu PUBLIC nlohmann_json::nlohmann_json)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(kadisredu PRIVATE -g -O0)
endif ()



