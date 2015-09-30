set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
add_executable (SM Data.cpp MyDistribution.cpp MyModel.cpp main.cpp)
set_target_properties (SM PROPERTIES OUTPUT_NAME main)
target_link_libraries (SM rjobject ${RJOBJ_DEPS})
