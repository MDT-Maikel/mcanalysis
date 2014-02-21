##############################
### CMake file for scripts ###
##############################

## Add custum target to trigger copying of scripts
add_custom_target(copy ALL COMMENT "Copying scripts")

## Copy the scripts to the build directory
add_custom_command(TARGET copy COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/scripts ${PROJECT_BINARY_DIR}/scripts)


