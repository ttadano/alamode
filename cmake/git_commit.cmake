# Git commit of the source tree, shown in the program banners and stored in the HDF5 provenance attributes.
#
# include()d from a CMakeLists.txt, this file defines alamode_git_commit(<targets...>): on every build it
# regenerates alamode_git_commit.h in the current binary directory and puts that directory on the targets'
# include path (consumed by git_commit.cpp). Run with -P by that build step, it writes the header, rewriting
# it only when the commit or the dirty state changes so that an unchanged tree recompiles nothing.

if (CMAKE_SCRIPT_MODE_FILE)
    get_filename_component(source_root "${CMAKE_SCRIPT_MODE_FILE}/../.." ABSOLUTE)
    set(commit "unknown")
    find_package(Git QUIET)
    # Require the source root itself to be a work tree, not e.g. a release tarball unpacked inside another repository.
    if (GIT_FOUND AND EXISTS "${source_root}/.git")
        execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
                        WORKING_DIRECTORY "${source_root}"
                        OUTPUT_VARIABLE head
                        RESULT_VARIABLE status
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        ERROR_QUIET)
        if (status EQUAL 0)
            set(commit "${head}")
            # Exit code 1 = tracked files differ from HEAD; untracked files are ignored.
            execute_process(COMMAND "${GIT_EXECUTABLE}" diff --quiet HEAD --
                            WORKING_DIRECTORY "${source_root}"
                            RESULT_VARIABLE dirty
                            ERROR_QUIET)
            if (dirty EQUAL 1)
                string(APPEND commit "-dirty")
            endif ()
        endif ()
    endif ()

    set(content "#define ALAMODE_GIT_COMMIT_ID \"${commit}\"\n")
    set(previous "")
    if (EXISTS "${OUTPUT}")
        file(READ "${OUTPUT}" previous)
    endif ()
    if (NOT content STREQUAL previous)
        file(WRITE "${OUTPUT}" "${content}")
    endif ()
else ()
    function(alamode_git_commit)
        set(header "${CMAKE_CURRENT_BINARY_DIR}/alamode_git_commit.h")
        if (NOT TARGET ${PROJECT_NAME}_git_commit)
            add_custom_target(${PROJECT_NAME}_git_commit
                              COMMAND "${CMAKE_COMMAND}" "-DOUTPUT=${header}" -P "${CMAKE_CURRENT_FUNCTION_LIST_FILE}"
                              BYPRODUCTS "${header}"
                              VERBATIM)
        endif ()
        foreach (target IN LISTS ARGN)
            add_dependencies(${target} ${PROJECT_NAME}_git_commit)
            target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
        endforeach ()
    endfunction ()
endif ()
