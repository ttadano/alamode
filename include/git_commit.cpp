/*
 git_commit.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

// Shared by alm and anphon. The commit lives in its own translation unit so that a new commit recompiles
// only this file. CMake builds generate alamode_git_commit.h at build time (cmake/git_commit.cmake);
// other builds report "unknown".

#include "version.h"

#if __has_include("alamode_git_commit.h")
#include "alamode_git_commit.h"
#else
#define ALAMODE_GIT_COMMIT_ID "unknown"
#endif

const char *const ALAMODE_GIT_COMMIT = ALAMODE_GIT_COMMIT_ID;
