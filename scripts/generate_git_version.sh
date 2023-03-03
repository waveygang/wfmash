INC_DIR=$1

GIT_VERSION=$(git describe --always --tags --long)

echo "#define WFMASH_GIT_VERSION" \"$GIT_VERSION\" > "$INC_DIR"/wfmash_git_version.hpp
