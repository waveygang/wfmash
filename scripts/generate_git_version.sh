INC_DIR=$1
WFLIGN_DIR=$2

# Go to the directory where the script is
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"

GIT_VERSION=$(git describe --always --tags --long)

# Write main wfmash version header
echo "#define WFMASH_GIT_VERSION" \"$GIT_VERSION\" > "$INC_DIR"/wfmash_git_version.hpp

# Write wflign version header 
echo "#define WFLIGN_GIT_VERSION" \"$GIT_VERSION\" > "$WFLIGN_DIR"/wflign_git_version.hpp
