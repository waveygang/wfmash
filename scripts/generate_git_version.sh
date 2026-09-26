INC_DIR=$1
WFLIGN_DIR=$2

# The VERSION file is the authoritative base version (e.g. 0.24.2);
# it lives at the repository root, one level above scripts/.
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
VERSION_FILE="$SCRIPT_DIR/../VERSION"
BASE_VERSION=$(cat "$VERSION_FILE" 2>/dev/null)

# git describe gives the tag-relative state, e.g. v0.24.2-4-g<sha> or a
# bare sha when the tree has no tags. It needs git on the PATH and a git
# checkout; tarball builds and minimal environments will not have it.
GIT_DESCRIBE=$(git -C "$SCRIPT_DIR/.." describe --always --tags --long 2>/dev/null)

if [ -n "$BASE_VERSION" ]; then
    if [ -n "$GIT_DESCRIBE" ]; then
        # e.g. 0.24.2-4-g6f1747c when the tree is ahead of the VERSION tag
        WFMASH_VERSION="${GIT_DESCRIBE#v}"
    else
        WFMASH_VERSION="$BASE_VERSION"
    fi
else
    WFMASH_VERSION="${GIT_DESCRIBE:-unknown}"
fi

# Write main wfmash version header
echo "#define WFMASH_GIT_VERSION" \"$WFMASH_VERSION\" > "$INC_DIR"/wfmash_git_version.hpp

# Write wflign version header
echo "#define WFLIGN_GIT_VERSION" \"$WFMASH_VERSION\" > "$WFLIGN_DIR"/wflign_git_version.hpp
