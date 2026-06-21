#!/usr/bin/env bash
# Enter or source the Guix environment used to build wfmash from this checkout.
#
# Start an initialized shell:
#   ./env.sh
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDISABLE_LTO=ON
#   cmake --build build --parallel
#
# Or add the Guix profile to the current shell:
#   source ./env.sh
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDISABLE_LTO=ON
#   cmake --build build --parallel
#
# Static project/vendor-library build:
#   source ./env.sh
#   cmake -S . -B build-static -DCMAKE_BUILD_TYPE=Release -DVENDOR_EVERYTHING=ON -DDISABLE_LTO=ON
#   cmake --build build-static --parallel
#
set -euo pipefail

_wfmash_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_guix="${GUIX:-${HOME}/.guix-profile/bin/guix}"
_manifest="$_wfmash_dir/manifest.scm"

if [[ ! -x "$_guix" ]]; then
    echo "error: Guix command not found or not executable: $_guix" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    exec bash --noprofile --norc -c 'source "$1"; exec bash --noprofile --norc -i' bash "$_wfmash_dir/env.sh"
fi

set +u
_search_paths="$("$_guix" shell -m "$_manifest" --search-paths)"
if [[ -z "$_search_paths" ]]; then
    echo "error: Guix did not return a usable environment" >&2
    return 1 2>/dev/null || exit 1
fi
eval "$_search_paths"
set -u

_profile_bin="${PATH%%:*}"
_profile="${_profile_bin%/bin}"

export GUIX_PROFILE="$_profile"
export CC=gcc
export CXX=g++

_rpath_flag="-Wl,-rpath,$_profile/lib"
case " ${LDFLAGS-} " in
    *" $_rpath_flag "*) ;;
    *) export LDFLAGS="$_rpath_flag${LDFLAGS:+ $LDFLAGS}" ;;
esac

_cmake_version="$(cmake --version | { read -r _ _ _version; printf '%s' "$_version"; })"
echo "wfmash Guix build env ready: gcc $(gcc -dumpfullversion -dumpversion), cmake $_cmake_version"

unset _wfmash_dir _guix _manifest _search_paths _profile_bin _profile _rpath_flag _cmake_version _version
