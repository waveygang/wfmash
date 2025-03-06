#!/bin/bash

# run from root of the repository

mkdir source-tarball
cd source-tarball
git clone --recursive https://github.com/waveygang/wfmash
cd wfmash
git fetch --tags origin
LATEST_TAG="$(git describe --tags `git rev-list --tags --max-count=1`)"
git checkout "${LATEST_TAG}"
git submodule update --init --recursive
bash scripts/generate_git_version.sh $PWD/src $PWD/src/common/wflign/src
sed 's/execute_process(COMMAND bash/#execute_process(COMMAND bash/g' CMakeLists.txt -i
rm -Rf .git
find src/common -name ".git" -exec rm -Rf "{}" \;
cd ..
mv wfmash "wfmash-${LATEST_TAG}"
tar -czf "wfmash-${LATEST_TAG}.tar.gz" "wfmash-${LATEST_TAG}"
rm -Rf "wfmash-${LATEST_TAG}"
