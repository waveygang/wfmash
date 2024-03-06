{ pkgs ? import <nixpkgs> {} }:

pkgs.callPackage ./wfmash.nix {
  inherit (pkgs) stdenv fetchFromGitHub cmake gsl gmp makeWrapper jemalloc htslib git zlib pkg-config;
}