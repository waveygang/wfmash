{ pkgs ? import <nixpkgs> { }
, pkgsLinux ? import <nixpkgs> { system = "x86_64-linux"; }
}:

let
  wfmash = pkgs.callPackage ./wfmash.nix { };
in
pkgs.dockerTools.buildImage {
  name = "wfmash-docker";
  tag = "latest";
  copyToRoot = [ wfmash ];
  config = {
    Entrypoint = [ "${wfmash}/bin/wfmash" ];
  };
}
