{
  description = "A flake for building wfmash and a Docker image for it";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }: let
    system = "x86_64-linux";
    pkgs = import nixpkgs { inherit system; };
  in {
    packages.${system}.wfmash = pkgs.stdenv.mkDerivation rec {
      pname = "wfmash";
      version = "0.14.0";

      src = pkgs.fetchFromGitHub {
        owner = "waveygang";
        repo = "wfmash";
        rev = "7376468d0d1f67ad58ca3fc5d07e888745b22c06";
        sha256 = "sha256-0vThKJflmasPoVhGdFwm1l9VUoxZkEwZruYhlzF1Ehw=";
      };

      nativeBuildInputs = [ pkgs.cmake pkgs.makeWrapper ];

      buildInputs = [
        pkgs.gsl
        pkgs.gmp
        pkgs.jemalloc
        pkgs.htslib
        pkgs.git
        pkgs.zlib
        pkgs.pkg-config
      ];

      # Define custom attributes
      enableOptimizations = true;
      reproducibleBuild = false;

      # Use custom attributes to set compiler flags
      CFLAGS = if enableOptimizations then "-Ofast -march=x86-64-v3" else "";
      CXXFLAGS = if enableOptimizations then "-Ofast -march=x86-64-v3" else "";

      postPatch = ''
        mkdir -p include
        echo "#define WFMASH_GIT_VERSION \"${version}\"" > include/wfmash_git_version.hpp
      '';

      postInstall = ''
        wrapProgram $out/bin/wfmash --prefix PATH : ${pkgs.lib.makeBinPath [ pkgs.gsl pkgs.gmp ]}
      '';

      meta = with pkgs.lib; {
        description = "Base-accurate DNA sequence alignments using WFA and mashmap2";
        homepage = "https://github.com/ekg/wfmash";
        license = licenses.mit;
        platforms = platforms.linux;
        maintainers = [ maintainers.bzizou ];
      };
    };

    dockerImage = pkgs.dockerTools.buildImage {
      name = "wfmash-docker";
      tag = "latest";
      copyToRoot = [ self.packages.${system}.wfmash ];
      config = {
        Entrypoint = [ "${self.packages.${system}.wfmash}/bin/wfmash" ];
      };
    };
  };
}
