{ lib, stdenv, fetchFromGitHub, cmake, gsl, gmp, makeWrapper, jemalloc, htslib, git, zlib, pkgconfig }:

stdenv.mkDerivation rec {
  pname = "wfmash";
  version = "0.12.0";

  src = fetchFromGitHub {
    owner = "waveygang"; 
    repo = "wfmash";
    rev = "01f812e5100d302b593bb3f2687635d987a63edd";
    sha256 = "sha256-X74Bk4+Am/ubhvAD/oH791Dw2QHA4ph2CGsx6sVBu3M=";
  };

  nativeBuildInputs = [ cmake makeWrapper ];

  buildInputs = [ 
    gsl 
    gmp
    jemalloc
    htslib
    git
    zlib
    pkgconfig
  ];

  postPatch = ''
    mkdir -p include
    echo "#define WFMASH_GIT_VERSION \"${version}\"" > include/wfmash_git_version.hpp
  '';

  postInstall = ''
    wrapProgram $out/bin/wfmash --prefix PATH : ${lib.makeBinPath [ gsl gmp ]}
  '';

  meta = with lib; {
    description = "Base-accurate DNA sequence alignments using WFA and mashmap2";
    homepage = "https://github.com/ekg/wfmash";
    license = licenses.mit;
    platforms = platforms.linux;
    maintainers = [ maintainers.bzizou ];
  };
}
