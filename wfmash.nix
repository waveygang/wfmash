{ lib, stdenv, fetchFromGitHub, cmake, gsl, gmp, makeWrapper, jemalloc, htslib, git, zlib, pkgconfig }:

stdenv.mkDerivation rec {
  pname = "wfmash";
  version = "0.10.5";

  src = fetchFromGitHub {
    owner = "waveygang"; 
    repo = "wfmash";
    rev = "3b86528b7f83ff23a74ad14d9302b8f73c246f91";
    sha256 = "0n7j1x3jd9zy25pqsal8qwnpgnqbdq7xnppsbiw8m4zbhgn1gc9g";
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
