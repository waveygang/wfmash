;; To use this file to build a version of wfmash using git HEAD:
;;
;;  rm -rf build/*
;;  guix build -f guix.scm
;;
;; To get a development container using a recent guix (see `guix pull`)
;;
;;   guix shell -C -D -f guix.scm
;;
;; and inside the container
;;
;;   mkdir build
;;   cd build
;;   cmake ..
;;   make -j 12
;;
;; For the tests you may need /usr/bin/env. Inside the container:
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;
;; by Pjotr Prins & Andrea Guarracino (c) 2023-2024

(use-modules
 ((guix licenses) #:prefix license:)
 (guix build-system cmake)
 (guix gexp)
 (guix git-download)
 (guix packages)
 (guix utils)
 (gnu packages algebra)
 (gnu packages base)
 (gnu packages bioinformatics)
 (gnu packages build-tools)
 (gnu packages compression)
 (gnu packages gcc)
 (gnu packages jemalloc)
 (gnu packages llvm)
 (gnu packages maths)
 (gnu packages multiprecision)
 (gnu packages pkg-config)
 (gnu packages python)
 (gnu packages version-control)
 (srfi srfi-1)
 (ice-9 popen)
 (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public wfmash-git
  (package
    (name "wfmash-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (arguments
     `(#:tests? #f)) ; disable tests until I fix finding the binary wfmash
    (inputs
     `(
       ("bzip2" ,bzip2)
       ("coreutils" ,coreutils) ; for echo and env in tests
       ("gcc" ,gcc-12)
       ("git" ,git)
       ("gmp" ,gmp)
       ("gsl" ,gsl)
       ("htslib" ,htslib)
       ("libdeflate" ,libdeflate)
       ("make" ,gnu-make)
       ("pkg-config" ,pkg-config)
       ("xz" ,xz)
       ("zlib" ,zlib)))
     (synopsis "wfmash")
     (description
      "wfmash is an aligner for pangenomes that combines efficient homology
mapping with base-level alignment. It uses MashMap to find approximate
mappings between sequences, then applies WFA (Wave Front Alignment) to
obtain base-level alignments.")
     (home-page "https://github.com/waveygang/wfmash")
     (license license:expat)))

wfmash-git
