;; To use this file to build a version of wfmash using git HEAD:
;;
;;  rm -rf build/*
;;  guix build -f guix.scm
;;
;; To get a development container using a recent guix (see `guix pull`)
;;
;;   guix shell -C -D -F -f guix.scm
;;
;; and inside the container
;;
;;   mkdir build
;;   cd build
;;   cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_OPTIMIZED=1 ..
;;   make -j 12 VERBOSE=1
;;   ctest . --verbose
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
 (gnu packages bash)
 (gnu packages bioinformatics)
 (gnu packages build-tools)
 (gnu packages compression)
 (gnu packages gcc)
 (gnu packages jemalloc)
 (gnu packages linux) ; for util-linux column
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

(define-public wfmash-base-git
  (package
    (name "wfmash-base-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (inputs
     `(
       ("bash" ,bash) ; for testing
       ("bedtools" ,bedtools) ; for testing
       ("util-linux" ,util-linux) ; for testing
       ("samtools" ,samtools) ; for testing
       ("bzip2" ,bzip2)
       ("coreutils" ,coreutils) ; for echo and env in tests
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

(define-public wfmash-gcc-git
  (package
    (inherit wfmash-base-git)
    (name "wfmash-gcc-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append gcc-14
                 )))
    ))

(define-public wfmash-clang-git
  (package
    (inherit wfmash-base-git)
    (name "wfmash-clang-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append clang-toolchain-17
                 lld
                 libomp
                 )))
    ))


wfmash-gcc-git
;; wfmash-clang-git
