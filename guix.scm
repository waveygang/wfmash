;; To use this file to build a static version of wfmash using git HEAD:
;;
;;   guix build -f guix.scm
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
;;   make
;;
;; For the tests you may need /usr/bin/env. Inside the container:
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;
;; by Pjotr Prins (c) 2023

(use-modules
 ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system cmake)
  ; (guix gexp)
  (guix utils)
  (gnu packages algebra)
  (gnu packages base)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages compression)
  ; (gnu packages curl)
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

;; A minimal version of htslib that does not depend on curl and openssl. This
;; reduces the number of higher order dependencies in static linking.
(define-public htslib-static
  (package
    (inherit htslib)
    (name "htslib-static")
    (arguments
     (substitute-keyword-arguments (package-arguments htslib)
       ((#:configure-flags flags ''())
        ''())))
    (inputs
     (list bzip2 xz))))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public wfmash-git
  (package
    (name "wfmash-git")
    (version (git-version "0.10.7" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (inputs
     `(
       ;; ("clang" ,clang)      ; add this to test clang builds
       ;; ("lld" ,lld)          ; add this to test clang builds
       ("gcc" ,gcc-12)
       ("gsl" ,gsl)
       ("gmp" ,gmp)
       ("make" ,gnu-make)
       ("pkg-config" ,pkg-config)
       ("jemalloc" ,jemalloc)
       ("htslib" ,htslib)
       ("git" ,git)
       ; ("bc" ,bc)               ; for tests
       ("coreutils" ,coreutils) ; for echo and env in tests
       ; ("curl" ,curl)
       ; ("parallel" ,parallel) ; for wfmash-parallel
       ("bzip2" ,bzip2)    ; libz2 part of htslib
       ("xz" ,xz)     ;
       ("zlib" ,zlib)))
     (synopsis "wfmash")
     (description
      "wfmash.")
     (home-page "https://github.com/waveygang/wfmash")
     (license license:expat)))

wfmash-git