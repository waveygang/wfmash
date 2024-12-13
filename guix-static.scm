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
;;   cmake -DBUILD_STATIC=1 ..
;;   make
;;
;; For the tests you may need /usr/bin/env. Inside the container:
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;
;; by Pjotr Prins (c) 2023-2024

(use-modules
 ((guix licenses) #:prefix license:)
  (guix build-system cmake)
  (guix download)
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

;; Guix does not come with a static version of libdeflate
(define-public libdeflate-static
  (package
    (inherit libdeflate)
    (name "libdeflate-static")
    (version "1.19")
    (arguments
     (list #:configure-flags
           #~(list "-DLIBDEFLATE_BUILD_STATIC_LIB=YES"
                   "-DLIBDEFLATE_BUILD_TESTS=YES")))))

;; A minimal static version of htslib that does not depend on curl and openssl. This
;; reduces the number of higher order dependencies in static linking.
(define-public htslib-static
  (package
    (inherit htslib)
    (name "htslib-static")
    (version "1.19")
    (source (origin
            (method url-fetch)
            (uri (string-append
                  "https://github.com/samtools/htslib/releases/download/"
                  version "/htslib-" version ".tar.bz2"))
            (sha256
             (base32
              "0dh79lwpspwwfbkmllrrhbk8nkvlfc5b5ib4d0xg5ld79w6c8lc7"))))
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
    (name "wfmash-static-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (arguments
     `(;; #:tests? #f
       #:configure-flags
       ,#~(list
           "-DBUILD_STATIC=ON"
           "-DBUILD_OPTIMIZED=ON"
           "-DCMAKE_BUILD_TYPE=Release"
           "-DCMAKE_INSTALL_RPATH="))) ; force cmake static build and do not rewrite RPATH
    (inputs
     `(
       ("bzip2-static" ,bzip2 "static")    ; libz2 part of htslib for static builds
       ("coreutils" ,coreutils) ; for echo and env in tests
       ("gcc" ,gcc-12)
       ("git" ,git)
       ("gmp" ,gmp)
       ("gsl-static" ,gsl "static")
       ("gsl" ,gsl)
       ("htslib-static" ,htslib-static)
       ("jemalloc" ,jemalloc)
       ("libdeflate-static" ,libdeflate-static)
       ("make" ,gnu-make)
       ("pkg-config" ,pkg-config)
       ("xz-static" ,xz "static")     ; for static builds
       ("zlib-static" ,zlib "static")))
     (synopsis "wfmash")
     (description
      "wfmash is an aligner for pangenomes that combines efficient homology
mapping with base-level alignment. It uses MashMap to find approximate
mappings between sequences, then applies WFA (Wave Front Alignment) to
obtain base-level alignments.")
     (home-page "https://github.com/waveygang/wfmash")
     (license license:expat)))

wfmash-git
