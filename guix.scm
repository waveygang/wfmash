;; To use this file to build a version of wfmash using git HEAD:
;;
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
;;   rm -rf ../build/*
;;   cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_OPTIMIZED=1 ..
;;   make -j 12 VERBOSE=1
;;   ctest . --verbose
;;
;; alternative builds are
;;
;;   cmake -DBUILD_STATIC=1 ..
;;
;; by Pjotr Prins & Andrea Guarracino (c) 2023-2024

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
       ("jemalloc" ,jemalloc)
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

;; ==== The following is for static binary builds using gcc - used mostly for deployment ===

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

(define-public wfmash-static-gcc-git
  (package
    (inherit wfmash-base-git)
    (name "wfmash-static-gcc-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (arguments
     `(;; #:tests? #f
       #:configure-flags
       ,#~(list
           "-DBUILD_STATIC=ON"
           "-DBUILD_OPTIMIZED=ON"
           "-DCMAKE_BUILD_TYPE=Release"
           "-DCMAKE_INSTALL_RPATH="))) ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
                    (prepend
                     `(,bzip2 "static")
                     `(,zlib "static")
                     `(,gsl "static")
                     `(,xz "static")
                     libdeflate-static
                     htslib-static
                     )))))

wfmash-static-gcc-git ;; default deployment build
;; wfmash-gcc-git
;; wfmash-clang-git
