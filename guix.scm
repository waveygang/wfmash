;; To use this file to build a version of wfmash using git HEAD:
;;
;;   guix build -f guix.scm                  # default build
;;   guix build -L . wfmash-gcc-git          # stanard gcc build
;;   guix build -L . wfmash-gcc-debug-git    # gcc build with debug and ASAN
;;   guix build -L . wfmash-gcc-profile-git  # run the profiler!
;;   guix build -L . wfmash-static-gcc-git   # gcc static build (default)
;;   guix build -L . wfmash-clang-git        # clang build
;;
;; To get a development container using a recent guix (see `guix pull`)
;;
;;   guix shell -C -D -F -f guix.scm         # default build
;;   guix shell -L . -C -D -F wfmash-gcc-git # preferred development container
;;   guix shell -L . -C -D -F wfmash-gcc-static-git
;;   guix shell -L . -C -D -F wfmash-clang-git
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
;; alternative builds
;;
;;   cmake -DCMAKE_BUILD_TYPE=Debug ..           # for development (use wfmash-gcc-git)
;;   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..  # for distros including Debian (use wfmash-gcc-git)
;;   cmake -DBUILD_STATIC=1 ..                   # static binary (use wfmash-gcc-static-git)
;;
;; list packages
;;
;;   guix package -L . -A|grep wfm
;;
;; Installing guix (note that Debian comes with guix). Once installed update as a normal user with:
;;
;;   mkdir ~/opt
;;   guix pull -p ~/opt/guix # update guix takes a while - don't do this often!
;;
;; Use the update guix to build wfmash:
;;
;;   ~/opt/guix/bin/guix build -f guix.scm
;;
;; Or get a shell
;;
;;   ~/opt/guix/gin/guix build -f guix.scm
;;
;; If things do not work you may also have to update the guix-daemon in systemd. Guix mostly downloads binary
;; substitutes. If it wants to build a lot of software you probably have substitutes misconfigured.

;; by Pjotr Prins & Andrea Guarracino (c) 2023-2024

(define-module (guix)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system cargo)
  #:use-module (guix build-system cmake)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:use-module (guix utils)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bash)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages build-tools)
  #:use-module (gnu packages certs)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages cpp)
  #:use-module (gnu packages crates-io)
  #:use-module (gnu packages curl)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages jemalloc)
  #:use-module (gnu packages linux) ; for util-linux column
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages perl)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages python)
  #:use-module (gnu packages rust)
  #:use-module (gnu packages rust-apps) ; for cargo
  #:use-module (gnu packages tls)
  #:use-module (gnu packages version-control)
  #:use-module (srfi srfi-1)
  #:use-module (ice-9 popen)
  #:use-module (ice-9 rdelim)
  )

(define-public pafcheck-github
  (package
    (name "pafcheck-github")
    (version "0.1.0")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/ekg/pafcheck")
             (commit (string-append "v" version))))
       (file-name (git-file-name name version))
       (sha256
        (base32
         "0prlkq8am3sskg55x7b8vr4j54dmkjqldyl50isq5qyy9pff3xxs"))))
    (build-system cargo-build-system)
    (inputs (list curl gnutls lzip openssl pkg-config zlib xz)) ;; mostly for htslib
    (arguments
     `(#:cargo-inputs (("rust-addr" ,rust-addr-0.14)
                       ("rust-anyhow" ,rust-anyhow-1)
                       ("rust-cssparser" ,rust-cssparser-0.28)
                       ("rust-clap" ,rust-clap-4)
                       ("rust-lifeguard" ,rust-lifeguard-0.6)
                       ("rust-rmp-serde" ,rust-rmp-serde-1)
                       ("rust-rust-htslib" ,rust-rust-htslib-0.38)
                       ("rust-tempfile" ,rust-tempfile-3)
                       ("rust-thiserror" ,rust-thiserror-1))
       ;; #:cargo-development-inputs ()))
       #:cargo-package-flags '("--no-metadata" "--no-verify" "--allow-dirty")
     ))
    (synopsis "pafcheck")
    (description
     "Tool for validating PAF (Pairwise Alignment Format) files against their corresponding FASTA sequences. It ensures that the alignments described in the PAF file match the actual sequences in the FASTA files.")
    (home-page "https://github.com/ekg/pafcheck")
    (license license:expat)))

(define-public pafcheck-shell-git
  "Shell version to use 'cargo build'"
  (package
    (inherit pafcheck-github)
    (name "pafcheck-shell-git")
    ;; (version (git-version "0.21" "HEAD" %git-commit))
    (inputs
     (modify-inputs (package-inputs pafcheck-github)
         (append binutils coreutils-minimal ;; for the shell
                 )))
    (propagated-inputs (list cmake rust rust-cargo nss-certs openssl perl gnu-make-4.2
                             coreutils-minimal which perl binutils gcc-toolchain pkg-config zlib
                             )) ;; to run cargo build in the shell
    (arguments
     `(
       #:cargo-inputs (("rust-anyhow" ,rust-anyhow-1)
                       ("rust-clap" ,rust-clap-4)
                       ("rust-rust-htslib" ,rust-rust-htslib-0.38)
                       ("rust-tempfile" ,rust-tempfile-3)
                       ("rust-thiserror" ,rust-thiserror-1)
                       )
       ;; #:cargo-development-inputs ()))
       #:cargo-package-flags '("--no-metadata" "--no-verify" "--allow-dirty")
       #:phases (modify-phases %standard-phases
                               (delete 'configure)
                               (delete 'build)
                               (delete 'package)
                               (delete 'check)
                               (delete 'install)
                               )
     ))
    ))


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
       ("pafcheck" ,pafcheck-shell-git)
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
  "Default build with gcc - as is used in distros"
  (package
    (inherit wfmash-base-git)
    (name "wfmash-gcc-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append gcc-14
                 )))
    ))

(define-public wfmash-gcc-debug-git
  "Build with debug options"
  (package
    (inherit wfmash-gcc-git)
    (name "wfmash-gcc-debug-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (arguments
     `(;; #:tests? #f ;; skip tests, this is mostly to run a shell
       #:configure-flags
       ,#~(list
           "-DASAN=ON"
           "-DDISABLE_LTO=ON"
           "-DCMAKE_BUILD_TYPE=Debug"))) ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
         (append gperftools
                 )))
    (arguments
     `(#:phases (modify-phases %standard-phases
                               (delete 'configure)
                               (delete 'build)
                               (delete 'package)
                               (delete 'check)
                               (delete 'install)
                               )
     ))
    ))

(define-public wfmash-clang-git
  "Clang+LLVM build"
  (package
    (inherit wfmash-base-git)
    (name "wfmash-clang-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append clang-toolchain-18
                 lld
                 libomp
                 )))
    ))

(define-public wfmash-gcc-profile-git
  "Build wfmash optimally and automatically run profiler on all tests"
  (package
    (inherit wfmash-gcc-git)
    (name "wfmash-gcc-profile-git")
    (version (git-version "0.21" "HEAD" %git-commit))
    (arguments
     `(#:tests? #f ;; running tests as profiler
       #:configure-flags
         ,#~(list
             "-DCMAKE_BUILD_TYPE=Release"
             "-DBUILD_OPTIMIZED=ON"
             "-DPROFILER=ON")
       #:phases
         ,#~(modify-phases %standard-phases
            (add-after 'install 'run-profiler-on-all2all
                       (lambda* (#:key outputs #:allow-other-keys)
                         (invoke "ctest" "--verbose" "-R" "all2all") ; only run all2all test
                         (invoke "ls" "-l" "bin/wfmash")
                         (invoke "ls" "-l")
                         (invoke "pprof" "--text" "bin/wfmash" "wfmash.prof")
                         (mkdir-p (string-append #$output:doc "/share")))))))
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
                    (append gperftools
                            coreutils ;; for ls
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
  "Optimized for latest AMD architecture build and static deployment.
These binaries can be copied to HPC."
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

wfmash-static-gcc-git ;; default optimized static deployment build
