;; To use this file to build a version of wfmash using git HEAD:
;;
;;   guix build -f guix.scm                  # default build
;;   guix build -L . wfmash-gcc-git          # standard gcc build
;;   guix build -L . wfmash-gcc-debug-git    # gcc build with debug and ASAN
;;   guix build -L . wfmash-gcc-profile-git  # run the profiler!
;;   guix build -L . wfmash-gcc-static-git --without-tests=wfmash-gcc-static-git # gcc static build (default)
;;   guix build -L . wfmash-clang-git        # clang build
;;
;; Note that the build happens in a fresh isolated container, even though it uses the current git repo.
;; Note that for abov build commands testing may be disabled with guix option --without-tests=target.
;;
;; To get a development container using a recent guix (see `guix pull`)
;;
;;   guix shell --share=$HOME/.cargo -C -D -F -v 3 -f guix.scm            # default build
;;   guix shell --share=$HOME/.cargo -L . -C -D -F wfmash-gcc-git         # preferred development container
;;   guix shell --share=$HOME/.cargo -L . -C -D -F wfmash-gcc-static-git
;;   guix shell --share=$HOME/.cargo -L . -C -D -F wfmash-clang-git
;;
;; To include a prebuilt wgatools binary:
;;
;;   guix shell --share=$HOME/.cargo -C -D -F -v 3 --expose=../wgatools=/wgatools -L . wfmash-gcc-debug-git
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
;; For tests to work build wgatools in the next directory and add --expose=../wgatools=/wgatools to the shell options above.
;; Inside the container you should be able to run:
;;
;; env LD_LIBRARY_PATH=$GUIX_PROFILE/lib /wgatools/target/release/wgatools
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

;; by Pjotr Prins & Andrea Guarracino (c) 2023-2025

(define-module (guix)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system cargo)
  #:use-module (guix build-system cmake)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix gexp)
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
  #:use-module (gnu packages crates-check) ; for cargo
  #:use-module (gnu packages crates-io)
  #:use-module (gnu packages crates-web)
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
  (let ((commit "6aad6af378d91c9ee1ecb7865cf2f8ecbe67cb90"))
    (package
     (name "pafcheck-github")
     (version (string-append "0.1.0-" (string-take commit 7)))
     (source
      (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/pjotrp/pafcheck")
             (commit commit)))
       (file-name (string-append name "-" version "-checkout"))
       (sha256
        (base32
         "1dndbym1a0pz5asa7xs6369dszmxpyxanb0naknlx6b62f6r8j1c"
        ))))
     (build-system cargo-build-system)
     (inputs (list curl gnutls lzip openssl pkg-config zlib xz)) ;; mostly for htslib
     (arguments
      `(#:cargo-inputs (("rust-addr" ,rust-addr-0.14)
                        ("rust-adblock" ,rust-adblock-0.7)
                        ("rust-anyhow" ,rust-anyhow-1)
                        ("rust-cssparser" ,rust-cssparser-0.28)
                        ("rust-clap" ,rust-clap-4)
                        ("rust-criterion" ,rust-criterion-0.4)
                        ("rust-criterion" ,rust-criterion-0.5)
                        ("rust-reqwest" ,rust-reqwest-0.11)
                        ("rust-mock-instant" ,rust-mock-instant-0.2)
                        ("rust-lifeguard" ,rust-lifeguard-0.6)
                        ("rust-rmp-serde" ,rust-rmp-serde-1)
                        ("rust-rust-htslib" ,rust-rust-htslib-0.38)
                        ("rust-tempfile" ,rust-tempfile-3)
                        ("rust-thiserror" ,rust-thiserror-1))
        ;; #:cargo-development-inputs ()))
        #:cargo-package-flags '("--no-metadata" "--no-verify" "--allow-dirty")
        #:tests? #f))
     (synopsis "pafcheck")
     (description
      "Tool for validating PAF (Pairwise Alignment Format) files against their corresponding FASTA sequences. It ensures that the alignments described in the PAF file match the actual sequences in the FASTA files.")
     (home-page "https://github.com/ekg/pafcheck")
     (license license:expat))))

(define-public pafcheck-shell-git
  "Shell version to use 'cargo build'"
  (package
    (inherit pafcheck-github)
    (name "pafcheck-shell-git")
    (inputs
     (modify-inputs (package-inputs pafcheck-github)
         (append binutils coreutils-minimal ;; for the shell
                 )))
    (propagated-inputs (list cmake rust rust-cargo nss-certs openssl perl gnu-make-4.2
                             coreutils-minimal which perl binutils gcc-toolchain pkg-config zlib
                             )) ;; to run cargo build in the shell
    ))


(define %source-dir (dirname (current-filename)))

(define %version
  (read-string (open-pipe "git describe --always --tags --long|tr -d $'\n'" OPEN_READ)))

(define-public wfmash-base-git
  (package
    (name "wfmash-base-git")
    (version %version)
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (properties '((tunable? . #t)))
    (inputs (list
             bash bedtools util-linux samtools ; for testing
             coreutils bzip2 git gmp gsl htslib libdeflate gnu-make pkg-config xz zlib))
    (propagated-inputs (list
                        pafcheck-github))
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
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append gcc)))))

(define-public wfmash-gcc-debug-git
  "Build with debug options"
  (package
    (inherit wfmash-gcc-git)
    (name "wfmash-gcc-debug-git")
    (arguments
     `(;; #:tests? #f ;; skip tests, this is mostly to run a shell
       #:configure-flags
       ,#~(list
           "-DASAN=ON"
           "-DDISABLE_LTO=ON"
           "-DCMAKE_BUILD_TYPE=Debug"))) ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
                    (prepend bzip2)
                    (append gperftools)))
    (propagated-inputs (list
                        pafcheck-github))
    (arguments
     `(#:phases (modify-phases %standard-phases
                               (delete 'configure)
                               (delete 'build)
                               (delete 'package)
                               (delete 'check)
                               (delete 'install))))))

(define-public wfmash-clang-git
  "Clang+LLVM build"
  (package
    (inherit wfmash-base-git)
    (name "wfmash-clang-git")
    (inputs
     (modify-inputs (package-inputs wfmash-base-git)
         (append clang-toolchain-18
                 lld
                 libomp)))))

(define-public wfmash-gcc-profile-git
  "Build wfmash optimally and automatically run profiler on all tests"
  (package
    (inherit wfmash-gcc-git)
    (name "wfmash-gcc-profile-git")
    (arguments
     `(#:tests? #f ;; running tests as profiler
       #:configure-flags
         ,#~(list
             "-DCMAKE_BUILD_TYPE=Generic"
             ;; "-DBUILD_OPTIMIZED=ON" -- use --tune switch
             "-DPROFILER=ON")
       #:phases
         ,#~(modify-phases %standard-phases
            (add-after 'install 'run-profiler
                       (lambda* (#:key outputs #:allow-other-keys)
                         (invoke "ctest" "--verbose" "-R" "wfmash-time-LPA")
                         (invoke "ctest" "--verbose" "-R" "wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF")
                         (invoke "ls" "-l" "bin/wfmash")
                         (invoke "ls" "-l")
                         (invoke "pprof" "--text" "bin/wfmash" "wfmash.prof")
                         (mkdir-p (string-append #$output:doc "/share")))))))
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
                    (append gperftools
                            coreutils ;; for ls
                 )))))

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

(define-public wfmash-gcc-static-git
  "Optimized for latest AMD architecture build and static deployment.
These binaries can be copied to HPC."
  (package
    (inherit wfmash-base-git)
    (name "wfmash-gcc-static-git")
    (arguments
     `(#:tests? #f  ;; no Rust tools
       #:configure-flags
       ,#~(list
           "-DBUILD_STATIC=ON"
           ;; "-DBUILD_OPTIMIZED=ON"    ;; we don't use the standard cmake optimizations
           "-DCMAKE_BUILD_TYPE=Generic" ;; to optimize use guix --tune=march-type (e.g. --tune=native)
           "-DCMAKE_INSTALL_RPATH=")))   ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs wfmash-gcc-git)
                    (delete pafcheck-github)
                    (prepend
                     `(,bzip2 "static")
                     `(,zlib "static")
                     `(,gsl "static")
                     `(,xz "static")
                     libdeflate-static
                     htslib-static)))))

wfmash-gcc-static-git ;; default optimized static deployment build
