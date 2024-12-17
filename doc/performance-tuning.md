# Performance tuning

Tests were done on my desktop with AMD Ryzen 7 3700X 8-Core Processor that has AVX2.

The default builds for the dynamic version:

```
guix build -f guix.scm
time /gnu/store/27x7aviy4l9f6ldyd4k0i105rpnys9nj-wfmash-git-0.21-HEAD.45c34b8/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m10.271s
user    1m22.931s
sys     0m0.192s
```

The default builds for the static version:


```
guix build -f guix-static.scm
time /gnu/store/zf8gyy6v5xnnaan6z4fhdwkahwbhagx5-wfmash-git-0.21-HEAD.45c34b8/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m9.352s
user    1m21.728s
sys     0m0.596s
```

Performance went up using EXTRA_FLAGS, essentially using ` -fopenmp  -Ofast -march=x86-64-v3 -flto -DNDEBUG`. We can combine the command into

```
GUIX_PATH=$(guix build -f guix-static.scm) && echo "*** Running $GUIX_PATH" && ls -l $GUIX_PATH/bin/wfmash && time $GUIX_PATH/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m8.947s
user    1m1.943s
sys     0m0.212s
```

# Testing build options

We'll run a build with different options. Note that the following commands are running in a container:

```
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_OPTIMIZED=ON .. && make -j 12 VERBOSE=1 && time ./bin/wfmash ../data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
```

```
-fopenmp -DNDEBUG -O3 -mavx2 -funroll-all-loops -fPIC -MD -MT
real    0m11.353s
user    1m8.639s
sys     0m0.288s
-fopenmp -DNDEBUG -O3 -march=x86-64-v3 -funroll-all-loops
real    0m10.460s
user    1m5.980s
sys     0m0.228s
-fopenmp -DNDEBUG -Ofast -march=x86-64-v3 -funroll-all-loops
real    0m10.860s
user    1m5.510s
sys     0m0.608s
-fopenmp -DNDEBUG -O3 -march=native -funroll-all-loops
real    0m10.268s
user    1m3.388s
sys     0m0.549s
-DNDEBUG -O3 -march=native -funroll-all-loops # try w.o. openmp
real    0m10.858s
user    1m4.946s
sys     0m0.232s
-fopenmp -DNDEBUG -O3 -march=native -funroll-all-loops -flto -flto=auto -fno-fat-lto-objects # try lto
real    0m10.358s
user    1m3.675s
sys     0m0.613s
-fopenmp  -DNDEBUG -O3 -march=native -funroll-all-loops -fprofile-use=/export/local/home/wrk/iwrk/opensource/code/pangenome/wfmash/build/../pgo -flto=auto -fno-fat-lto-objects
real    0m10.959s
user    1m6.400s
sys     0m0.668s
```

Next the static build, using guix-static.scm

```
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON -DBUILD_OPTIMIZED=ON ..
```

```
-fopenmp  -DNDEBUG -Ofast -march=x86-64-v3
real    0m10.147s
user    1m11.245s
sys     0m0.220s
-fopenmp -DNDEBUG -O3 -march=native -funroll-all-loops -flto=auto -fno-fat-lto-objects
real    0m9.949s
user    1m3.712s
sys     0m0.284s
-fopenmp -DNDEBUG -Ofast -march=x86-64-v3 -flto=auto -fno-fat-lto-objects
real    0m9.651s
user    1m3.924s
sys     0m0.381s
```

Running the last outside the container:

```
-fopenmp -DNDEBUG -Ofast -march=x86-64-v3 -flto=auto -fno-fat-lto-objects
real    0m9.447s
user    1m3.644s
sys     0m0.320s
-fopenmp -DNDEBUG -O3 -march=native -funroll-all-loops -flto=auto -fno-fat-lto-objects
real    0m9.448s
user    1m2.373s
sys     0m0.276s
-fopenmp -DNDEBUG -Ofast -march=native -flto=auto -fno-fat-lto-objects
real    0m9.248s
user    1m0.113s
sys     0m0.684s
```

shows the difference is not too bad. O3 and Ofast appear to be similar. march=native does best on Ryzen.

# Profiling

You can use the guix shell with

```
guix shell -L . -C -D -F wfmash-gcc-debug-git
```

Running a profiler I added gperftools support to cmake with `-DPROFILER=ON`. Run CMake in Debug mode with that switch

```
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DPROFILER=ON ..
make -j 12
```

and run the profiler on a test case with

```
ctest --verbose -R wfmash-time-LPA
pprof --text ./bin/wfmash ../wfmash.prof
Using local file ./bin/wfmash.
Using local file ../wfmash.prof.
Total: 8435 samples
    3019  35.8%  35.8%     3019  35.8% wavefront_bialign_breakpoint_indel2indel
    2228  26.4%  62.2%     2228  26.4% wavefront_compute_affine2p_idm
     560   6.6%  68.8%      560   6.6% wavefront_bialign_breakpoint_m2m
     537   6.4%  75.2%      537   6.4% wavefront_extend_matches_packed_kernel (inline)
     496   5.9%  81.1%      872  10.3% wavefront_extend_matches_packed_end2end_max
     134   1.6%  82.7%      134   1.6% _mm_pause (inline)
     122   1.4%  84.1%      122   1.4% wavefront_compute_trim_ends
      99   1.2%  85.3%       99   1.2% wavefront_compute_init_ends_wf_higher
      78   0.9%  86.2%       78   0.9% wavefront_compute_init_ends_wf_lower
      76   0.9%  87.1%      215   2.5% errmod_cal@@HTSLIB_1.4
      76   0.9%  88.0%       76   0.9% skch::CommonFunc::makeUpperCaseAndValidDNA
```

# Conclusion

With a bit of tweaking a 10-20% speed gain is easily possible on my Ryzen. Native compilation, openmp, lto and the static build appears to have the largest impact. PGO is, somewhat surprisingly, detrimental. Running outside a container is faster than running inside a container.

From above profiler work we can see most of the time is spent in wavefront.

For future work:

1. DONE: We ought to run a profiler to validate all the effort is going into wfa's kernels.
1. DONE: We could look at clang+LLVM performance because their optimizations are a bit different from gcc. But I don't expect much difference.
1. The AMD Genoa's we have in Octopus have AVX512. WFA is not yet optimized for that target! There could be some gains.
1. Probably is worth trying the GPU version too - at least for the supercomputers.
1. As we are pumping data has anyone looked at mmap support?
