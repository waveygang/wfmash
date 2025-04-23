# Performance tuning

Tests were done on my desktop with AMD Ryzen 7 3700X 8-Core Processor that has AVX2.

The default builds for the dynamic version:

```
guix build -L . wfmash-gcc-git --without-tests=wfmash-gcc-git
 time /gnu/store/svwjmgvh1j21wl3dc318vfjp46crrc42-wfmash-gcc-git-v0.22.0-177-ga75a6917/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m9.368s
user    1m9.595s
sys     0m1.148s
```

which is faster than the earlier 2024-12-12 version

```
guix build -f guix.scm
time /gnu/store/27x7aviy4l9f6ldyd4k0i105rpnys9nj-wfmash-git-0.21-HEAD.45c34b8/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m10.271s
user    1m22.931s
sys     0m0.192s
```

Now alto trying with native (on my Ryzen)

```
guix build -L . wfmash-gcc-git --without-tests=wfmash-gcc-git --tune=native
real    0m7.604s
user    0m58.019s
sys     0m0.900s
```

The default builds for the static version:

```
guix build -L . wfmash-gcc-static-git --without-tests=wfmash-gcc-static-git --tune=native
real    0m5.971s
user    0m42.343s
sys     0m1.811s
```

Note the speed gain may be partly because of getting rid of position independent code.
This static build is almost a doubling of speed with the earlier 2024-12-12 version

```
guix build -f guix-static.scm
time /gnu/store/zf8gyy6v5xnnaan6z4fhdwkahwbhagx5-wfmash-git-0.21-HEAD.45c34b8/bin/wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > /dev/null
real    0m9.352s
user    1m21.728s
sys     0m0.596s
```

Nice!

On a GENOA AMD EPYC 9274F 24-Core Processor the generic version runs at

```
real    0m4.473s
user    0m31.940s
sys     0m1.044s
```

and the znver2 on that machine (for model 17)

```
real    0m3.350s
user    0m23.337s
sys     0m1.627s
```

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

And -- more relevant -- with optimizations `-fopenmp -g  -DNDEBUG -Ofast -march=native -flto=auto -fno-fat-lto-objects -fPIC -MD -MT`

```
cmake -DCMAKE_BUILD_TYPE=Release -DPROFILER=ON ..
make clean
make -j 12 VERBOSE=1
rm ../wfmash.prof
ctest --verbose -R wfmash-time-LPA
pprof --text ./bin/wfmash ../wfmash.prof
Using local file ./bin/wfmash.
Using local file ../wfmash.prof.
Total: 3754 samples
     888  23.7%  23.7%      888  23.7% wavefront_bialign_breakpoint_indel2indel.localalias
     735  19.6%  43.2%      735  19.6% align::Aligner::single_reader_thread
     414  11.0%  54.3%      414  11.0% wavefront_extend_matches_packed_end2end_max.localalias
     278   7.4%  61.7%      278   7.4% inflateBackEnd@@ZLIB_1.2.0
     277   7.4%  69.0%      818  21.8% errmod_cal@@HTSLIB_1.4
     191   5.1%  74.1%      191   5.1% wavefront_compute_affine2p_idm.localalias
     178   4.7%  78.9%      432  11.5% bgzf_getc@@HTSLIB_1.0
     167   4.4%  83.3%      167   4.4% wavefront_bialign_breakpoint_m2m.localalias
      81   2.2%  85.5%       81   2.2% wavefront_extend_matches_packed_end2end.localalias
      72   1.9%  87.4%     2091  55.7% align::Aligner::worker_thread [clone .isra.0]
      69   1.8%  89.2%       69   1.8% wavefront_compute_trim_ends.localalias.lto_priv.0
      39   1.0%  90.3%       39   1.0% wavefront_compute_init_ends.localalias.lto_priv.0
      33   0.9%  91.2%       33   0.9% crc32_z@@ZLIB_1.2.9
```

Next basic all2all test runs as `wfmash -t 8 data/scerevisiae8.fa.gz > all2all.paf`.
Optimizations `-fopenmp -g  -DNDEBUG -Ofast -march=native -flto=auto -fno-fat-lto-objects -fPIC -MD -MT`

```
cmake -DCMAKE_BUILD_TYPE=Release -DPROFILER=ON ..
make clean
make -j 12 VERBOSE=1
rm ../wfmash.prof # remove older runs
ctest --verbose -R all2all # writes to ../wfmash.prof
pprof --text ./bin/wfmash ../wfmash.prof
Using local file ./bin/wfmash.
Using local file ../wfmash.prof.
Total: 58878 samples
   18686  31.7%  31.7%    18693  31.7% align::Aligner::single_reader_thread
   16577  28.2%  59.9%    16577  28.2% wavefront_bialign_breakpoint_indel2indel.localalias
    6573  11.2%  71.1%     6573  11.2% wavefront_extend_matches_packed_end2end_max.localalias
    2871   4.9%  75.9%     2871   4.9% wavefront_compute_affine2p_idm.localalias
    2841   4.8%  80.8%     2841   4.8% wavefront_bialign_breakpoint_m2m.localalias
    1223   2.1%  82.8%     1223   2.1% wavefront_extend_matches_packed_end2end.localalias
    1042   1.8%  84.6%     1333   2.3% skch::CommonFunc::addMinmers
```

When we profile all tests together we get

```
ctest
Test project /export/local/home/wrk/iwrk/opensource/code/pangenome/wfmash/build
    Start 1: wfmash-time-LPA
1/7 Test #1: wfmash-time-LPA .......................................   Passed   10.17 sec
    Start 2: wfmash-subset-LPA-to-SAM
2/7 Test #2: wfmash-subset-LPA-to-SAM ..............................   Passed   14.14 sec
    Start 3: wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF
3/7 Test #3: wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF ...   Passed   29.08 sec
    Start 4: wfmash-short-reads-500bps-to-SAM
4/7 Test #4: wfmash-short-reads-500bps-to-SAM ......................   Passed   73.20 sec
    Start 5: wfmash-short-reads-255bps-to-PAF
5/7 Test #5: wfmash-short-reads-255bps-to-PAF ......................   Passed    0.92 sec
    Start 6: wfmash-input-mapping
6/7 Test #6: wfmash-input-mapping ..................................   Passed   11.21 sec
    Start 7: wfmash-all2all
7/7 Test #7: wfmash-all2all ........................................   Passed  131.95 sec

100% tests passed, 0 tests failed out of 7

Total Test time (real) = 270.68 sec
wrk@napoli /export/local/home/wrk/iwrk/opensource/code/pangenome/wfmash/build [env]$ pprof --text ./bin/wfmash ../wfmash.prof
Using local file ./bin/wfmash.
Using local file ../wfmash.prof.
Total: 52257 samples
   15850  30.3%  30.3%    15850  30.3% wavefront_bialign_breakpoint_indel2indel.localalias
    9844  18.8%  49.2%     9844  18.8% std::__atomic_base::load (inline)
    3804   7.3%  56.4%     6340  12.1% wavefront_extend_matches_packed_end2end_max.localalias
    3526   6.7%  63.2%     3526   6.7% wavefront_extend_matches_packed_kernel (inline)
    2846   5.4%  68.6%     2846   5.4% wavefront_bialign_breakpoint_m2m.localalias
    2753   5.3%  73.9%     2753   5.3% wavefront_compute_affine2p_idm.localalias
```

which is not that different from all2all.

## Automatically running the profiler

If you use Guix you can simply run the profiler for all2all from a checked out git repo with the command:

```
guix build -L . wfmash-gcc-profile-git
```

This will build wfmash with optimization flags and run the profiler automatically.
The tests run for a minute or two.
If you scroll a bit up in the output you can see

```
[wfmash] Performing all-vs-all mapping including self mappings.
1/1 Test #7: wfmash-all2all ...................   Passed  121.70 sec
(...)
Total: 45863 samples
   11809  25.7%  25.7%    11814  25.8% align::Aligner::single_reader_thread
   11529  25.1%  50.9%    11529  25.1% wavefront_bialign_breakpoint_indel2indel
    6363  13.9%  64.8%     6363  13.9% wavefront_extend_matches_packed_end2end_max
    2837   6.2%  70.9%     2837   6.2% wavefront_compute_affine2p_idm
    2040   4.4%  75.4%     2040   4.4% wavefront_bialign_breakpoint_m2m
    1275   2.8%  78.2%     1275   2.8% wavefront_extend_matches_packed_end2end
    1186   2.6%  80.8%     1485   3.2% skch::CommonFunc::addMinmers
     784   1.7%  82.5%      784   1.7% wavefront_compute_trim_ends
     627   1.4%  83.8%      627   1.4% inflateBackEnd@@ZLIB_1.2.0
```

Update the sources and simply run again! It is the coolest thing.

Above profiling was done for the all2all test. We switched to running

```scheme
(invoke "ctest" "--verbose" "-R" "wfmash-time-LPA")
(invoke "ctest" "--verbose" "-R" "wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF")
```

Now running with native and the 'Generic' build, i.e., `-DNDEBUG -O3 -funroll-all-loops -fopenmp -g -flto=auto -fno-fat-lto-objects -fPIE`

```
guix build -L . wfmash-gcc-profile-git --tune=native
Total: 3142 samples
    1178  37.5%  37.5%     1635  52.0% skch::CommonFunc::addMinmers
     289   9.2%  46.7%      289   9.2% MurmurHash3_x64_128 [clone .constprop.0]
     286   9.1%  55.8%     1158  36.9% skch::Map::mapSingleQueryFrag
     239   7.6%  63.4%      239   7.6% std::__adjust_heap [clone .constprop.0]@4a7f70
     213   6.8%  70.2%      213   6.8% std::__push_heap [clone .constprop.0]
phase `run-profiler' succeeded after 40.6 seconds
```

The non-native version ran in 42 seconds, not much of a difference there.

Above output differs when running with much slower 'Debug' statements:

```
Total: 28620 samples
    1193   4.2%   4.2%     1193   4.2% __gnu_cxx::__normal_iterator::__normal_iterator
    1162   4.1%   8.2%     1162   4.1% skch::CommonFunc::reverseComplement
     911   3.2%  11.4%      911   3.2% MurmurHash3_x64_128
     730   2.6%  14.0%     1136   4.0% __gnu_cxx::__normal_iterator::operator+@4f1e6c
     629   2.2%  16.2%      629   2.2% __gnu_cxx::__normal_iterator::operator*@4dcb72
     596   2.1%  18.2%     1370   4.8% std::__get_helper@49ced1
```

Funny thing is that here is no wavefront in the full list, even though it builds with

```
-g -fno-inline  -g -fPIC   -Wfatal-errors -MD -MT
```

Removing no-inline makes it less noisy

```
    1154   4.2%   4.2%     1154   4.2% skch::CommonFunc::reverseComplement
     881   3.2%   7.4%      881   3.2% MurmurHash3_x64_128
     651   2.4%   9.7%     1456   5.3% std::__get_helper@4a25f1
     637   2.3%  12.0%      962   3.5% __gnu_cxx::__normal_iterator::operator+@4f70aa
     597   2.2%  14.2%      597   2.2% __gnu_cxx::__normal_iterator::operator*@4e1f98
     579   2.1%  16.3%     2035   7.4% std::get@486b3b
```

Obviously reverseComplement and MurmurHash3 are slow, but it is puzzling we see no wavefront in the list.
Most of the time in the test goes to [wfmash::align]. `skch`, btw, is the jaccard mapper, part of MashMap.
The Jaccard Similarity Index is a measure of the similarity between two sets of data. For an [example](https://www.statology.org/jaccard-similarity/).
From the MashMap paper: one class of algorithms for fast, approximate mapping relies on ideas originally developed for ﬁnding
similarities between web documents. Broder (1997) proved that an unbiased estimate of the Jaccard
similarity coefﬁcient between two sets can be computed efﬁciently using a subset of hashed elements called
a sketch. Schleimer et al. (2003) proposed the winnowing algorithm, which picks a minimum hashed item
(also known as a minimizer; Roberts et al., 2004) from each consecutive window of text as a means to more
quickly estimate local similarity between web documents. These ideas have been used to develop new
mapping and assembly algorithms for long reads, such as the MinHash Alignment Process (Berlin et al.,
2015), minimap (Li, 2016), and BALAUR (Popic and Batzoglou, 2016).

Just running LPA:

```
1: [wfmash::mashmap] Processing 1 target subsets (≈2317910bp/subset)
1: [wfmash::mashmap] Processing subset 1/1 (mapping)
1: [wfmash::mashmap] indexing [0.0% complete, 0/100 units, 0s elapsed]
1: [wfmash::mashmap] indexing [100.0% complete, 100/100 units, 0s elapsed]
1: [wfmash::mashmap] Processed 8 sequences (0 skipped, 2317910 total bp), 20552 unique hashes, 258475 windows
1: [wfmash::mashmap] Filtered 179727/438202 k-mers occurring > 87 times (target: 0.02%)
1: [wfmash::mashmap] reference index computed in 1.3s
1: [wfmash::mashmap] building index data structures...
1: [wfmash::mashmap] mapping [0.0% complete, 0/2317910 units, 0s elapsed]
1: [wfmash::mashmap] mapping [2.1% complete, 49000/2317910 units, 0s elapsed]
1: [wfmash::mashmap] mapping [100.0% complete, 2317910/2317910 units, 1s elapsed]
1: [wfmash::mashmap] Mapped query in 3.0s, results saved to: /tmp/guix-build-wfmash-gcc-profile-git-v0.22.0-184-g8f6a1d64.drv-0/source/wfmash-GWRMuL
1: [wfmash::align] Reference = [data/LPA.subset.fa.gz]
1: [wfmash::align] Query = [data/LPA.subset.fa.gz]
1: [wfmash::align] Mapping file = /tmp/guix-build-wfmash-gcc-profile-git-v0.22.0-184-g8f6a1d64.drv-0/source/wfmash-GWRMuL
1: [wfmash::align] Alignment identity cutoff = 0.00%
1: [wfmash::align] Alignment output file = /dev/stdout
1: [wfmash::align] Alignment scoring: [mismatch, gap1(o,e), gap2(o,e)[] = [5, 8, 2, 24, 1]
1: [wfmash::align] Successfully loaded thread-safe FASTA indices
1: [wfmash::align] time spent loading the reference index: 0.00 sec
1: [wfmash::align] Found 861 mapping records for alignment (13012371 query bp, 14815481 target bp)
1: [wfmash::align] aligning [0.7% complete, 85000/13012371 units, 0s elapsed]
1: [wfmash::align] aligning [18.8% complete, 2446256/13012371 units, 10s elapsed]
1: [wfmash::align] aligning [47.7% complete, 6201280/13012371 units, 20s elapsed]
1: [wfmash::align] aligning [65.2% complete, 8486876/13012371 units, 30s elapsed]
1: [wfmash::align] aligning [86.7% complete, 11276414/13012371 units, 40s elapsed]
1: [wfmash::align] aligning [100.0% complete, 13012371/13012371 units, 47s elapsed]
1: [wfmash::align] total aligned records = 861, total aligned bp = 14294162, completed in 53 seconds
1: [wfmash::align] time spent computing the alignment: 53.1 sec
```

and we get wavefront output alright

```
    6639  45.6%  45.6%     6639  45.6% wavefront_compute_affine2p_idm
    2949  20.3%  65.9%     2949  20.3% wavefront_extend_matches_packed_end2end_max
    2491  17.1%  83.0%     2491  17.1% wavefront_bialign_breakpoint_indel2indel
     408   2.8%  85.8%      408   2.8% wavefront_compute_affine2p_idm_piggyback
     364   2.5%  88.3%      364   2.5% wavefront_bialign_breakpoint_m2m
     189   1.3%  89.6%      189   1.3% wavefront_compute_trim_ends
     135   0.9%  90.5%      135   0.9% wavefront_extend_matches_packed_end2end
```

Running yeast we get:

```
3: [wfmash::mashmap] indexing [0.0% complete, 0/100 units, 0s elapsed]
3: [wfmash::mashmap] indexing [100.0% complete, 100/100 units, 0s elapsed]
3: [wfmash::mashmap] Processed 136 sequences (0 skipped, 96255507 total bp), 670229 unique hashes, 7449528 windows
3: [wfmash::mashmap] Filtered 7245/7456773 k-mers occurring > 1491 times (target: 0.02%)
3: [wfmash::mashmap] reference index computed in 117.4s
3: [wfmash::mashmap] building index data structures...
3: [wfmash::mashmap] mapping [0.0% complete, 0/96255507 units, 0s elapsed]
3: [wfmash::mashmap] mapping [0.0% complete, 22000/96255507 units, 0s elapsed]
3: [wfmash::mashmap] mapping [5.6% complete, 5394000/96255507 units, 10s elapsed]
3: [wfmash::mashmap] mapping [11.3% complete, 10876000/96255507 units, 20s elapsed]
3: [wfmash::mashmap] mapping [17.9% complete, 17193000/96255507 units, 30s elapsed]
3: [wfmash::mashmap] mapping [23.7% complete, 22805000/96255507 units, 40s elapsed]
(...)
3: [wfmash::mashmap] mapping [100.0% complete, 96255507/96255507 units, 158s elapsed]
3: [wfmash::mashmap] Mapped query in 278.9s, results saved to: /dev/stdout
```

```
    1140   4.1%   4.1%     1140   4.1% skch::CommonFunc::reverseComplement
     835   3.0%   7.1%      835   3.0% MurmurHash3_x64_128
     705   2.5%   9.6%     1465   5.2% std::__get_helper@4a260b
     684   2.4%  12.0%     1102   3.9% __gnu_cxx::__normal_iterator::operator+@4f70c4
     606   2.2%  14.2%      606   2.2% __gnu_cxx::__normal_iterator::__normal_iterator@4e1e2c
     600   2.1%  16.3%     4203  15.0% std::__tuple_compare::__less@486b6f
(...)
     213   0.8%  47.6%    11305  40.4% skch::CommonFunc::addMinmers
```

I toyed a bit with skch::CommonFunc::reverseComplement and it turned out that, when taking out the comparisons, the speed for yeast tests on 4 threads was improved.

```
old: 3/6 Test #3: wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF ...   Passed   27.42 sec
new: 3/6 Test #3: wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF ...   Passed   22.01 sec
```

And the function is no longer a bottleneck with Yeast:

```
     833   3.1%   3.1%      833   3.1% MurmurHash3_x64_128
     678   2.5%   5.6%     1451   5.3% std::__get_helper@4a26df
     656   2.4%   8.0%      980   3.6% __gnu_cxx::__normal_iterator::operator+@4f7198
     606   2.2%  10.2%      606   2.2% __gnu_cxx::__normal_iterator::operator*@4e2086
     537   2.0%  12.2%     1988   7.3% std::get@486c29
```

The profiler turns out to be correct.

For more information on installing and running Guix see the header of [guix.scm](../guix.scm). Guix is distribution agnostic and runs on, for example, Debian. That is what we do.

## Tuning addMinmers

Above profiling suggests that addMinmers is taking a chunk of time. We need to validate if the profiler has got it right.
First I am going to slow down mapSingleQueryFrag and see if that comes on top. Next we probably need to de-inline functions in addMinMers. Also I should check the build is correct for the other libs, so performance gets tracked there too.

# Conclusion

With a bit of tweaking a 10-20% speed gain is easily possible on my Ryzen. Native compilation, openmp, lto and the static build appears to have the largest impact. PGO is, somewhat surprisingly, detrimental. Running outside a container is faster than running inside a container.

From above profiler work we can see most of the time is spent in wavefront.

For future work:

1. DONE: We ought to run a profiler to validate all the effort is going into wfa's kernels.
1. DONE: We could look at clang+LLVM performance because their optimizations are a bit different from gcc. But I don't expect much difference.
1. DONE: The AMD Genoa's we have in Octopus have AVX512. WFA is not yet optimized for that target! There could be some gains.
1. Probably is worth trying the GPU version too - at least for the supercomputers.
1. As we are pumping data has anyone looked at mmap support?

# Addendum

Some older exercises:

## Testing build options (2024-12-12)

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
