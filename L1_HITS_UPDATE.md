# L1 Hits Default Value Update

## Summary
Updated the default value for L1-hits (-H) parameter from 1 to 3 in wfmash.

## Changes Made

### File: `src/interface/parse_args.hpp`

1. **Line 82**: Updated help text
   ```cpp
   // Before:
   args::ValueFlag<int> min_hits(filtering_opts, "INT", "min hits for L1 filtering [1]", {'H', "l1-hits"});
   
   // After:
   args::ValueFlag<int> min_hits(filtering_opts, "INT", "min hits for L1 filtering [3]", {'H', "l1-hits"});
   ```

2. **Line 641**: Updated default value
   ```cpp
   // Before:
   map_parameters.minimum_hits = 1; // default minimum
   
   // After:
   map_parameters.minimum_hits = 3; // default minimum
   ```

## Impact
- The L1 filtering will now require at least 3 hits by default (instead of 1)
- This makes the filtering more stringent, reducing spurious mappings
- Users can still override this with `-H` parameter if needed

## Verification
```bash
$ ./bin/wfmash --help | grep l1-hits
-H[INT], --l1-hits=[INT]         min hits for L1 filtering [3]
```

The program now uses 3 as the default minimum hits for L1 filtering.