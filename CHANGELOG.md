# Changelog

## [Unreleased]

- Adds more targets for SIMD acceleration and more complete use of instruction sets
- Reduce memory usage almost back down to [v2.2.0]
- Reduce runtime for 1.17 by about 20% and 1.18+ by 5%

## [v2.3.1] - 2026-09-04

- Fixes incorrect SIMD instruction sets being selected, fastest should now be prioritized

## [v2.3.0] - 2026-09-04

- Accelerates search with SIMD
    - This can use AVX2 (1.3-1.5x faster) and AVX-512 (1.6-2.6x faster) or fallback to SSE on x86
    - Also has support for NEON on ARM, though I have not tested performance changes there
    - Unsupported architectures will need to compile this program natively themselves
- Improves ETA latency
- Memory usage was increased (1.6x) to facilitate new optimizations
- Adds actions workflow and removes some Nix outputs

## [v2.2.0] - 2026-08-22

- Allows setting `--loaded-radius` down to 0 to allow searching for smaller regions (requested from siestol)
- Micro-optimizations from `geode-viber` experiment by @worldInColors
- Preparations for 26.3+ support

## [v2.1.1] - 2026-08-16

- Constrains search area within the world border during `--estimate` when given center coordinates

## [v2.1.0] - 2026-07-28

This updates finally adds multi-threading so you can perform your world-size searches a bit faster (maybe). I also added some minor features, fixed a few bugs, and reworked some stuff under the hood.

- Added `--threads` option. Setting this option will enable geode-finder to use as many cores as specified to accelerate both phases of the search. Set it to 0 to use all available cores. **Note that by default, geode-finder will only run with 1 core**.
- Slight optimizations to performance should net about a 25% improvement 1.17 and 1% improvement on 1.18+ on larger phase 1 searches.
- All errors should be properly emitted to the user instead of some being consumed now (thanks Undec)
- Adds a percentage completion to the progress bar. Due to changes from adding multi-threading, the ETA is now more jumpy and less accurate, so this was added to alleviate those problems.
- `geodegen.sc` script now works with the new output format and also computes
- Removed some unnecessary dependencies to speed up builds
- Adds a Nix flake

## [v2.0.0] - 2026-03-22

This update is partly a massive refactor/rewrite to make it easier to port to new versions should the geode placement configuration ever update in the future. The other part is a streamlined geode cluster searching algorithm that greatly improves the speed and provides flexibility for a parallelized approach in a future release.

- Several options have been renamed (`-h`)
    - Default settings have been tweaked (`-h`)
    - Limits have been added to prevent searching outside the MC world limit
- 1.17 and 1.18+ first pass performance improved by **3.2x** and **2.1x**
- Final results are saved in `output.json`. Use the `--output-path` option to configure the location.
- `--estimate` option calculates a rough estimate of how many clusters will be found
- Extra unit tests and Fabric mod for validating parity

## [v1.1.0] - 2025-07-14

- Added `--random-tick-radius` argument for versions 1.21.5+ (or additional cases where you wish to have a different chunk radius). Fox example, if you wish to load all geodes within a 16 simulation distance area, pass 16 to this option. Thanks Elyrio for the reminder!

## [v1.0.4] - 2024-01-23

- Fix for negative starting chunk coordinates

## [v1.0.3] - 2024-01-22

- Fixes issues with negative seed input
- Adds support for 1.20 (mostly negligible differences in budding amethyst placement)
- Adds options to set starting chunk coordinates

## [v1.0.0] - 2023-01-24

- Adds searching for the largest cluster of geodes and budding amethyst in a Minecraft world
