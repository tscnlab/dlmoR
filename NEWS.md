# dlmoR (development version)

- `calculate_dlmo()` is now the single entry point for all DLMO methods.
  Choose the method with `method = "hockeystick"` (default) or `method = "3k"`,
  and pass the settings of that method as named arguments (for example
  `threshold` for the hockey-stick method, `n_baseline` for the 3k method).
  Existing hockey-stick calls using named arguments, such as
  `calculate_dlmo(file_path = f, threshold = 5)`, remain compatible.
  Positional calls are not: the third argument is now `method`, so
  `calculate_dlmo(profile, NULL, 3)` must become
  `calculate_dlmo(profile, threshold = 3)`. Settings the chosen method does
  not use now give an informative error.
- New argument `dlmo_off` (default `FALSE`) to request the dim-light melatonin
  offset (DLMOoff) instead of the onset. DLMOoff is not implemented yet, so
  `dlmo_off = TRUE` stops with an error for now.
- Source files in `R/` are grouped by method with file-name prefixes:
  `hockeystick-*.R`, `3k-*.R` and `shared-*.R` (code used by several methods).
  Offset functions will go to `<method>-offset-*.R`.
- Added the 3k threshold method for DLMO estimation (Voultsios et al., 1997;
  Molina & Burgess, 2011) with `calculate_dlmo(method = "3k")`.
  - Individual threshold: mean of the first `n_baseline` samples (default 3) plus
    `k_sd` standard deviations (default 2); optional `min_threshold` floor.
  - DLMO: first linearly interpolated crossing above the threshold that stays above
    it for `min_duration` hours (default 2).
  - Edge cases stop with an error: invalid input (too few samples, duplicated time
    stamps, negative values, invalid arguments) and profiles without a DLMO
    (`no_crossing`, `not_sustained`, `truncated`). The "no DLMO" errors have class
    `dlmo_3k_error` and carry the partial result, including the plot, as `result`.
  - Warnings for a threshold that may lie within assay noise (baseline SD of 0,
    baseline values at the assay floor) and for long sampling gaps around the
    DLMO. All issues are also returned in `issues` and `flags`.
- New building blocks: `check_3k_input()`, `compute_3k_threshold()`,
  `find_threshold_intervals()`, `select_dlmo_3k()`, `check_3k_baseline()` and
  `plot_profile_3k()`.
- Added unit tests (testthat) for the 3k method.

# dlmoR 0.2.0

- Improved repeated threshold-crossing handling via `interval_limit`.
  - Single-point above-threshold excursions are treated as blips and ignored.
  - Multi-point candidate rises are compared using `interval_limit`.
  - Later candidate rises within `interval_limit` can replace earlier candidates.
  - Later candidate rises outside `interval_limit` do not replace the earlier selected rise.
- Updated numeric `interval_limit` handling so fractional hour values, such as `0.5`, work as documented.
- Added `threshold` and `interval_limit` to the object returned by `calculate_dlmo()`.
- Improved handling of sparse or poorly matched profiles by giving more informative errors when too few fitting points remain.
- Stabilized parallelogram-based ascending-segment truncation and plotting.
- Added raw-profile plotting support with `plot_raw_profile()` and `plot_profile(..., mode = "raw")`.
- Added a tutorial-style pkgdown vignette covering raw plots, thresholds, `interval_limit`, residual explorer interpretation, parallelogram truncation, and batch workflows.
- Added anonymized example profiles and residual explorer figures used in the vignette.
- Expanded documentation for `calculate_dlmo()`, `plot_profile()`, `plot_raw_profile()`, and repeated-crossing behavior.

# dlmoR 0.1.0

- Initial public release of the `dlmoR` R package.
- Implements the Danilenko et al. (2014) hockey-stick method for estimating dim-light melatonin onset (DLMO).
- Includes plotting utilities and example data.
