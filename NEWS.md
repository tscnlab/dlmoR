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
