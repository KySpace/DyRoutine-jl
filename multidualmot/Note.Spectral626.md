# Spectrum 626 workflow

[`spectral626.jl`](spectral626.jl) converts the scope captures to JLD2 caches and produces the CH1 cycle, aligned CH2, time-to-frequency calibration, and calibrated spectrum figures. The default data directory is the `spectrum 626` folder configured near the top of the script. Its `config.yaml` controls the aligned display indices, tagged calibrations, and result traces/peak annotations.

## Processing sequence

1. **Cache scope CSV files as JLD2.** The script scans for files named `SDS#####.csv`. It reads the scope metadata and the three numeric columns (time, CH1, CH2), checks the sample count against `Record Length`, and saves the numeric vectors and metadata to a same-basename `.jld2`. An existing JLD2 is treated as the cache and its CSV is skipped.
2. **Load and smooth the channels.** Analysis loads `ch1_v`, `ch2_v`, and metadata from JLD2. Each channel is smoothed with a centered, symmetric 11-sample Gaussian kernel (σ = 2.5 samples, spanning 10 sample intervals). Replicate padding handles the capture edges. The smoothing is applied in memory after loading; it does not replace the stored arrays in the JLD2 cache.
3. **Construct the time axis.** Scope metadata supplies horizontal units and horizontal scale. The script uses the 14 horizontal divisions to map samples across the full capture duration in seconds.
4. **Find CH1 cycles and align traces.** CH1 minima and maxima are found using 2-second windows/dead-time, with the sample extremum refined within the window. Complete periods are minimum-to-minimum cycles; incomplete edge cycles are omitted from the period table. For aligned plots and calibration, the time origin is the interpolated rising CH1 zero crossing of the selected cycle (or the corresponding alignment crossing for a whole capture).
5. **Find CH2 peaks for the aligned comparison.** In the cycle spanning raw time zero, the CH2 baseline is the mean after excluding the top 1% of values. Local maxima must exceed this baseline by 5 mV before the 50 ms dead-time rule is applied; this order prevents small candidates from suppressing real peaks. The strongest candidate in each dead-time window is retained.
6. **Build each tagged calibration using YAML references.** Each `calibration` entry has a unique `tag`, a `model` (`p4` or `cubicspline`), and a `source` list. For every source, the script detects a CH2 local maximum within ±1 ms of each configured aligned time. With `model: p4`, a detected 0 MHz peak fixes that source's $t_0$; when none is supplied, $t_0$ is fitted. References within one tag are fitted together using the existing unrestricted model

   \[
   f(u)=\alpha_i u\left(1+p_2u+p_3u^3+p_4u^4\right),\qquad u=t-t_{0,i}.
   \]

   The shape parameters $p_2,p_3,p_4$ are shared across sources in that tag; the scale αᵢ is fitted independently for each source. With `model: cubicspline`, each source gets its own shape-preserving monotone cubic Hermite interpolant through its detected references, ordered by aligned time. Its endpoint slopes use one-sided, monotonicity-limited derivatives. Set `extrapolation: linear` (the default) for endpoint-tangent lines or `extrapolation: quadratic` for quadratic tails determined by the first and last three references. Quadratic tails stop at their first turning point so the time–frequency mapping remains one-to-one; they require at least three references per source. The displayed $t_0$ is the time where this curve crosses 0 MHz. Each tagged calibration has its own figure. Solid model segments span the reference points; dashed segments extrapolate toward ±800 MHz or the first turn. The results and figures remain available as `calibration_results` and `fig_calibrations` arrays in an interactive Julia session, with `calibrations_by_tag` for lookup.
7. **Reconstruct result spectra and annotate configured peaks.** Each `result` entry selects a calibration by its `calibration` tag. The script loads and smooths its JLD2 channels and aligns time to CH1. For `p4`, a calibration source reuses its fitted $t_0$ and α; another source fits them from its `ref: true` peaks while retaining the calibration's shape parameters. The p4 trace extrapolates up to the result's ±720 MHz plot limits or its first turning point. For `cubicspline`, a calibration source reuses its interpolant directly; another source builds a new interpolant from its `ref: true` peaks, using the selected extrapolation method. The spline trace stops at ±720 MHz or a quadratic-tail turning point. The initial frequency-axis major ticks are spaced by 200 MHz; zooming in selects denser ticks automatically. The spectrum uses the smoothed CH2 voltage as read from JLD2, without the temporary vertical offsets used in comparison plots and without height normalization. Non-reference peaks are display annotations only: if a frequency is configured, the script searches for the strongest peak within ±5 MHz on the plotted trace, whether or not a time is also present; otherwise it searches near the configured time. Labels include the configured tag and detected time, frequency, and voltage.

## Generated files

By default, outputs are written beside the data files:

- `SDS00014_CH1_vs_time.png` — high-resolution CH1 capture with detected extrema.
- `SDS00014_CH1_periods.csv` — complete CH1 periods and uphill/downhill durations and sample counts.
- `SDS<display indices>_CH2_aligned_periods.png` — one cycle spanning zero per capture listed in a `display` entry, aligned to the rising CH1 zero crossing. CH2 traces have temporary 0.5 V vertical offsets; CH1 is scaled by 1/8 before using the same offset and a lighter color. These display-only adjustments are not used in calibration or result reconstruction.
- `spectrum626_time_frequency_calibration_<tag>.png` — one calibration and reference-residual figure per calibration tag.
- `spectrum626_calibrated_results.png` — YAML-selected calibrated result traces and peak annotations, with default frequency limits −720 to +720 MHz and 200 MHz major ticks.

## Running

```powershell
julia --project=. multidualmot/spectral626.jl
```

An optional first argument selects another data folder. `--skip-conversion` uses only existing JLD2 caches and avoids scanning/reading CSV files. `--calibration-only` also skips the SDS00014-specific CH1 plot, period CSV, and aligned multi-file cycle plot, but still produces the calibration and result figures. These flags can be combined.
