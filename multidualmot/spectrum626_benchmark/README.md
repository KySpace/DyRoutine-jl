# Spectrum 626 calibration benchmark

Run from the MultiDualMOT worktree:

```powershell
julia --project=. multidualmot/benchmark_spectral626.jl
```

The three `config.yaml` versions are calibration-only snapshots of the live data configuration:

| Folder | Method |
| --- | --- |
| `joint_offsets` | Existing shared cubic B-spline fit, source frequency offsets, and exact natural cubic residual interpolation |
| `joint_no_offsets` | The same joint fit with every first-stage frequency offset fixed to zero, followed by the same residual interpolation |
| `individual_natural` | Independent natural cubic interpolation directly through each source's detected references |

Both `common` (sources 1, 2, 3, 4, 7, 8) and `endtoend` (sources 1, 3, 4, 7, 8) are evaluated. Source order, assigned frequencies, search times and detection logic are identical across variants. The joint variants share the current six knots and regularisation 0.002. Direct natural interpolation has no regularisation parameter or additional common knots. `naturalcubicspline` is a benchmark-only model name, separate from the production monotone Hermite `cubicspline` model; these calibration-only files are inputs to the benchmark runner, not replacements for the full spectrum configuration.

## Criterion and range

At the same CH1-aligned time `t`, define disagreement as

```text
D(t) = max_s f_s(t) − min_s f_s(t)
score = maximum D(t) over the selected time interval
```

Lower scores mean greater consistency among sources. They do not establish absolute calibration accuracy. The final full curves include offsets and residuals; offsets are not removed before scoring, and optical zero times are not realigned separately.

The selected range is −750 to +750 MHz. Its time endpoints are fixed by independent natural interpolation of source 1, with linear tails. This same reference axis and time interval are used for every method, preventing a method from selecting a more favorable evaluation interval. Each tag derives its reference from its own configured source-1 peaks, which are identical in the supplied snapshots.

The score is the continuous maximum to floating-point precision, rather than the largest plotting-grid value: union all spline breakpoints, then check every interval endpoint and stationary point of every pairwise cubic difference. `plot_samples` only controls exported/visualized curve sampling and does not affect the score. A maximum marker is drawn at the analytic extremum even if it lies between exported samples.

Each component keeps the existing linear endpoint-tangent extrapolation. The score includes extrapolation where the chosen range exceeds a source's reference coverage. `maximizing_pair_extrapolated` explicitly flags maxima where either of the two responsible sources is outside its own reference-time span. Source 2 has only the +441.6, +220.8 and 0 MHz references, so much of its selected range requires extrapolation. `endtoend` provides a separate comparison without that sparsely anchored source.

“Without offsets” fixes the explicit offsets in the joint first stage. Residual interpolation remains unrestricted and can itself contain constant corrections. Direct natural interpolation needs no separate offset-on variant: adding a constant offset and interpolating the remaining residual values gives exactly the same natural curve.

## Configuration and outputs

Edit `benchmark.frequency_range_mhz` in all three configs to select another common range, and `benchmark.path_data` if the captures move. Source definitions must remain identical across variants for each tag; the runner rejects mismatched reference sets, range settings, and duplicate variant names. You may independently change the joint variants' knots and regularisation for later experiments.

The optional first command-line argument selects another root containing variant folders, each with `config.yaml`. The optional second argument sets the output directory:

```powershell
julia --project=. multidualmot/benchmark_spectral626.jl path/to/config-root path/to/output
```

Default outputs are in `results/`:

- `summary.csv`: continuous maximum, reference-axis frequency and aligned time at that maximum, upper/lower source indices, extrapolation flag, and maximum reference interpolation error.
- `disagreement.png`: all three methods on the same axis, separately for each source set.
- `disagreement_curves.csv`: sampled disagreement curves with the fixed reference-frequency axis.
- `detected_references.csv`: the actual detected peak times used by every method.
- `calibration_curves.jld2`: sampled individual final calibration curves and summary values, stored as plain arrays and metadata.
- `config_snapshots/<variant>/config.yaml`: exact configuration copies for the run.

The live data-folder `config.yaml`, original report, and production spectrum/result exports are not modified by a benchmark run.

## Reusing the Julia session

```julia
include("multidualmot/benchmark_spectral626.jl")
# Edit the YAMLs, then repeat just the benchmark:
benchmark_data = run_spectral626_benchmark(path_benchmark_configs;
    path_output=path_benchmark_output, source_cache=benchmark_source_cache)
```

The first pass loads, smooths and detects peaks from each unique capture/reference specification once. The supplied configs require six captures total; subsequent variants and reruns reuse compact detected-reference records. `benchmark_data.cases[(tag, variant)]` retains each fitted model and analytic maximum. If a capture file changes in place, run `empty!(benchmark_source_cache)` before rerunning; changing configured reference times/frequencies creates a new cache key automatically.

The production helper now accepts `fit_source_offsets=false`, also read from a calibration entry's YAML. Omission defaults to `true`, preserving the current fit. `spectral626.jl --definitions-only` loads its definitions without running workflow exports; the benchmark uses this mode in a separate module.

The analytic maximum helper has standalone tests:

```powershell
julia --project=. script/test_calibration_benchmark.jl
```
