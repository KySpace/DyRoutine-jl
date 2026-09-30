# Multidualmot workflow overview

## Overall workflow

The project treats each dataset as:

```text
Data/
└── Dataset name/
    └── 162-164/
        ├── config.yaml
        ├── liferes MMDD runXX.mat
        ├── ...
        └── generated SVG/PNG figures
```

A dataset-specific runner defines how the YAML variables should be interpreted, then delegates common reading, reshaping, statistics, and plotting to shared scripts.

The main entry point is [`multidualmot_runner.jl`](multidualmot_runner.jl). It currently runs all analyses, regenerates the combined SVG/PNG table, and creates a new OneNote page.

## `config.yaml` structure

Each `config.yaml` is a top-level sequence of processing entries:

```yaml
- tag: "t-balanced"
  data:
    - tbiasmot: [0]
      t_hold: [50, 100, 200, 500, 1000, 2000]
      loadcfg: [DDM, DIS]
      istp: [162, 164]
      vars: [rep, tbiasmot, t_hold, loadcfg, istp]
      source:
        - "liferes 0827 run54.mat"
```

The hierarchy is:

```text
config.yaml
└── processing entry
    ├── tag
    └── data
        ├── rectangular data block 1
        ├── rectangular data block 2
        └── ...
```

Important rules:

- Every processing entry needs a unique, nonempty `tag`.
- Every processing entry needs a nonempty `data` sequence.
- Each `data` entry describes one rectangular acquisition block.
- `vars` lists the actual acquisition nesting order exactly once.
- Earlier variables in `vars` vary more slowly; later variables vary faster.
- `rep` may appear anywhere in `vars`, although it is usually first.
- `rep` is inferred from the MAT sample count and is not assigned values in YAML.
- `source` contains filenames local to the isotope-pair folder.
- Source filenames must follow `liferes MMDD runXX.mat`.

### One `source` list versus multiple `data` entries

These mean different things.

Multiple files in one source list:

```yaml
source:
  - "liferes 0712 run08.mat"
  - "liferes 0712 run10.mat"
```

are concatenated in that order into one flat acquisition stream, then reshaped together.

Separate `data` entries:

```yaml
data:
  - ...
    source: ["liferes 0827 run54.mat"]
  - ...
    source: ["liferes 0827 run55.mat"]
```

are each reshaped and analyzed independently first. Their multidimensional arrays are then aligned by variable values and appended along the repetition axis. Combinations absent from one block remain `missing`.

Use separate `data` entries when runs have different rectangular variable coverage.

## Runner-defined variable schemas

The shared schemas are in [`dualmotcommons.jl`](dualmotcommons.jl).

Examples:

| Dataset | Canonical variables |
|---|---|
| MOT/CMOT lifetime | `β_MOT`, `t_hold`, `loadcfg`, `istp` |
| MOT loading 421 | `β_MOT`, `t_load`, `loadcfg`, `istp` |
| MOT loading 626 | `t_load`, `loadcfg`, `istp` |
| MOT loading balance | `β_MOT`, `loadcfg`, `istp` |
| ODT B field | `ib`, `loadcfg`, `istp` |

The YAML name `tbiasmot` is converted into the Julia field `β_MOT`.

Isotopes and load configurations are converted to symbols. Isotope symbols must be constructed as:

```julia
Symbol("162")
```

not `:162`, because Julia interprets `:162` as an integer expression rather than a `Symbol`.

## How `liferes.mat` is read

The central reader is `read_cres_fields` in [`dualmotcommons.jl`](dualmotcommons.jl).

For each MAT file it reads:

```julia
file["liferes"]["cres"]
```

If `liferes["baddata"]` is present, its linear indices are applied to the
selected number field (atomnum or pixsum) before files are concatenated. Those
samples remain in the acquisition stream as `missing`.

For number-evolution analysis, it extracts the selected number field and these
fit-quality fields from every `cres` entry:

```text
atomnum (default; lifetime runners use atomnum)
sigmax
sigmay
```

Each field must be a numeric scalar. Values from all filenames in the block's `source` list are concatenated in source-list order.

The reader handles `cres` stored as:

- a MATLAB struct array;
- a single dictionary-like struct;
- or an ordinary vector of structures.

### Selectors

A `data` entry may include a `selector` sequence to mark complete repetitions
or variable conditions as missing after its sources have been combined and
reshaped, but before separate `data` entries are combined. Each selector maps
one configured variable to either a `value` or `index` predicate. Predicates
are Julia functions written as strings, for example:

```yaml
selector:
  - rep:
      value: "a -> a <= 4"
```

`value` receives the configured axis values; `index` receives their 1-based
positions. Samples outside the predicate remain present as `missing`, keeping
all acquisition axes aligned.

Other shot rejection is determined from:

- `atomnum`;
- `sigmax`;
- `sigmay`;
- and the bounds configured by the dataset runner.

## Reshaping the flat acquisition stream

For one data block:

1. Compute the number of non-repetition conditions:

   ```julia
   n_variation = product(length of every configured variable except rep)
   ```

2. Infer repetitions:

   ```julia
   n_rep = sample_count ÷ n_variation
   ```

3. Validate that the sample count is rectangular.

4. Reshape according to the YAML `vars` acquisition order.

5. Permute the result into the runner's canonical order, normally:

   ```text
   rep × canonical variable 1 × canonical variable 2 × ...
   ```

This work happens in `calc_num_evol_block`.

Most runners require a perfectly rectangular sample count. The balance runner has an opt-in partial-repetition mode that can pad an incomplete final repetition with `missing`.

## Rejected shots and missing values

A shot becomes `missing` if the selected number field is nonfinite or outside
`0:num_max_num`, or if either fitted size is nonfinite or outside its configured
bounds:

```julia
bounds_sigmax_num[1] ≤ sigmax ≤ bounds_sigmax_num[2]
bounds_sigmay_num[1] ≤ sigmay ≤ bounds_sigmay_num[2]
0 ≤ selected number field ≤ num_max_num
```

Nonfinite fitted sizes fail the size masks.

Rejected shots are not deleted from the flattened acquisition. They are kept as `missing` so all acquisition axes remain aligned.

For each condition, the code computes over nonmissing repetitions:

```text
mean atom number
sample standard deviation
number of valid repetitions
```

With fewer than two valid repetitions, the sample deviation is `NaN`.

## Combining multiple rectangular blocks

`combine_num_evol_blocks`:

1. Verifies that every block has the same set of axes.
2. Builds an ordered union of each non-repetition axis's values.
3. Sums the repetition counts.
4. Allocates a combined `Union{Missing,Float64}` array.
5. Aligns each block by its configured variable values.
6. Places each block into a separate range on the combined repetition axis.
7. Leaves unavailable combinations as `missing`.

This permits, for example, one source to contain both t- and n-balanced panels while another contains only the n-balanced panel.

## Shared number-evolution pipeline

The main stages are:

```text
Dataset runner
    ↓
read_num_evol_runinfos
    ↓
read config.yaml and validate source definitions
    ↓
anlz_num_evol.jl
    ↓
calc_num_evol_block for each data entry
    ↓
combine_num_evol_blocks
    ↓
means, deviations, valid repetition counts
    ↓
anlz_num_evol_output.jl
    ↓
linear/log SVG and 4× PNG figures
```

Key files:

- [`dualmotcommons.jl`](dualmotcommons.jl): YAML/MAT readers, validation, reshaping, styles, and fits.
- [`anlz_num_evol.jl`](anlz_num_evol.jl): common block processing and statistics.
- [`anlz_num_evol_output.jl`](anlz_num_evol_output.jl): standard evolution figures.
- [`anlz_num_evol_decay.jl`](anlz_num_evol_decay.jl): CMOT/MOT fitting figures and fit records.

A runner is responsible for defining:

```text
data root
pair order
variable schema
x variable
size/number rejection bounds
plot labels and transformations
fit model/configuration, when relevant
```

It then loops over the parsed `runinfo` entries and includes the shared processing scripts.

## Special analysis paths

### MOT loading 626

This uses three independent rectangular blocks:

```text
DCS for both isotopes
SCS for isotope 1
SCS for isotope 2
```

Each block is processed independently before curves and ratios are combined. It uses [`anlz_cmpr_loadcfg.jl`](anlz_cmpr_loadcfg.jl), rather than the ordinary combined-block plotting path.

### Lifetime fits

CMOT uses a κ-only model:

```text
τ → ∞
```

and currently fits only:

```julia
t_hold < 2.0 s
```

Later CMOT points remain visible as hollow markers, while the fitted curve is extended across the complete displayed time range.

MOT uses a τ-only model:

```text
κ = 0
```

and currently fits all valid points.

Each fit figure contains:

- a log-y axis;
- an equal-sized linear-y axis below it;
- one parameter listing on the right.

Fit records are saved as timestamped JLD2 files under:

```text
CMOT lifetime/Fit results/
MOT lifetime/Fit results/
```

CMOT and MOT lifetime runners use `atomnum` as their number field. Their decay
axes and κ units name the selected signal. Pairwise CMOT decay comparisons use
`1/κ` from those atom-number fits, while MOT comparisons use the fitted `τ`
values. Pairwise fitted `N₀` values are labeled in atom-number units.

The CMOT/MOT pairwise lifetime comparison exports the fitted `N₀` and the
`1/κ` or `τ` comparisons as separate figures. In OneNote, the two figures share
one comparison cell and are stacked vertically. Their workbook data use
separate sheets. Each per-pair CMOT/MOT lifetime cell also
includes a `size` figure with repetition means and sample deviations of fitted
`sigmax` and `sigmay`, using the same accepted shots and condition alignment as
the number curves. The size figure occupies the final inner-grid column.

The isotope-pair lifetime comparisons read the newest JLD2 file in each folder.

### Isotope-pair comparisons

The comparison runner writes separate t-balanced and n-balanced value and ratio
figures for MOT loading 421, CMOT decay, and MOT decay. It uses 30-second
loading values for 421 and the fitted `1/κ` and `τ` values for CMOT and MOT.
Balance biases come from `MOT loading balance/balance.csv`; β_MOT = 0 supplies
t-balance and the pair-specific balance bias supplies n-balance. For 162–164,
both comparisons use the same zero-bias data. MOT loading 626 has no balance
axis, so its 30-second values and ratios are repeated in both comparison rows.
Unavailable pair/variant data are left blank. The concluding OneNote page has
four comparison rows: t-balanced values, n-balanced values, t-balanced ratios,
and n-balanced ratios.

For MOT loading 421 comparisons, a sole tagged group is used directly. If a
pair has multiple tagged groups, exactly one must set `for_load421: true`; all
421 comparison points, including t-balanced DCS, are selected from that group.
Missing or duplicate selections raise an error.

Each MOT loading 421 pair folder also receives two log-y ratio-evolution
figures for each balance in the selected group: DIS/DCS and DDM/DIS. The two
isotopes share each figure, distinguished by isotope color. Ratio means are
ratios of the per-isotope means; uncertainty is propagated from the repetition
standard deviations assuming independent samples. In the combined table, 421
tag/balance combinations occupy separate subrows, with the inner columns
ordered linear loading, log loading, DDM/DIS, and DIS/DCS.

The t-balanced 421 number figure also includes DCS points selected from 30-second,
zero-bias data, using the standard DCS circle marker. Its y-axis limits are
controlled by `limits_y_numbers` in `anlz_isotope_pair_comparison_runner.jl`.
The t-balanced loading comparison includes a DIS/DCS ratio figure with the same
axis limits and styling as the DDM/DIS ratio figure; both ratio panels are shown
under MOT loading 421 in the OneNote comparison table.

### ODT / CMOT number comparison

[`anlz_odt_cmot_comparison.jl`](anlz_odt_cmot_comparison.jl) compares each
pair's ODT peak number with its MOT loading 421 number at the source-table
`tbiasmot` and selected `t_load`. ODT points are the maximum mean over `ib`
within each direction, load configuration, and isotope; available x and z
maxima are both retained. CMOT points use the matching 421 bias and loading
time. The log-number figure distinguishes ODT (filled) and CMOT (open) markers;
the linear figure shows ODT/CMOT with independent standard-deviation error
propagation. Pair labels include bias and loading time. The figures are added
to the ODT column at the top of the OneNote comparison rows.

The ODT source table supplies the pair-specific `tmotload` used to select the
matching 421 loading point: 20 s for all pairs except 162-163 and 163-164,
which use 15 s. These values are also recorded as metadata in each ODT config.
When a 421 config has one processing tag, that entry is used. When it has
multiple tags, exactly one must set `for_odt: true`; those selections currently
apply to 162-164 and 161-163.
The standalone number and ratio keys are generated by
[`make_odt_cmot_comparison_legend.jl`](make_odt_cmot_comparison_legend.jl) as
`[ODT.CMOT.comparison].[legend].[numbers]` and
`[ODT.CMOT.comparison].[legend].[ratios]` SVG/PNG files in the comparison folder.

### ODT B field

Configured currents remain unchanged in the data arrays. Display values are transformed only for plotting:

```text
Bz = 1.72j - (-0.35)
Bx = 1.48j - 0.28
```

The plotted units are gauss.

## Plotting conventions

The shared style encodes:

- isotope by color;
- load configuration by marker shape.

```text
DDM → square
DIS → upward triangle
DCS → circle
SCS → diamond
```

Lines are drawn first. Markers and error bars are drawn afterward.

Standard outputs use names such as:

```text
[MOT.loading.421].[lin].[t-balanced].png
[CMOT.lifetime].[fit.log].[n-balanced].svg
[ODT.BField].[log].[z].png
```

The three bracketed fields are:

```text
test tag
style tag
subvariant tag
```

The table generators rely on this naming convention.

## Aggregate output

Figure data are exported to Excel alongside the plotted SVG/PNG outputs. Each
pair-specific analysis folder contains one workbook per figure family, with a
sheet for each processing tag and balance panel; logarithmic and linear views
of the same data share a sheet. Sheets retain the configured main-variable
values, plotted x values, means, sample standard deviations, valid repetition
counts, and contributing `MMDD runXX` sources. Empty/masked values are written
as `NaN`.

The isotope-pair comparison folder has one workbook per OneNote comparison
column (`MOT loading 421`, `MOT loading 626`, `CMOT lifetime`, `MOT lifetime`,
and `ODT`). Sheets correspond to each balance/configuration/ratio figure.
Ratio sheets keep numerator and denominator source names and repetition counts
in separate columns. Fit comparison sheets use the number of time points in
the fit as `n_rep`.

[`make_multi_dual_mot_table.jl`](../helpers/make_multi_dual_mot_table.jl) collects individual SVGs into the combined table and can render the combined PNG.

[`make_multi_dual_mot_onenote.ps1`](../helpers/make_multi_dual_mot_onenote.ps1) embeds the individual PNG files into a OneNote section. It also parses `config.yaml` to make compact source captions such as:

```text
0712-run08, 10, 11
```

Each OneNote invocation adds a new page to the existing `.one` section.

The current [`multidualmot_runner.jl`](multidualmot_runner.jl) ends with:

```julia
main_multi_dual_mot_table(["--formats=svg,png,one"])
```

Therefore a complete run now also creates a new OneNote page.

The separate four-panel Result figure is generated by [`draw_selected_results.jl`](draw_selected_results.jl) and is not part of the main runner.
