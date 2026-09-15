using CairoMakie
using Printf: @sprintf, @printf

using XLSX
using LsqFit: curve_fit, coef

function read_cmot_decay(path::AbstractString)
    samples = Tuple{Float64, Float64}[]
    XLSX.openxlsx(path) do workbook
        for sheet_name in XLSX.sheetnames(workbook)
            data = workbook[sheet_name][:]
            size(data, 2) >= 5 || throw(DimensionMismatch("$path, sheet $sheet_name: expected at least five columns"))
            for row in axes(data, 1)
                t, num = data[row, 2], data[row, 5]
                # Skip the text headers and empty rows, retaining every loop measurement.
                t isa Real && num isa Real || continue
                isfinite(t) && isfinite(num) && t >= 0 ||
                    throw(ArgumentError("$path, sheet $sheet_name, row $row: invalid (t, N) = ($t, $num)"))
                push!(samples, (Float64(t), Float64(num)))
            end
        end
    end
    isempty(samples) && throw(ArgumentError("No numeric (t, N) rows found in $path"))
    samples
end

function model_cmot_decay(t::AbstractVector, params::AbstractVector)
    num_initial, tau_decay, background = params
    @. num_initial * exp(-t / tau_decay) + background
end

function fit_cmot_decay(samples::AbstractVector{<:Tuple{Real, Real}}; scale_t::Real=1e-3)
    isfinite(scale_t) && scale_t > 0 || throw(ArgumentError("scale_t must be finite and positive"))
    t = first.(samples) .* scale_t
    num = last.(samples)
    length(unique(t)) >= 4 || throw(ArgumentError("Need at least four distinct hold times for a three-parameter fit"))
    scale_num = maximum(abs, num)
    scale_num > 0 || throw(ArgumentError("Cannot fit an all-zero decay"))
    num_scaled = num ./ scale_num
    span_t = maximum(t) - minimum(t)
    params_initial = [maximum(num_scaled) - minimum(num_scaled), span_t / 4, minimum(num_scaled)]
    fit = curve_fit(model_cmot_decay, t, num_scaled, params_initial;
        lower=[0.0, eps(Float64), -Inf])
    params = coef(fit) .* [scale_num, 1.0, scale_num]
    (; params, num_initial=params[1], tau_decay=params[2], background=params[3])
end


path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\CMOT lifetime"
val_case = [:SW162, :SW164, :MX162, :MX164]
scale_t = 1e-3 # Workbook hold time in ms; fits and plots use seconds.
scale_num_plot = 1e7

# Each dictionary entry is a vector of (t_hold_ms, N), pooled over all sheets.
decays = Dict(label => read_cmot_decay(joinpath(path_root, "CMOT lifetime $(label).xlsx")) for label in val_case)
decay_fits = Dict(label => fit_cmot_decay(decays[label]; scale_t) for label in val_case)
num_initial = Dict(label => decay_fits[label].num_initial for label in val_case)
tau_decay = Dict(label => decay_fits[label].tau_decay for label in val_case)
background = Dict(label => decay_fits[label].background for label in val_case)

# Rerun from here to replot the already loaded data and fits.
fig_decay = Figure(size=(1000, 720))
clr_type = Dict(:MX => :darkviolet, :SW => :darkblue)
num_max_plot = maximum(maximum(last, decay) for decay in values(decays)) / scale_num_plot
for (idx_row, isotope) in enumerate((162, 164)), (idx_col, type) in enumerate((:MX, :SW))
    label = Symbol(type, isotope)
    decay = decays[label]
    fit = decay_fits[label]
    t = first.(decay) .* scale_t
    t_plot = range(minimum(t), maximum(t); length=400)
    ax = Axis(fig_decay[idx_row, idx_col];
        title=string(label), xlabel="t_hold (sec)", ylabel="Atom number (×10⁷)",
        limits=(nothing, (-0.02 * num_max_plot, 1.22 * num_max_plot)))
    scatter!(ax, t, last.(decay) ./ scale_num_plot; color=(clr_type[type], 0.55), markersize=7)
    lines!(ax, t_plot, model_cmot_decay(t_plot, fit.params) ./ scale_num_plot;
        color=clr_type[type], linewidth=2)
    text!(ax, 0.97, 0.97; space=:relative, align=(:right, :top), fontsize=13,
        text=@sprintf("τ_decay = %.2f sec\nnum_initial = %.3g\nbackground = %.3g", fit.tau_decay, fit.num_initial, fit.background))
    @printf("%s: %d samples, tau = %.4f sec, N0 = %.6g, background = %.6g\n",
        string(label), length(decay), fit.tau_decay, fit.num_initial, fit.background)
end

path_svg = joinpath(path_root, "anlz_cmot_lifetime_dual.svg")
save(path_svg, fig_decay)
fig_decay
