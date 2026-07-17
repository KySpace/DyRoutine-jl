using Random
using Statistics
using CairoMakie

Random.seed!(1234)

# ------------------------------------------------------------
# Fourier grid
# ------------------------------------------------------------

Nk = 101
kmax = 4.0

kx = range(0, kmax, length=Nk)
ky = range(-kmax/2, kmax/2, length=51)

# Matrices whose first dimension is kx and second is ky
KX = reshape(collect(kx), :, 1)
KY = reshape(collect(ky), 1, :)

Δkx = step(kx)
Δky = step(ky)
ΔA  = Δkx * Δky

# ------------------------------------------------------------
# Spectrum parameters
# ------------------------------------------------------------

k0 = (2.0, 0.0)        # mean fringe peak position

σmain   = 0.20         # central Fourier peak width
σfringe = 0.10         # width of each individual fringe peak

Amain   = 1.0
Afringe = 0.08

Nsamples = 3000

# Fixed integration regions
r_main   = 0.80
r_fringe = 1.2

main_mask =
    KX.^2 .+ KY.^2 .<= r_main^2

fringe_mask =
    (KX .- k0[1]).^2 .+
    (KY .- k0[2]).^2 .<= r_fringe^2

# ------------------------------------------------------------
# Gaussian peak
# ------------------------------------------------------------

gaussian2d(KX, KY, cx, cy, σ) =
    exp.(-((KX .- cx).^2 .+ (KY .- cy).^2) ./ (2σ^2))

main_peak =
    Amain .* gaussian2d(KX, KY, 0.0, 0.0, σmain)

# ------------------------------------------------------------
# Generate coherent average directly
#
# σf:     wavevector-position spread
# σphase: phase spread in radians
# ------------------------------------------------------------

function coherent_average(
    KX,
    KY,
    σf,
    σphase;
    Nsamples=300,
    k0=(2.0, 0.6),
    σfringe=0.10,
    Afringe=0.08,
    main_peak
)
    Fcoh = zeros(ComplexF64, size(KX, 1), size(KY, 2))

    for _ in 1:Nsamples
        # Isotropic 2D wavevector jitter
        δkx = σf * randn()
        δky = σf * randn()

        kxi = k0[1] + δkx
        kyi = k0[2] + δky

        ϕ = σphase * randn()

        fringe =
            Afringe .*
            gaussian2d(KX, KY, kxi, kyi, σfringe) .*
            cis(ϕ)

        Fcoh .+= main_peak .+ fringe
    end

    return Fcoh ./ Nsamples
end

# ------------------------------------------------------------
# Coefficients normalized by the central peak
# ------------------------------------------------------------

function coherence_coefficients(
    Fcoh,
    main_mask,
    fringe_mask;
    ΔA=1.0
)
    # ΔA cancels in these ratios for a uniform grid,
    # but is included to emphasize that these are integrals.
    fringe_L1 = ΔA * sum(abs,  @view Fcoh[fringe_mask])
    main_L1   = ΔA * sum(abs,  @view Fcoh[main_mask])

    fringe_L2 = ΔA * sum(abs2, @view Fcoh[fringe_mask])
    main_L2   = ΔA * sum(abs2, @view Fcoh[main_mask])

    C1 = fringe_L1 / (main_L1 + fringe_L1)
    C2 = fringe_L2 / (main_L2 + fringe_L2)

    return (; C1, C2)
end

# ------------------------------------------------------------
# Parameter scan
# ------------------------------------------------------------

σf_values =
    range(0.0, 0.7, length=31)

σphase_values =
    range(0.0, 2pi, length=31)

C1map =
    zeros(length(σf_values), length(σphase_values))

C2map =
    similar(C1map)

for (ifreq, σf) in enumerate(σf_values)
    for (iphase, σphase) in enumerate(σphase_values)

        Fcoh = coherent_average(
            KX,
            KY,
            σf,
            σphase;
            Nsamples,
            k0,
            σfringe,
            Afringe,
            main_peak
        )

        C = coherence_coefficients(
            Fcoh,
            main_mask,
            fringe_mask;
            ΔA
        )

        C1map[ifreq, iphase] = C.C1
        C2map[ifreq, iphase] = C.C2
    end
end

# Normalize each estimator to its value at σf = σphase = 0.
# This makes their relative degradation directly comparable.

C1relative = C1map ./ C1map[1, 1]
C2relative = C2map ./ C2map[1, 1]

##
# ------------------------------------------------------------
# Plot parameter maps
# ------------------------------------------------------------

fig = Figure(size = (1050, 440))

ax1 = Axis(
    fig[1, 1],
    xlabel = "phase spread sigma_phi",
    ylabel = "frequency spread sigma_f",
    title = "Sum of absolute coherent amplitude"
)

hm1 = heatmap!(
    ax1,
    σphase_values,
    σf_values,
    C1relative',
    colorrange = (0, 0.1)
)

Colorbar(
    fig[1, 2],
    hm1,
    label = "relative coefficient"
)

ax2 = Axis(
    fig[1, 3],
    xlabel = "phase spread sigma_phi",
    ylabel = "frequency spread sigma_f",
    title = "Sum of squared coherent amplitude"
)

hm2 = heatmap!(
    ax2,
    σphase_values,
    σf_values,
    C2relative',
    colorrange = (0, 0.005)
)

Colorbar(
    fig[1, 4],
    hm2,
    label = "relative coefficient",
    
)

path_output = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\AnlzRoutine\50.MeanAbsl2D.VaryMask.Squared"
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
fig |> f -> save(joinpath(path_output, "square_or_amplitude.pdf"), f)