using HDF5, DataFrames, CSV, JLD
using StatsBase, Polynomials, CurveFit
using ExoplanetsSysSim

kep_filename = joinpath(dirname(pathof(ExoplanetsSysSim)), "../data/inputs/", "q1_q17_dr25_stellar.csv")
gaia_filename = joinpath(dirname(pathof(ExoplanetsSysSim)), "../data", "gaiadr2_keplerdr25_crossref.csv")
mast_filename = joinpath(dirname(pathof(ExoplanetsSysSim)), "../data", "KeplerMAST_TargetProperties.csv")
stellar_catalog_file_out = joinpath(dirname(pathof(ExoplanetsSysSim)), "../data", "q1q17_dr25_gaia_fgk_interpolate_ebprp.jld2")

kep_df = CSV.read(kep_filename)
gaia_df = CSV.read(gaia_filename)
mast_df = CSV.read(mast_filename)

dup_gaiaid = findall(nonunique(DataFrame(x = gaia_df[:source_id])))
gaiaid_keep = trues(size(gaia_df,1))
gaiaid_keep[dup_gaiaid] .= false
gaia_df = gaia_df[gaiaid_keep,:]

println("Total crossref target stars = ", length(gaia_df[:kepid]))

mag_diff = gaia_df[:phot_g_mean_mag].-gaia_df[:kepmag]
quant_arr = quantile(mag_diff, [0.067,0.933])   # 1.5-sigma cut
mag_match = findall(x->quant_arr[1]<=x<=quant_arr[2], mag_diff)
gaia_df = gaia_df[mag_match,:]

gaia_col = [:kepid,:source_id,:parallax,:parallax_error,:astrometric_gof_al,:astrometric_excess_noise_sig,:bp_rp,:priam_flags,:teff_val,:teff_percentile_lower,:teff_percentile_upper,:radius_val,:radius_percentile_lower,:radius_percentile_upper,:lum_val,:lum_percentile_lower,:lum_percentile_upper,:e_bp_min_rp_val,:e_bp_min_rp_percentile_lower,:e_bp_min_rp_percentile_upper, :l, :b]
df = join(kep_df, gaia_df[:,gaia_col], on=:kepid)
df = join(df, mast_df, on=:kepid)
kep_df = nothing
gaia_df = nothing

println("Total target stars (KOIs) matching magnitude = ", length(df[:kepid]), " (", sum(df[:nkoi]),")")

df[:teff] = df[:teff_val]
df[:teff_err1] = df[:teff_percentile_upper].-df[:teff_val]
df[:teff_err2] = df[:teff_percentile_lower].-df[:teff_val]
deletecols!(df, :teff_val)
deletecols!(df, :teff_percentile_upper)
deletecols!(df, :teff_percentile_lower)
df[:radius_err1] = df[:radius_percentile_upper].-df[:radius_val]
df[:radius_err2] = df[:radius_percentile_lower].-df[:radius_val]
df[:radius] = df[:radius_val]
deletecols!(df, :radius_val)
deletecols!(df, :radius_percentile_upper)
deletecols!(df, :radius_percentile_lower)

not_binary_suspect = (df[:astrometric_gof_al] .<= 20) .& (df[:astrometric_excess_noise_sig] .<= 5)
astrometry_good = []
for x in 1:length(df[:kepid])
    if !(ismissing(df[x,:priam_flags]))
        pflag = string(df[x,:priam_flags])
         if (pflag[2] == '0') & (pflag[3] == '0') # WARNING: Assumes flag had first '0' removed by crossref script
             push!(astrometry_good, true)
         else
             push!(astrometry_good, false)
         end
     else
         push!(astrometry_good, false)
     end
end
astrometry_good = astrometry_good .& (df[:parallax_error] .< 0.1*df[:parallax])
# planet_search = df[:kepmag] .<= 16.

has_mass = .! (ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
has_radius = .! (ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))
#has_dens = .! (ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))

has_cdpp = .! (ismissing.(df[:rrmscdpp01p5]) .| ismissing.(df[:rrmscdpp02p0]) .| ismissing.(df[:rrmscdpp02p5]) .| ismissing.(df[:rrmscdpp03p0]) .| ismissing.(df[:rrmscdpp03p5]) .| ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:rrmscdpp05p0]) .| ismissing.(df[:rrmscdpp06p0]) .| ismissing.(df[:rrmscdpp07p5]) .| ismissing.(df[:rrmscdpp09p0]) .| ismissing.(df[:rrmscdpp10p5]) .| ismissing.(df[:rrmscdpp12p0]) .| ismissing.(df[:rrmscdpp12p5]) .| ismissing.(df[:rrmscdpp15p0]))
has_rest = .! (ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))
has_limbd = .! (ismissing.(df[:limbdark_coeff1]) .| ismissing.(df[:limbdark_coeff2]) .| ismissing.(df[:limbdark_coeff3]) .| ismissing.(df[:limbdark_coeff4]))

mast_cut =.&(df[:numLCEXqtrs].>0,df[:numLCqtrs].>4)

# is_usable = .&(has_radius, has_mass, has_rest, has_cdpp, has_dens, astro_good, not_binary, planet_search, mast_cut)
is_usable = .&(has_radius, has_mass, has_rest, has_cdpp, mast_cut, astrometry_good, not_binary_suspect)

df = df[findall(is_usable),:]
println("Total stars (KOIs) with valid parameters = ", length(df[:kepid]), " (", sum(df[:nkoi]),")")



##### To make a new copy of 'df' for the purposes of fitting an FGK main sequence (using uncorrected 'bp_rp'), and then fitting+interpolating 'e_bp_min_rp_val', so that we can make new cuts based on corrected 'bp_rp' (using interpolated 'e_bp_min_rp_val'):

df_fit = deepcopy(df)
has_ebprp = .!isnan.(df_fit[:e_bp_min_rp_val])
df_fit = df_fit[findall(has_ebprp),:]

#fgk_color = (0.5 .<= df_fit[:bp_rp] .- df_fit[:e_bp_min_rp_val] .<= 1.7)
fgk_color = (0.5 .<= df_fit[:bp_rp] .<= 1.7)
# FGK = findall(fgk_color)
df_fit = df_fit[findall(fgk_color),:]
near_FGK_MS = []
coeff = [2.5,-3.6,0.9]
for i in 1:6
    #global near_FGK_MS = (log10.(df_fit[:lum_val]).< map(x->polyval(Poly(coeff),x),df_fit[:bp_rp] .-  df_fit[:e_bp_min_rp_val]) .+ log10(1.75))
    #global coeff = poly_fit(df_fit[near_FGK_MS,:bp_rp] .- df_fit[near_FGK_MS,:e_bp_min_rp_val],log10.(df_fit[near_FGK_MS,:lum_val]),2)
    global near_FGK_MS = (log10.(df_fit[:lum_val]).< map(x->polyval(Poly(coeff),x),df_fit[:bp_rp]) .+ log10(1.75))
    global coeff = poly_fit(df_fit[near_FGK_MS,:bp_rp], log10.(df_fit[near_FGK_MS,:lum_val]),2)
    println(i, ": ",coeff," : ", length(findall(near_FGK_MS)))
end
FGK = findall(near_FGK_MS)

println("Total FGK stars (KOIs) with valid parameters = ", length(FGK), " (", sum(df_fit[FGK, :nkoi]),")")

df_fit = df_fit[FGK,:]


##### To interpolate 'e_bp_min_rp_val':

using Statistics
using Dierckx

logL = log10.(df_fit[:lum_val])
teff = df_fit[:teff]
bprp = df_fit[:bp_rp]
ebprp = df_fit[:e_bp_min_rp_val]
dist = 1 ./df_fit[:parallax]
cos_b = cos.(deg2rad.(df_fit[:b]))

ebprp_norm_d = ebprp ./ dist

# To bin 'bprp', compute the median 'ebprp_norm_d' in each bin, and interpolate 'ebprp_norm_d' as a function of 'bprp':
bins = 20
bin_sizes = Int(ceil(size(df_fit,1)/bins))
idx_sort_bprp = sortperm(bprp)

bprp_med_bins = zeros(bins)
ebprp_norm_d_med_bins = zeros(bins)
for i in 1:bins
    i_start, i_stop = (i-1)*bin_sizes+1, min(i*bin_sizes+1, size(df_fit,1))
    idx = idx_sort_bprp[i_start:i_stop]
    bprp_med_bins[i] = median(bprp[idx])
    ebprp_norm_d_med_bins[i] = median(ebprp_norm_d[idx])
end

spl = Spline1D(bprp_med_bins, ebprp_norm_d_med_bins)
ebprp_norm_d_interp = spl(bprp) # interpolated values of 'ebprp_norm_d'
ebprp_interp = ebprp_norm_d_interp .* dist

# Insert the interpolated values of E(Bp-Rp) ('ebprp_interp_refit') as a column before the column of real E(Bp-Rp) ('e_bp_min_rp_val'):
insert!(df_fit, (1:length(names(df_fit)))[names(df_fit) .== :e_bp_min_rp_val][1], ebprp_interp, :e_bp_min_rp_interp)


##### To make some plots checking the stellar parameters:
#=
using PyPlot

ids = range(1, stop=length(bprp), step=10)

# Color histograms:
hist([bprp, bprp .- ebprp, bprp .- ebprp_interp], bins=100, histtype="step", label=["Bp-Rp", "Bp-Rp - E(Bp-Rp)", "Bp-Rp - E*(Bp-Rp)"])
xlabel("Bp-Rp", fontsize=20); ylabel("Counts", fontsize=20)
legend(fontsize=20)

# Extinction histograms:
hist([ebprp, ebprp_interp], bins=100, histtype="step", label=["E(Bp-Rp)", "E*(Bp-Rp)"])
xlim([0, 0.8])
xlabel("E(Bp-Rp)", fontsize=20); ylabel("Counts", fontsize=20)
legend(fontsize=20)

# Color-magnitude (HR) diagram:
PyPlot.scatter(bprp[ids], logL[ids], s=1, label="Bp-Rp (uncorrected)")
PyPlot.scatter((bprp .- ebprp)[ids], logL[ids], s=1, label="Bp-Rp - E(Bp-Rp)")
PyPlot.scatter((bprp .- ebprp_interp)[ids], logL[ids], s=1, label="Bp-Rp - E*(Bp-Rp)")
xlabel("Bp-Rp", fontsize=20); ylabel("log(L)", fontsize=20)
legend(fontsize=20)

# Color-extinction diagram:
PyPlot.scatter(bprp, ebprp, s=1, label="E(Bp-Rp)")
PyPlot.scatter(bprp, ebprp_norm_d, s=1, label="E(Bp-Rp)/dist")
PyPlot.scatter(bprp_med_bins, ebprp_norm_d_med_bins, marker="x", s=100, label="Median per bin")
PyPlot.scatter(bprp, ebprp_norm_d_interp, s=0.1, label="E*(Bp-Rp)/dist")
xlim([0.5, 1.7])
ylim([0, 1])
xlabel("Bp-Rp", fontsize=20); ylabel("E(Bp-Rp)", fontsize=20)
legend(fontsize=20)

# Color-temperature diagram:
PyPlot.scatter(bprp[ids], teff[ids], s=1, label="Bp-Rp (uncorrected)")
PyPlot.scatter((bprp .- ebprp)[ids], teff[ids], s=1, label="Bp-Rp - E(Bp-Rp)")
PyPlot.scatter((bprp .- ebprp_interp)[ids], teff[ids], s=1, label="Bp-Rp - E*(Bp-Rp)")
xlabel("Bp-Rp", fontsize=20); ylabel("T_eff (K)", fontsize=20)
legend(fontsize=20)
=#





##### To re-do cuts and FGK main sequence fitting (on 'df') using 'bp_rp'-'e_bp_min_rp_val' (interpolated):

ebprp_interp = spl(df[:bp_rp]) ./ df[:parallax]
fgk_color = (0.5 .<= df[:bp_rp] .- ebprp_interp .<= 1.7)
df = df[findall(fgk_color),:]
near_FGK_MS = []
coeff = [2.5,-3.6,0.9]
for i in 1:6
    global ebprp_interp = spl(df[:bp_rp]) ./df[:parallax]
    global near_FGK_MS = (log10.(df[:lum_val]) .< map(x->polyval(Poly(coeff),x),df[:bp_rp] .- ebprp_interp) .+ log10(1.75))
    global ebprp_interp = spl(df[near_FGK_MS,:bp_rp]) ./df[near_FGK_MS,:parallax]
    global coeff = poly_fit(df[near_FGK_MS,:bp_rp] .- ebprp_interp, log10.(df[near_FGK_MS,:lum_val]),2)
    println(i, ": ",coeff," : ", length(findall(near_FGK_MS)))
end
FGK = findall(near_FGK_MS)

println("Total FGK stars (KOIs) after re-fitting using (interpolated) corrected colors = ", length(FGK), " (", sum(df[FGK, :nkoi]),")")

df = df[FGK,:]



logL_refit = log10.(df[:lum_val])
teff_refit = df[:teff]
bprp_refit = df[:bp_rp]
dist_refit = 1 ./df[:parallax]
cos_b_refit = cos.(deg2rad.(df[:b]))

ebprp_refit = df[:e_bp_min_rp_val] # NOTE: will have NaNs!
ebprp_norm_d_interp_refit = spl(bprp_refit)
ebprp_interp_refit = ebprp_norm_d_interp_refit .* dist_refit

# Insert the interpolated values of E(Bp-Rp) ('ebprp_interp_refit') as a column before the column of real E(Bp-Rp) ('e_bp_min_rp_val'):
insert!(df, (1:length(names(df)))[names(df) .== :e_bp_min_rp_val][1], ebprp_interp_refit, :e_bp_min_rp_interp)


##### To make some plots checking the stellar parameters:
#=
ids = range(1, stop=length(bprp_refit), step=10)

# Color histograms:
hist([bprp_refit, bprp_refit .- ebprp_refit, bprp_refit .- ebprp_interp_refit], bins=100, histtype="step", label=["Bp-Rp", "Bp-Rp - E(Bp-Rp)", "Bp-Rp - E*(Bp-Rp)"])
xlabel("Bp-Rp", fontsize=20); ylabel("Counts", fontsize=20)
legend(fontsize=20)

# Extinction histograms:
hist([ebprp_refit, ebprp_interp_refit], bins=100, histtype="step", label=["E(Bp-Rp)", "E*(Bp-Rp)"])
xlim([0, 0.8])
xlabel("E(Bp-Rp)", fontsize=20); ylabel("Counts", fontsize=20)
legend(fontsize=20)

# Color-magnitude (HR) diagram:
PyPlot.scatter(bprp_refit[ids], logL_refit[ids], s=1, label="Bp-Rp (uncorrected)")
PyPlot.scatter((bprp_refit .- ebprp_refit)[ids], logL_refit[ids], s=1, label="Bp-Rp - E(Bp-Rp)")
PyPlot.scatter((bprp_refit .- ebprp_interp_refit)[ids], logL_refit[ids], s=1, label="Bp-Rp - E*(Bp-Rp)")
xlabel("Bp-Rp", fontsize=20); ylabel("log(L)", fontsize=20)
legend(fontsize=20)

# Color-extinction diagram:
PyPlot.scatter(bprp_refit, ebprp_refit, s=1, label="E(Bp-Rp)")
PyPlot.scatter(bprp_refit, ebprp_refit ./ dist_refit, s=1, label="E(Bp-Rp)/dist")
PyPlot.scatter(bprp_med_bins, ebprp_norm_d_med_bins, marker="x", s=100, label="Median per bin")
PyPlot.scatter(bprp_refit, ebprp_norm_d_interp_refit, s=0.1, label="E*(Bp-Rp)/dist")
#xlim([0.5, 1.7])
ylim([0, 1])
xlabel("Bp-Rp", fontsize=20); ylabel("E(Bp-Rp)", fontsize=20)
legend(fontsize=20)

# Color-temperature diagram:
PyPlot.scatter(bprp_refit[ids], teff_refit[ids], s=1, label="Bp-Rp (uncorrected)")
PyPlot.scatter((bprp_refit .- ebprp_refit)[ids], teff_refit[ids], s=1, label="Bp-Rp - E(Bp-Rp)")
PyPlot.scatter((bprp_refit .- ebprp_interp_refit)[ids], teff_refit[ids], s=1, label="Bp-Rp - E*(Bp-Rp)")
xlabel("Bp-Rp", fontsize=20); ylabel("T_eff (K)", fontsize=20)
legend(fontsize=20)
=#





# See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
# TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
symbols_to_keep = [ :kepid, :source_id, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :teff, :bp_rp, :lum_val, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :dataspan, :dutycycle, :limbdark_coeff1, :limbdark_coeff2, :limbdark_coeff3, :limbdark_coeff4, :contam, :e_bp_min_rp_interp, :e_bp_min_rp_val, :e_bp_min_rp_percentile_lower, :e_bp_min_rp_percentile_upper]
#df = df[FGK, symbols_to_keep]
#df = df_fit[:, symbols_to_keep]
df = df[:, symbols_to_keep]
tmp_df = DataFrame()
for col in names(df)
    tmp_df[col] = collect(skipmissing(df[col]))
end
df = tmp_df

save(stellar_catalog_file_out,"stellar_catalog", df)
