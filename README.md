# Instructions

The SysSimDataInputScripts repository should be automatically installed as a submodule of SysSimData (https://github.com/ExoJulia/SysSimData) in a directory named inputs.
Run scripts/makes.jl from the SysSimData directory (typically named data).

To generate the binary file KeplerMAST_TargetProperties.jld2 from inputs/KeplerMAST_TargetProperties.csv, run
'julia scripts/make.jl --make-target-properties'.

To download the CSV version of the stellar catalogs (useful for remaking binary files yourself), you'll want to run
'julia scripts/make.jl --download-stellar-catalog'.

## Other scripts
There are some other scripts that aren't yet fully automated, but still may be useful:
- download_stellar_tables.jl: Downloads CSV stellar tables (can be called by 'scripts/make.jl')
- make_target_properties.jl: Makes 'KeplerMAST_TargetProperties.jld2' from 'inputs/KeplerMAST_TargetProperties.csv' (called by 'scripts/make.jl')
- combine_dr25_gaia.jl: For combining 'inputs/q1_q17_dr25_stellar.csv', 'gaiadr2_keplerdr25_crossref.csv', and 'KeplerMAST_TargetProperties.csv' into 'q1q17_dr25_gaia_fgk.jld'.
- convert_mcmc_quant_to_jld2.jl:  For updating planet properties in 'inputs/q1_q17_dr25_koi.csv' with quantiles from MCMC posteriors stored in 'dr25_koi_mcmc_quant.csv'.
- make_small_osds.jl:  Used to make 'dr25fgk_small_osds.jld2' from 'allosds.jld' for testing on laptops.
- parse_args.jl: Used by multiple scripts to parse command line arguments.  Use 'julia scripts/make.jl -h' to see the command-line options.

## Downloading files manually
If you would rather download and copy files manually, key files can be found in a Penn State's ScholarSphere [Collection of SysSim Data Files](https://scholarsphere.psu.edu/resources/175d7036-555b-4fc9-b57e-db8af8bb4004a).

### One sigma depth functions
- [FGK Kepler target stars selected for Hsu et al. 2019 and He et al. 2020 studies](https://scholarsphere.psu.edu/resources/460cb19b-86e9-4c04-ac17-69950444437f/downloads/13899)
- [All Kepler target stars](https://scholarsphere.psu.edu/resources/ef562cd1-7a32-40fd-b59b-36f06c47db1c)

### Stellar tables (in case NExScI changes things)
- https://psu.box.com/v/KeplerStarTableDR24csv
- https://psu.box.com/v/KeplerStarTableDR25csv


