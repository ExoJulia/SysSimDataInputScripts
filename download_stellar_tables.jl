if !isdefined(Main,:parsed_args)
  include("parse_args.jl")
end

begin
    nothing_to_download = true

if !isfile("q1_q17_dr24_stellar.csv") && !isfile("inputs/q1_q17_dr24_stellar.csv")
   println("# Downloading q1_q17_dr24_stellar.csv")
   input_dir = parsed_args["input-path"] != nothing ? parsed_args["input-path"] : "inputs"
   output_dir = parsed_args["output-path"] != nothing ? parsed_args["output-path"] : input_dir
   download("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=q1_q17_dr24_stellar&format=csv&select=*",joinpath(output_dir,"q1_q17_dr24_stellar.csv"))
   nothing_to_download = false
end

if !isfile("q1_q17_dr25_stellar.csv") && !isfile("inputs/q1_q17_dr25_stellar.csv")
   println("# Downloading q1_q17_dr25_stellar.csv")
   input_dir = parsed_args["input-path"] != nothing ? parsed_args["input-path"] : "inputs"
   output_dir = parsed_args["output-path"] != nothing ? parsed_args["output-path"] : input_dir
   download("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=q1_q17_dr25_stellar&format=csv&select=*",joinpath(output_dir,"q1_q17_dr25_stellar.csv"))
   nothing_to_download = false
end

if !isfile("dr25fgk_relaxcut_osds.jld2")
   println("# Downloading dr25fgk_relaxcut_osds.jld2")
   output_dir = parsed_args["output-path"] != nothing ? parsed_args["output-path"] : "."
   download("https://scholarsphere.psu.edu/resources/460cb19b-86e9-4c04-ac17-69950444437f/downloads/13899",joinpath(output_dir,"dr25fgk_relaxcut_osds.jld2"))
   nothing_to_download = false
end

if !isfile("allosds.jld") && parsed_args["download-large-files"]
   println("# Downloading allosds.jld (this is a large file)")
   output_dir = parsed_args["output-path"] != nothing ? parsed_args["output-path"] : "."
   download("https://scholarsphere.psu.edu/resources/ef562cd1-7a32-40fd-b59b-36f06c47db1c/downloads/13890",joinpath(output_dir,"allosds.jld"))
   nothing_to_download = false
end
if nothing_to_download
    println("# All stellar files that need to be downloaded are already on disk.")
end

end # begin
