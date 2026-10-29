using JLD2
using CSV
using DataFrames

length(ARGS) == 2 || error("Usage: julia --project=REPO SCRIPT CN_PARAM_JLD2 OUTPUT_CSV")
const INPUT_FILE = ARGS[1]
const OUTPUT_FILE = ARGS[2]
const OUTPUT_DIR = dirname(abspath(OUTPUT_FILE))

const SAMPLE_LABELS = Dict(
    "G6" => "U1",
    "G6A7first" => "U1.1",
    "G6B11first" => "U1.2",
    "G6C7first" => "U1.3",
    "G6A7E1" => "U1.1.1",
    "G6A7E11" => "U1.1.2",
    "G6B11C1" => "U1.2.1",
    "G6B11C6" => "U1.2.2",
    "G6C7E1" => "U1.3.1",
    "G6C7E9" => "U1.3.2",
)

isfile(INPUT_FILE) || error("Missing JLD2 parameter file: $INPUT_FILE")
mkpath(OUTPUT_DIR)
parameters = load(INPUT_FILE, "cn_param")

rows = NamedTuple[]
for parameter in parameters
    raw_name = String(parameter.sample_name)
    haskey(SAMPLE_LABELS, raw_name) || error("No U1 label mapping for $raw_name")
    lengths = (length(parameter.copy_labels), length(parameter.means),
               length(parameter.stds), length(parameter.weights))
    length(unique(lengths)) == 1 || error("Unequal GMM parameter lengths for $raw_name")
    for index in eachindex(parameter.means)
        push!(rows, (
            raw_sample_name = raw_name,
            sample_name = SAMPLE_LABELS[raw_name],
            copy_number = Int(parameter.copy_labels[index]),
            mean = Float64(parameter.means[index]),
            std = Float64(parameter.stds[index]),
            weight = Float64(parameter.weights[index]),
        ))
    end
end

table = DataFrame(rows)
CSV.write(OUTPUT_FILE, table)
println("Wrote ", nrow(table), " peak rows to: ", OUTPUT_FILE)

