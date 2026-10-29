using JLD2, CSV, DataFrames

# Export saved parameters only; do not refit the GMM or re-call CNVs.
length(ARGS) == 2 || error("Usage: julia --project=REPO export_public_copy_number_peaks.jl CN_PARAM_JLD2 OUTPUT_CSV")
input_file, output_file = ARGS
parameters = load(input_file, "cn_param")
rows = NamedTuple[]
for p in parameters
    println("sample=", p.sample_name, " fields=", propertynames(p))
    n = length(p.copy_labels)
    all(length(x)==n for x in [p.means, p.stds, p.weights]) || error("Parameter lengths differ")
    for i in eachindex(p.copy_labels)
        push!(rows, (sample=String(p.sample_name), cn=Int(p.copy_labels[i]),
            mean=Float64(p.means[i]), std=Float64(p.stds[i]), weight=Float64(p.weights[i])))
    end
end
mkpath(dirname(abspath(output_file)))
CSV.write(output_file, DataFrame(rows))
println("Exported ", length(rows), " saved peak rows to ", output_file)
