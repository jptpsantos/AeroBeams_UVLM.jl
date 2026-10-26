# Filesystem helpers only; safe to load without a solver or plotting backend.
module PazyLCOOutput

export pazy_lco_angle_directory, pazy_lco_new_run_directory,
    pazy_lco_latest_output_directory, pazy_lco_prepare_output_directory

"""Directory for a root incidence specified in degrees (e.g. AoA3 or AoA3p5)."""
function pazy_lco_angle_directory(base_directory, angle_of_attack_deg)
    isfinite(angle_of_attack_deg) || throw(ArgumentError("Root incidence must be finite"))
    # Remove only conversion roundoff, e.g. rad2deg(deg2rad(3)) = 3.0000000000000004.
    degrees = round(Float64(angle_of_attack_deg); digits=10)
    iszero(degrees) && (degrees = 0.0)
    label = isinteger(degrees) ? string(round(Int, degrees)) : string(degrees)
    label = replace(label, "-" => "m", "." => "p")
    return joinpath(base_directory, "AoA$label")
end

"""Atomically reserve a fresh numbered run; never reuse an existing directory."""
function pazy_lco_new_run_directory(base_directory, angle_of_attack_deg)
    angle_directory = pazy_lco_angle_directory(base_directory, angle_of_attack_deg)
    mkpath(angle_directory)
    index = 1
    for name in readdir(angle_directory)
        matched = match(r"^run_(\d+)$", name)
        isnothing(matched) || (index = max(index, parse(Int, matched.captures[1])+1))
    end
    while true
        directory = joinpath(angle_directory, "run_" * lpad(string(index), 3, '0'))
        if !ispath(directory)
            try
                mkdir(directory)
                return directory
            catch
                # Another process can reserve the same number between the check and mkdir.
                ispath(directory) || rethrow()
            end
        end
        index += 1
    end
end

"""Use an explicit empty directory or allocate a new run under the angle folder."""
function pazy_lco_prepare_output_directory(output_directory;
    base_directory, angle_of_attack_deg)
    isnothing(output_directory) &&
        return pazy_lco_new_run_directory(base_directory, angle_of_attack_deg)
    if ispath(output_directory)
        isdir(output_directory) && isempty(readdir(output_directory)) ||
            throw(ArgumentError("Refusing to overwrite existing LCO output: $output_directory. " *
                "Use an empty directory, or omit output_directory to create a fresh angle/run folder."))
    else
        mkpath(output_directory)
    end
    return output_directory
end

"""Find the newest numbered run, then legacy angle/flat data, without writing.

Legacy flat results are accepted only when their saved incidence matches the
requested angle. Thus requesting 3 or 7 degrees never falls back to 5-degree data.
The newest reserved run is selected even if incomplete, so missing branch data
is reported by the plotter rather than silently replaced with an older run.
"""
function pazy_lco_latest_output_directory(base_directory, angle_of_attack_deg)
    angle_directory = pazy_lco_angle_directory(base_directory, angle_of_attack_deg)
    if isdir(angle_directory)
        runs = Tuple{Int,String}[]
        for name in readdir(angle_directory)
            matched = match(r"^run_(\d+)$", name)
            directory = joinpath(angle_directory, name)
            if !isnothing(matched) && isdir(directory)
                push!(runs, (parse(Int, matched.captures[1]), directory))
            end
        end
        isempty(runs) || return last(sort!(runs; by=first))[2]
    end
    for directory in (angle_directory, base_directory)
        settings = joinpath(directory, "pazy_linear_LCO_settings.txt")
        isfile(settings) || continue
        matched = match(r"Root incidence \[deg\]:\s*([-+0-9.eE]+)", read(settings, String))
        isnothing(matched) && continue
        saved_angle = parse(Float64, matched.captures[1])
        isapprox(saved_angle, angle_of_attack_deg; atol=1e-9, rtol=0) && return directory
    end
    throw(ArgumentError("No saved LCO run found for $angle_of_attack_deg degrees in $base_directory"))
end

end # module
