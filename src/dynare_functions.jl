module DFunctions

using RuntimeGeneratedFunctions
using StatsFuns
using Suppressor
using TimeDataFrames
using Tokenize

export DynareFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

struct DynareFunctions
    dynamic!::Function
    dynamic_tmp_nbr::Vector{Int64}
    static!::Function
    static_tmp_nbr::Vector{Int64}
    set_auxiliary_variables!::Function
    set_dynamic_auxiliary_variables!::Function
    steady_state!::Function
    function DynareFunctions(compileoption, modfileinfo, modfilename, orig_maximum_lag, orig_maximum_lead)
        dynamic_tmp_nbr = zeros(4)
        static_tmp_nbr = zeros(4)
        if isfile(modfilename * "Dynamic.jl")
            dynamic_tmp_nbr .= dynamic_functions(modfilename)
            dynamic! = DFunctions.dynamic!
        else
            dynamic! = ()->nothing
        end
        if isfile(modfilename * "Static.jl")
            static_tmp_nbr .= static_functions(modfilename)
            static! = DFunctions.static!
        else
            static! = ()->nothing
        end
        if modfileinfo.has_auxiliary_variables
            set_dynamic_auxiliary_variables! =
                DFunctions.load_set_dynamic_auxiliary_variables(modfilename)
            set_auxiliary_variables! =
                load_dynare_function2(modfilename * "SetAuxiliaryVariables")
        else
            # no auxiliary variables
            set_dynamic_auxiliary_variables! = (a, b, c) -> nothing
            set_auxiliary_variables! = (a, b, c) -> nothing
        end
        if modfileinfo.has_steadystate_file
            steady_state! = load_steady_state_function(modfilename * "SteadyState2", compileoption)
        else
            steady_state! = ()->nothing
        end
        new(dynamic!,
            dynamic_tmp_nbr,
            static!,
            static_tmp_nbr,
            set_auxiliary_variables!,
            set_dynamic_auxiliary_variables!,
            steady_state!)
    end
end

function make_function(file)
    expr = Meta.parse(join(file, "\n"))
    return @RuntimeGeneratedFunction(expr)
end

function read_dynare_files!(files::Dict{String, Vector{String}},
                            tmp_nbr_str::Vector{String},
                            fname::String)
    infunction = false
    content = []
    filename = ""
    for line in readlines(fname)
        if infunction
            push!(content, line)
            if startswith(line, "end")
                !infunction && error("end before function")
                files[filename] = content
                content = Vector{String}(undef, 0)
                infunction = false
            end
        elseif startswith(line, "function")
            infunction && error("function in function")
            tokens = collect(tokenize(line))
            filename = untokenize(tokens[3])
            infunction = true
            push!(content, line)
        elseif startswith(line, "tmp")
            push!(tmp_nbr_str, line)
        else
            continue
        end
    end
    return files
end

function dynamic_functions(filename)
    files = Dict{String, Vector{String}}()
    tmp_nbr_str = Vector{String}(undef, 0)
    read_dynare_files!(files, tmp_nbr_str, "$(filename)Dynamic.jl")

    global dynamicResid! = make_function(files["dynamicResid!"])
    global dynamicG1! = make_function(files["dynamicG1!"])
    global dynamicG2! = make_function(files["dynamicG2!"])
    global dynamicG3! = make_function(files["dynamicG3!"])
    global dynamicResidTT! = make_function(files["dynamicResidTT!"])
    global dynamicG1TT! = make_function(files["dynamicG1TT!"])
    global dynamicG2TT! = make_function(files["dynamicG2TT!"])
    global dynamicG3TT! = make_function(files["dynamicG3TT!"])

    return eval(Meta.parse(join(tmp_nbr_str, "; ")))
end

function dynamic!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    dynamicResid!(T, residual, y, x, params, steady_state, it_, true)
    return nothing
end

function dynamic!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    dynamicG1!(T, g1, y, x, params, steady_state, it_, true)
    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing
end

function dynamic!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    dynamicG2!(T, g2, y, x, params, steady_state, it_, true)
    dynamicG1!(T, g1, y, x, params, steady_state, it_, false)
    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing 
end

function dynamic!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real}, g3::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    dynamicG3!(T, g3, y, x, params, steady_state, it_, true)
    dynamicG2!(T, g2, y, x, params, steady_state, it_, false)
    dynamicG1!(T, g1, y, x, params, steady_state, it_, false)
    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing
end

function static_functions(filename)
    files = Dict{String, Vector{String}}()
    tmp_nbr_str = Vector{String}(undef, 0)
    read_dynare_files!(files, tmp_nbr_str, "$(filename)Static.jl")

    global staticResid! = make_function(files["staticResid!"])
    global staticG1! = make_function(files["staticG1!"])
    global staticG2! = make_function(files["staticG2!"])
    global staticG3! = make_function(files["staticG3!"])
    global staticResidTT! = make_function(files["staticResidTT!"])
    global staticG1TT! = make_function(files["staticG1TT!"])
    global staticG2TT! = make_function(files["staticG2TT!"])
    global staticG3TT! = make_function(files["staticG3TT!"])

    return  eval(Meta.parse(join(tmp_nbr_str, "; ")))
end


function static!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    staticResid!(T, residual, y, x, params, steady_state, it_, true)
    return nothing
end

function static!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    staticG1!(T, g1, y, x, params, steady_state, it_, true)
    staticResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing
end

function static!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    staticG2!(T, g2, y, x, params, steady_state, it_, true)
    staticG1!(T, g1, y, x, params, steady_state, it_, false)
    staticResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing
end

function static!(T::AbstractVector{<: Real}, residual::AbstractVector{<: Real}, g1::AbstractMatrix{<: Real}, g2::AbstractMatrix{<: Real}, g3::AbstractMatrix{<: Real},
                  y::AbstractVector{<: Real}, x::AbstractMatrix{<: Real}, params::AbstractVector{<: Real}, steady_state::AbstractVector{<: Real}, it_::Int)
    staticG3!(T, g3, y, x, params, steady_state, it_, true)
    staticG2!(T, g2, y, x, params, steady_state, it_, false)
    staticG1!(T, g1, y, x, params, steady_state, it_, false)
    staticResid!(T, residual, y, x, params, steady_state, it_, false)
    return nothing
end

nearbyint(x::Real) = (abs((x)-floor(x)) < abs((x)-ceil(x)) ? floor(x) : ceil(x))

function get_power_deriv(x::Real, p::Real, k::Int64)
    if (abs(x) < 1e-12 && p > 0 && k > p && abs(p-nearbyint(p)) < 1e-12 )
        return 0.0
    else
        dxp = x^(p-k)
        for i = 1:k
            dxp *= p
            p -= 1
        end
        return dxp
    end
end

function load_set_dynamic_auxiliary_variables(modelname::String)
    source = []
    functionstart = false
    for line in readlines("$(modelname)DynamicSetAuxiliarySeries.jl", keep=true)
        if startswith(line, "function")
            functionstart = true
        end
        if functionstart
            push!(source, line)
            if startswith(line, "end")
#                push!(source, "end")
                break
            end
        end
    end
    exp1 = Meta.parse(join(source,"\n"))
#    convert_expression(exp1)
    return @RuntimeGeneratedFunction(exp1)
end
     
function load_dynare_function2(modname::String)::Function
    fun = readlines(modname * ".jl")
    return @RuntimeGeneratedFunction(Meta.parse(join(fun[3:(end-1)], "\n")))
end

function load_steady_state_function(modname::String, compileoption::Bool)
    fun = readlines(modname * ".jl")
    fun[8] = "function steady_state!(ys_::AbstractVector{<: Real}, exo_::AbstractVector{<: Real}, params::AbstractVector{<: Real})"
    return @RuntimeGeneratedFunction(Meta.parse(join(fun[8:length(fun)-1], "\n")))
end

end #end module

