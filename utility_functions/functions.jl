
using Plots
using DSP
using Statistics
using JSON
using CSV
using DataFrames

include("DetectDecroticNotch.jl")
include("readbin.jl")
include("slide_mean_filter.jl")
include("Variabilities.jl")
include("Sensitivity_and_PPV.jl")
include("Derivate.jl")
include("CalculateStrokeVolume.jl")

