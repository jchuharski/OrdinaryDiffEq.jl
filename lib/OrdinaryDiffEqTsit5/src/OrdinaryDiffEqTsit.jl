module OrdinaryDiffEqTsit

import OrdinaryDiffEq: alg_order, alg_stability_size, explicit_rk_docstring, trivial_limiter!,
                       OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache, alg_cache,
                       constvalue, uses_uprev, @unpack, perform_step!, calculate_residuals,
                       calculate_residuals!
import Static: False
import MuladdMacro: @muladd
import FastBroadcast: @..
import RecursiveArrayTools: recursivefill!

include("algorithms.jl")
include("alg_utils.jl")
include("tsit_caches.jl")
include("tsit_tableaus.jl")
include("interp_func.jl")
include("tsit_perform_step.jl")

end