include("../src/utilities.jl")
include("../src/hamiltonians.jl")
include("../src/init_states.jl")

using SpecialFunctions

#System parameters relative to w_m in the rotating frame
w_q = 0.0
w_m = 1.0 
g_1 = 0.0
w_d = w_m 
A = 2.405*w_d/2.0

dt = 0.01

#Params (except g_2 which will be variable) in the rotating frame.
base_params= Dict{String, Float64}([
    ("w_q", 0.0),
    ("w_b", w_m),
    ("g_1", g_1),
    ("w_d", w_d),
    ("A", A),
])

#Params (except g_2) with the rotating wave approximation.
base_params_RWA = Dict{String, Float64}([
    ("w_q", 0.0),
    ("w_b", 0.0),
    ("g_1", 0.0),
    ("w_d", 0.0),
    ("A", 0.0),
])

g_2_vec = [0.1*w_m, 0.01*w_m, 0.001*w_m]

for i in 1:length(g_2_vec)
    #Set params for this iteration
    params = copy(base_params)
    params["g_2"] = g_2_vec[i]
    params_RWA = copy(base_params_RWA)
    params_RWA["g_2"] = g_2_vec[i]*besselj(2, 2.0*A/w_d)

    #Adjust t_end to fit with coupling 
    t_end = 1.2/params["g_2"]

    # Full Hamiltonian
    t, psi = run_system_dynamic(state_init_boson_coherent_spin_minus(0.0), ham_driven_sunmesh_interaction_frame(params), dt=dt, t_end=t_end)   
    # Save state 
    t_c, psi_c = coarse_grain_state(t, psi, 1000)
    filename = "state_$(i)_fidelity_rwa.jld2"
    save_state_jld(t_c, psi_c, filename, sys_params=params, dt=dt, t_end=t_end)
    
    #RWA 
    t_rwa, psi_rwa = run_system(state_init_boson_coherent_spin_minus(0.0), ham_conditional_squeezing(params_RWA), dt=dt, t_end=t_end)
    #Save state
    t_rwa_c, psi_rwa_c = coarse_grain_state(t_rwa, psi_rwa, 1000)
    filename = "state_$(i)_cs_fidelity_rwa.jld2"
    save_state_jld(t_rwa_c, psi_rwa_c, filename, sys_params=params_RWA, dt=dt, t_end=t_end)
end
