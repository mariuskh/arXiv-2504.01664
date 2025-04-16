include("../src/utilities.jl")
include("../src/hamiltonians.jl")
include("../src/init_states.jl")

using SpecialFunctions

#System parameters relative to w_m (in rotating frame)
w_q = 0.0
w_m = 1.0
g_2 = 0.001*w_m 
g_1 = 0.0
w_d = w_m 
A = 2.405*w_d/2.0

sys_params= Dict{String, Float64}([
    ("w_q", w_q),
    ("w_b", w_m),
    ("g_1", g_1),
    ("g_2", g_2),
    ("w_d", w_d),
    ("A", A),
])

t_end=1.2/g_2 # So squeezing strength is approximately 1 at the end.  
dt = 0.01

#Dissipation parameters
n_th_m = 1.0 
gamma_m = 0.01*g_2

#Iterate over different qubit decoherence rates
gamma_q_vec = [0.1*g_2, 1.0*g_2]
gamma_phi_vec = [0.1*g_2, 1.0*g_2]

#Decay rates except qubit, which will be variable
base_decay_rates= Dict{String, Float64}([
    ("boson_loss", (n_th_m + 1)*gamma_m),
    ("boson_gain", n_th_m*gamma_m),
])

base_decay_params= Dict{String, Float64}([
    ("boson_loss", sqrt(base_decay_rates["boson_loss"])),
    ("boson_gain", sqrt(base_decay_rates["boson_gain"])),
])

for i in 1:length(gamma_q_vec)
    for j in 1:length(gamma_phi_vec)
        #Set decay rates for this iteration
        decay_rates = copy(base_decay_rates)
        decay_rates["qubit_decay"] = gamma_q_vec[i]
        decay_rates["qubit_dephase"] = gamma_phi_vec[j]

        decay_params = copy(base_decay_params)
        decay_params["qubit_decay"] = sqrt(decay_rates["qubit_decay"])
        decay_params["qubit_dephase"] = sqrt(decay_rates["qubit_dephase"]/2.0)
        
        #Run system
        t, rho = run_system_dynamic_master(state_init_boson_coherent_spin_minus(0.0), ham_driven_sunmesh_interaction_frame(sys_params), rates=decay_params, dt=dt, t_end=t_end)

        #Save state 
        t_c, rho_c = coarse_grain_state(t, rho, 1000)
        filename = "state_$(i)_$(j)_fidelity_open_system.jld2"
        save_state_jld(t_c, rho_c, filename, sys_params=sys_params, decay_rates=decay_rates, dt=dt, t_end=t_end)
    end
end

#Run a single instance of the RWA (conditional squeezing) Hamiltonian to compare fidelity. 
rwa_params = Dict{String, Float64}([
    ("w_q", 0.0),
    ("w_b", 0.0),
    ("g_1", 0.0),
    ("g_2", besselj(2, 2.405)*g_2),
    ("w_d", 0.0),
    ("A", 0.0),
])

#Conditional squeezing state evolution 
t_cs, psi_cs = run_system(state_init_boson_coherent_spin_minus(0.0), ham_conditional_squeezing(rwa_params), dt=dt, t_end=t_end)
t_cs_c, psi_cs_c = coarse_grain_state(t_cs, psi_cs, 1000)
filename = "state_cs_fidelity_open_system.jld2"
save_state_jld(t_cs_c, psi_cs_c, filename, sys_params=sys_params, dt=dt, t_end=t_end)