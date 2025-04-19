include("../src/utilities.jl")
include("../src/hamiltonians.jl")
include("../src/init_states.jl");


function init_state()
    # Wrapper for initial state \propto (|xi> + |-xi>) \otimes (|e> + |g>)
    xi = 1.0
    psi_plus = (squeeze(b_fock, xi) + squeeze(b_fock, -xi))*coherentstate(b_fock, 0.0)
    
    boson_init = psi_plus/sqrt(dagger(psi_plus)*psi_plus)
    qubit_init = (spinup(b_spin) + spindown(b_spin))/sqrt(2)
    return tensor(boson_init, qubit_init)
end

#System parameters without driving relative to w_m in the rotating frame
w_q = 0.0
w_m = 1.0 
g_1 = 0.0
g_2 = 0.001*w_m
w_d = 0.0
A = 0.0

sys_params= Dict{String, Float64}([
    ("w_q", w_q),
    ("w_b", w_m),
    ("g_1", g_1),
    ("g_2", g_2),
    ("w_d", w_d),
    ("A", A),
])

#Dissipation rates
n_th_m = 1.0
gamma_m = 0.01*g_2 
gamma_q = 1.0*g_2
gamma_phi= 1.0*g_2

decay_rates= Dict{String, Float64}([
    ("qubit_dephase", gamma_phi),
    ("qubit_decay", gamma_q),
    ("boson_loss", (n_th_m + 1)*gamma_m),
    ("boson_gain", n_th_m*gamma_m),
])
decay_params= Dict{String, Float64}([
    ("qubit_dephase", sqrt(decay_rates["qubit_dephase"]/2.0)),
    ("qubit_decay", sqrt(decay_rates["qubit_decay"])),
    ("boson_loss", sqrt(decay_rates["boson_loss"])),
    ("boson_gain", sqrt(decay_rates["boson_gain"])),
])

dt = 0.01
t_end = 1.0/g_2 #Corresponds to 1 qubit lifetime

# Case 1: No measurement 
t, rho = run_system_dynamic_master(init_state(), ham_driven_sunmesh_interaction_frame(sys_params), rates=decay_params, dt=dt, t_end=t_end)
t_c, rho_c = coarse_grain_state(t, rho, 1000);
filename = "state_1_no_drive_no_sz.jld2"
save_state_jld(t_c, rho_c, filename, sys_params=sys_params, decay_rates=decay_rates, dt=dt, t_end=t_end)


# Case 2: Project on sigma_z up 
init_state_up = project_qubit("up", init_state())
t, rho = run_system_dynamic_master(init_state_up, ham_driven_sunmesh_interaction_frame(sys_params), rates=decay_params, dt=dt, t_end=t_end)
t_c, rho_c = coarse_grain_state(t, rho, 1000);
filename = "state_2_no_drive_sz_up.jld2"
save_state_jld(t_c, rho_c, filename, sys_params=sys_params, decay_rates=decay_rates, dt=dt, t_end=t_end)

# Case 3: Project on sigma_z down
init_state_down = project_qubit("down", init_state())
t, rho = run_system_dynamic_master(init_state_down, ham_driven_sunmesh_interaction_frame(sys_params), rates=decay_params, dt=dt, t_end=t_end)
t_c, rho_c = coarse_grain_state(t, rho, 1000);
filename = "state_3_no_drive_sz_down.jld2"
save_state_jld(t_c, rho_c, filename, sys_params=sys_params, decay_rates=decay_rates, dt=dt, t_end=t_end)

# Case 4: Repeated measurements of sigma_z up 
num_meas = 100 #100 evenly spaced measurements.
t = [0.0]
rho = [init_state_up]
t_end_tmp = t_end/num_meas
for i in 1:num_meas
    tmp_init = project_qubit("up", rho[end])
    t_tmp, rho_tmp = run_system_dynamic_master(tmp_init, ham_driven_sunmesh_interaction_frame(sys_params), rates=decay_params, dt=dt, t_end=t_end_tmp)
    t_tmp = t[end] .+ t_tmp
    global t = [t; t_tmp]
    global rho = [rho; rho_tmp]
    
end

t_c, rho_c = coarse_grain_state(t, rho, 1000);
filename = "state_4_no_drive_sz_up_repeated.jld2"
save_state_jld(t_c, rho_c, filename, sys_params=sys_params, decay_rates=decay_rates, dt=dt, t_end=t[end])
