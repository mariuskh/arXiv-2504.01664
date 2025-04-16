using QuantumOptics
using Plots
using LaTeXStrings
using JLD2
include("./parameters.jl")

function run_system(state_init::Ket, hamiltonian::Operator; dt::Float64=0.01, t_end::Float64=20.0)
    #Simulate (Schroedinger) time evolution of system,
    #with a given initial state and Hamiltonian.
    #Return time vector and state vector.
    t = [0:dt:t_end;]
    t_out, psi_t = timeevolution.schroedinger(t, state_init, hamiltonian)
    return t_out, psi_t
end

function run_system_dynamic(state_init::Ket, hamiltonian::Function; dt::Float64=0.01, t_end::Float64=20.0)
    #Simulate (Schroedinger) time evolution of system,
    #with a given initial state and Hamiltonian.
    #Return time vector and state vector.
    t = [0:dt:t_end;]
    t_out, psi_t = timeevolution.schroedinger_dynamic(t, state_init, hamiltonian)
    return t_out, psi_t
end

function run_system_master(state_init::Ket, hamiltonian::Operator; rates=nothing, J=nothing, Jd=nothing, dt::Float64=0.01, t_end::Float64=20.0)
    #Simulate the time evolution of system using the Master equation approach,
    #with a given initial state and Hamiltonian.
    #Return time vector and state vector.
     
    if (!isnothing(rates) && isnothing(J))
        # Running with 'rates'
        J, Jd =  get_jump_operators(rates)
    elseif (isnothing(rates) && !isnothing(J))
        # Running with J
        if isnothing(Jd)
            Jd = dagger.(J)
        end
    else
        error("Incorrect input to 'run_system_dynamic_master'. Should either specify 'rates' or 'J' (and optionally 'Jd'), but not both.")
    end
    t = [0:dt:t_end;]

    t_out, rho_t = timeevolution.master(t, state_init, hamiltonian, J, Jdagger=Jd)
    return t_out, rho_t
end

function run_system_dynamic_master(state_init::Ket, hamiltonian::Function; rates=nothing, J=nothing, Jd=nothing, dt::Float64=0.01, t_end::Float64=20.0)
    #Simulate the time evolution of system using the Master equation approach.
    #Should either pass a list of Lindblad operators 'J' (and optionally their adjoints 'Jd'), 
    #or a dictionary with 'rates' (on format Dict{String, Float64}). This is passed to the function get_jump_operators which generates a list of operators.
 
    if (!isnothing(rates) && isnothing(J))
        # Running with 'rates'
        J, Jd =  get_jump_operators(rates)
    elseif (isnothing(rates) && !isnothing(J))
        # Running with J
        if isnothing(Jd)
            Jd = dagger.(J)
        end
    else
        error("Incorrect input to 'run_system_dynamic_master'. Should either specify 'rates' or 'J' (and optionally 'Jd'), but not both.")
    end
    t = [0:dt:t_end;]
    f = function (t, rho)
        return hamiltonian(t, rho), J, Jd
    end
    t_out, rho_t = timeevolution.master_dynamic(t, state_init, f)
    return t_out, rho_t
end

function run_system_dynamic_master(state_init::Operator, hamiltonian::Function; rates=nothing, J=nothing, Jd=nothing, dt::Float64=0.01, t_end::Float64=20.0)
    #Simulate the time evolution of system using the Master equation approach.
    #Should either pass a list of Lindblad operators 'J' (and optionally their adjoints 'Jd'), 
    #or a dictionary with 'rates' (on format Dict{String, Float64}). This is passed to the function get_jump_operators which generates a list of operators.
 
    if (!isnothing(rates) && isnothing(J))
        # Running with 'rates'
        J, Jd =  get_jump_operators(rates)
    elseif (isnothing(rates) && !isnothing(J))
        # Running with J
        if isnothing(Jd)
            Jd = dagger.(J)
        end
    else
        error("Incorrect input to 'run_system_dynamic_master'. Should either specify 'rates' or 'J' (and optionally 'Jd'), but not both.")
    end
    t = [0:dt:t_end;]
    f = function (t, rho)
        return hamiltonian(t, rho), J, Jd
    end
    t_out, rho_t = timeevolution.master_dynamic(t, state_init, f)
    return t_out, rho_t
end

function get_jump_operators(rates::Dict{String, Float64})
    ### Return a vector of jump operators and a vector of their adjoints. 
    ### Based on the rates passed as a parameter. 
    valid_keys = ["boson_loss", "boson_gain", "qubit_decay", "qubit_dephase"]
    #Initialize vectors for jump operators and adjoint.
    J = Vector(undef, length(rates))
    Jd = Vector(undef, length(rates))
    i = 1 #Index for vectors
    for key in eachindex(rates)
        if !(key in valid_keys)
            error("Decay rate '$(key)' not valid. Must be one of: $(valid_keys).")
        end
        J[i] = rates[key]*JUMP_OPERATORS[key]
        Jd[i] = rates[key]*JUMP_ADJOINT_OPERATORS[key]
        i = i+1
    end
    return J, Jd
end

function plot_wigner(t_vec::Vector, psi_t::Vector, t::Float64, grid::Vector, save_fig::Bool=false, fig_title::String="wigner_$(round(t, digits=2)).pdf")
    #Plot Wigner function of a state at a specified time t.
    #Takes as input a time vector and a state evolution vector.
    #Finds the index that corresponds to time t, and calls
    #plot_wigner to plot the state at that index.

    #Check that t_vec and psi_t has same length. 
    if length(t_vec) != length(psi_t)
        error("length(t_vec)=$(length(t_vec)) and length(psi_t)=$(length(psi_t)) must be the same.")
    end
    dt = t_vec[2]-t_vec[1]
    idx = floor(Int64, t/dt)+1 #Add +1 since Julia is 1-indexed.
    plot_title = "Wigner function at time: $(round(t, digits=2))."
    p = plot_wigner(psi_t[idx], grid, plot_title)    
    if save_fig
        savefig(p, fig_title)
    end
    return p 
end 

function plot_wigner(psi::Ket, grid::Vector, plot_title::String="")
    #Plot Wigner function of a given state psi.
    #The colorbar limits (-0.3, 0.3) are chosen based on observation.

    #Density matrix of bosonic system
    rho_b = ptrace(psi, 2)
    #Wigner function
    W = wigner(rho_b, grid, grid)
    return heatmap(grid, grid, W, clims=(-0.3, 0.3), title = plot_title)
end

function plot_wigner(rho::Operator, grid::Vector, plot_title::String="")
    #Plot Wigner function of a given state rho.
    #The colorbar limits (-0.3, 0.3) are chosen based on observation.

    #Density matrix of bosonic system
    rho_b = ptrace(rho, 2)
    #Wigner function
    W = wigner(rho_b, grid, grid)
    return heatmap(grid, grid, W, clims=(-0.3, 0.3), title = plot_title)
end

function animate_wigner(t_vec::Vector, psi_t::Vector, grid::Vector, t_end::Float64=last(t_vec))
    #Animate time evolution of Wigner function. 
    #Stops at t_end (by default the last element of t_vec)

    #Check if illegal value t_end
    if t_end < 0 || t_end > last(t_vec)
        error("Parameter t_end=$(t_end) must be in the range 0 < t_end <= $(last(t_vec)).")
    end

    anim = @animate for t in range(0, stop=t_end, length = 200)
        plot_wigner(t_vec, psi_t, t, grid)
    end
    gif(anim, fps=5)
end 

function state_purity(psi::Ket, subsystem::String="bosonic")
    #Calculate purity tr(rho^2) of subsystem of state psi.
    if subsystem=="bosonic"
        axis=2
    elseif subsystem=="qubit"
        axis=1
    else
        error("Invalid option for subsystem: '$(subsystem)'.")
    end
    rho = ptrace(psi, axis)
    return round(tr(rho*rho), digits = 4)

end

function state_purity(rho::Operator, subsystem::String="bosonic")
    #Calculate purity tr(rho^2) of subsystem of state psi.
    if subsystem=="bosonic"
        axis=2
    elseif subsystem=="qubit"
        axis=1
    else
        error("Invalid option for subsystem: '$(subsystem)'.")
    end
    rho_sub = ptrace(rho, axis)
    return round(tr(rho_sub*rho_sub), digits = 4)

end

function state_purity(psi_t::Vector, subsystem::String="bosonic")
    #Calculate vector of purity of states
    purity = zeros(length(psi_t))
    for i in 1:length(psi_t)
        purity[i] = state_purity(psi_t[i], subsystem)
    end
    return purity
end

function plot_state_purity(t_vec::Vector, psi_t::Vector, subsystem::String="bosonic", save_fig::Bool=false, fig_title::String="$(subsystem)_purity.pdf")
    #Plot state purity vs time.
    #Check that t_vec and psi_t has same length.
    if length(t_vec) != length(psi_t)
        error("length(t_vec)=$(length(t_vec)) and length(psi_t)=$(length(psi_t)) must be the same.")
    end
    purity = state_purity(psi_t, subsystem)
    p = plot(t_vec, purity, ylims = (0, 1.1), title="State purity of $(subsystem) system.", xlabel="t", ylabel=L"\gamma", label=L"\gamma (t)")
    if save_fig
        savefig(p, fig_title)
    end
    return p
end


function project_state(proj::Operator, psi::Ket)
    ### Project state 'psi' on the space spanned by 'proj'. 
    return proj*psi/sqrt(dagger(psi)*proj*psi)
end

function project_state(proj::Operator, rho::Operator)
    ### Project state 'rho' on the space spanned by 'proj'. 
    return proj*rho*proj/tr(proj*rho)
end

function project_state(proj::Operator, psi_t::Vector)
    ### Project a state over time 'psi_t' on the space spanned by 'proj'.
    psi_out = copy(psi_t)
    for i in 1:length(psi_t)
        psi_out[i] = project_state(proj, psi_t[i])
    end
    return psi_out
end

function project_state(proj::Function, t_vec::Vector, psi_t::Vector)
    ### Project a state over time 'psi_t' on the space spanned by 'proj'.
    ### Projector 'proj' is a function of time, so 't_vec' must be included.
    if length(t_vec) != length(psi_t)
        error("length(t_vec)=$(length(t_vec)) and length(psi_t)=$(length(psi_t)) must be the same.")
    end
    psi_out = copy(psi_t)
    for i in 1:length(psi_t)
        psi_out[i] = project_state(proj(t_vec[i]), psi_t[i])
    end
    return psi_out
end

function project_qubit(state::String, psi::Ket)
    ### Project a state on a qubit state. 
    ### Parameter 'state' is a string specifying the qubit state. 
    ### Must be one of the 'valid_states', or an error will be thrown. 
    valid_states = Dict{String, Ket}([
        ("up", spinup(b_spin)),
        ("down", spindown(b_spin)),
        ("x_plus", (spinup(b_spin) + spindown(b_spin))/sqrt(2)),
        ("x_minus", (spinup(b_spin) - spindown(b_spin))/sqrt(2)),
    ])
    if !(state in keys(valid_states))
        error("Parameter 'state' must be one of $(keys(valid_states)).\n 
         Was $(state)")
    end
    # Projector is |state> <state|
    qubit_projector = tensor(valid_states[state], dagger(valid_states[state]))
    #Add identityoperator on bosonic system
    composite_projector = tensor(identityoperator(b_fock), qubit_projector)
    return project_state(composite_projector, psi)
end

function project_qubit(state::String, rho::Operator)
    ### Project a state on a qubit state. 
    ### Parameter 'state' is a string specifying the qubit state. 
    ### Must be one of the 'valid_states', or an error will be thrown. 
    valid_states = Dict{String, Ket}([
        ("up", spinup(b_spin)),
        ("down", spindown(b_spin)),
        ("x_plus", (spinup(b_spin) + spindown(b_spin))/sqrt(2)),
        ("x_minus", (spinup(b_spin) - spindown(b_spin))/sqrt(2)),
    ])
    if !(state in keys(valid_states))
        error("Parameter 'state' must be one of $(keys(valid_states)).\n 
         Was $(state)")
    end
    # Projector is |state> <state|
    qubit_projector = tensor(valid_states[state], dagger(valid_states[state]))
    #Add identityoperator on bosonic system
    composite_projector = tensor(identityoperator(b_fock), qubit_projector)
    return project_state(composite_projector, rho)
end


function project_qubit(state::String, psi_t::Vector)
    ### Project a state over time on a qubit state.
    ### Parameter 'state' is a string specifying the qubit state. 
    ### Must be one of the 'valid_states', or an error will be thrown. 
    valid_states = Dict{String, Ket}([
        ("up", spinup(b_spin)),
        ("down", spindown(b_spin)),
        ("x_plus", (spinup(b_spin) + spindown(b_spin))/sqrt(2)),
        ("x_minus", (spinup(b_spin) - spindown(b_spin))/sqrt(2)),
    ])
    if !(state in keys(valid_states))
        error("Parameter 'state' must be one of $(keys(valid_states)). Was '$(state)'.")
    end
    # Projector is |state> <state|
    qubit_projector = tensor(valid_states[state], dagger(valid_states[state]))
    #Add identityoperator on bosonic system
    composite_projector = tensor(identityoperator(b_fock), qubit_projector)
    return project_state(composite_projector, psi_t)
end

function apply_boson_counter_rotation(t_vec::Vector, psi_t::Vector, w_b::Float64=0.1)
    # Apply boson counter rotation at frequency w_b according to |psi_rot> = exp(i*w_b*n*t)|psi>
    # Useful for looking at wave function relative to rotating frame
    if length(t_vec) != length(psi_t)
        error("length(t_vec)=$(length(t_vec)) and length(psi_t)=$(length(psi_t)) must be the same.")
    end
    psi_out = copy(psi_t)
    for i in 1:length(psi_t)
        psi_out[i] = exp(im*w_b*t_vec[i]*tensor(n, identityoperator(b_spin)))*psi_t[i]
    end
    return psi_out
end


function apply_qubit_counter_rotation(t_vec::Vector, psi_t::Vector, w_q::Float64=1.0, w_d::Float64=0.1, A::Float64 = 2.405*w_d/2.0)
    # Apply qubit counter rotation to see the qubit relative to the rotating frame. 
    # The format of this unitary is very much specific to the driven sunmesh Hamiltonian.
    if length(t_vec) != length(psi_t)
        error("length(t_vec)=$(length(t_vec)) and length(psi_t)=$(length(psi_t)) must be the same.")
    end
    psi_out = copy(psi_t)
    for i in 1:length(psi_t)
        psi_out[i] = exp(im*(A/w_d)*sin(w_d*t_vec[i])*tensor(identityoperator(b_fock), sx)
                    )*exp(im*(w_q/2.0)*t_vec[i]*tensor(identityoperator(b_fock), sz)
                    )*psi_t[i]
    end
    return psi_out
end

function fidelity_over_time(t::Vector, psi_1::Vector, psi_2::Vector, subsystem::String="boson")
    ## Since the vectors can get very large, coarse grain the output
    idx_step = Int(ceil(length(t)/1001))
    t_out = zeros(Int(ceil(length(t)/idx_step)))
    fid = zeros(length(t_out))
    valid_subsystems = Dict{String, Int}([
        ("qubit", 1),
        ("boson", 2),
        ("none", 0),
    ])
    if !(subsystem in keys(valid_subsystems))
        error("Parameter 'subsystem' must be one of $(keys(valid_subsystems)). Was '$(subsystem)'.")
    end
    
    #separate index for new vectors
    j = 1
    for i in 1:idx_step:length(t)
        t_out[j] = t[i]
        if !(subsystem=="none")
            fid[j] = round(real(fidelity(
                                        ptrace(psi_1[i], valid_subsystems[subsystem]
                                     ), ptrace(psi_2[i], valid_subsystems[subsystem]))), digits=5)
        else
            fid[j] = round(real(fidelity(
                                        projector(psi_1[i]
                                     ), projector(psi_2[i]))), digits=5)
        end
        j = j+1
    end
    return t_out, fid 
end

function coarse_grain_state(t_vec::Vector, state_t::Vector, number_of_points::Int=1000)
    if length(t_vec) != length(state_t)
        error("length(t_vec)=$(length(t_vec)) and length(state_t)=$(length(state_t)) must be the same.")
    end
    idx_step = Int(ceil(length(t_vec)/number_of_points))
    
    #Make sure last element is included, even if number of points does not exactly match the length. 
    t_coarse = push!(t_vec[begin:idx_step:end], t_vec[end])
    state_coarse = push!(state_t[begin:idx_step:end], state_t[end])
    return t_coarse, state_coarse
end

function save_state_jld(t::Vector, state::Vector, filename::String; params...)
    jldsave(filename; time=t, state=state, params)
end

function load_state_jld(filename::String)
    local time, state, params
    jldopen(filename) do file
        time = file["time"]
        state = file["state"]
        params = file["params"]
    end
    return time, state, params
end