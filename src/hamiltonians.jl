using QuantumOptics
include("./parameters.jl")

function ham_sunmesh(w_b::Float64, w_q::Float64, g_1::Float64, g_2::Float64)
    #SUNMESH Hamiltonian. With first (second) order term proportional to g_1 (g_2).
    H_boson = w_b*n
    H_qubit = w_q*sz/2
    H_int_1 = g_1*(tensor(a + at, sp + sm))
    H_int_2 = g_2*(tensor((a + at)*(a + at), sz))
    return tensor(H_boson, identityoperator(b_spin)) + tensor(identityoperator(b_fock), H_qubit) + H_int_1 + H_int_2
end

function ham_sunmesh(params::Dict{String, Float64})
    #SUNMESH Hamiltonian. With params passed as a dict for better transparency. 
    required_params = ["w_b", "w_q", "g_1", "g_2"]
    for i in 1:length(required_params)
        if !(required_params[i] in keys(params))
            error("Parameter '$(param)' missing from params.")
        end
    end
    return ham_sunmesh(params["w_b"], params["w_q"], params["g_1"], params["g_2"])
end

function ham_conditional_squeezing(w_b::Float64=1.0, w_q::Float64=2.0, g::Float64=1.0)
    #Hamiltonian with conditional squeezing. Coupling between sigma_z and squeezing.
    H_boson = w_b*n
    H_qubit = w_q*sz/2
    H_int = g*(tensor(a*a +at*at, sz))
    return tensor(H_boson, identityoperator(b_spin)) + tensor(identityoperator(b_fock), H_qubit) + H_int
end

function ham_conditional_squeezing(params::Dict{String, Float64})
    #Hamiltonian with conditional squeezing. With params passed as a dict for better transparency. 
    required_params = ["w_b", "w_q", "g_2"]
    for i in 1:length(required_params)
        if !(required_params[i] in keys(params))
            error("Parameter '$(param)' missing from params.")
        end
    end
    return ham_conditional_squeezing(params["w_b"], params["w_q"], params["g_2"])
end

function ham_driven_sunmesh(w_b::Float64=0.1, w_q::Float64=1.0, g_1::Float64=0.0, g_2::Float64=0.001, w_d::Float64=w_b, A::Float64=2.405*w_d/2.0)
   #Driven Hamiltonian. Must return a Function of (t, psi) instead of an operator.

   #Time independent terms can be taken outside
    H_sunmesh = ham_sunmesh(w_b, w_q, g_1, g_2)
    
    H_t = function (t, psi)
        #Time dependent term
        H_drive = A*cos(w_d*t)*(cos(w_q*t)*sx + sin(w_q*t)*sy)
        return H_sunmesh + tensor(identityoperator(b_fock), H_drive) 
    end
    return H_t
end

function ham_driven_sunmesh(params::Dict{String, Float64})
    #Driven Hamiltonian. With params passed as a dict for better transparency. 
    required_params = ["w_b", "w_q", "g_1", "g_2", "w_d", "A"]
    for i in 1:length(required_params)
        if !(required_params[i] in keys(params))
            error("Parameter '$(required_params[i])' missing from params.")
        end
    end
    return ham_driven_sunmesh(params["w_b"], params["w_q"], params["g_1"], params["g_2"], params["w_d"], params["A"])
 end

 function ham_driven_sunmesh_interaction_frame(w_b::Float64=0.1, w_q::Float64=1.0, g_1::Float64=0.0, g_2::Float64=0.001, w_d::Float64=w_b, A::Float64=2.405*w_d/2.0)
    #Driven SUNMESH Hamiltonian in the frame where w_q = 0. Driving term is different in this frame.
    #Time independent terms can be taken outside
     
    H_t = function (t, psi)
        #Time dependent coupling term 
        H_c = g_2*tensor((at*at*exp(im*2*w_b*t) + a*a*exp(-im*2*w_b*t) + 2*at*a + one(b_fock)), sz)
        #Time dependent driving term
        H_drive = tensor(identityoperator(b_fock), A*cos(w_d*t)*sx)
        return H_c + H_drive
    end
    return H_t
 end

 function ham_driven_sunmesh_interaction_frame(params::Dict{String, Float64})
    #Driven SUNMESH Hamiltonian in the frame where w_q = 0. Driving term is different in this frame.
    #Time independent terms can be taken outside
    required_params = ["w_b", "w_q", "g_1", "g_2", "w_d", "A"]
    for i in 1:length(required_params)
        if !(required_params[i] in keys(params))
            error("Parameter '$(required_params[i])' missing from params.")
        end
    end
    return ham_driven_sunmesh_interaction_frame(params["w_b"], params["w_q"], params["g_1"], params["g_2"], params["w_d"], params["A"])
 end