using QuantumOptics
include("./parameters.jl")

function state_init_boson_coherent_spin_down(alpha::Float64=2.0)
    #Initial state where bosonic system is in coherent state 
    #and qubit is in spin down.
    boson_init = coherentstate(b_fock, alpha)
    qubit_init = spindown(b_spin)
    return tensor(boson_init, qubit_init)
end

function state_init_boson_coherent_spin_up(alpha::Float64=2.0)
    #Initial state where bosonic system is in coherent state 
    #and qubit is in spin up.
    boson_init = coherentstate(b_fock, alpha)
    qubit_init = spinup(b_spin)
    return tensor(boson_init, qubit_init)
end

function state_init_boson_coherent_spin_plus(alpha::Float64=2.0)
    #Initial state where bosonic system is in coherent state 
    #and qubit is in plus state (up + down)/sqrt(2).
    boson_init = coherentstate(b_fock, alpha)
    qubit_init = (spinup(b_spin) + spindown(b_spin))/sqrt(2)
    return tensor(boson_init, qubit_init)
end

function state_init_boson_coherent_spin_minus(alpha::Float64=2.0)
    #Initial state where bosonic system is in coherent state 
    #and qubit is in minus state.
    boson_init = coherentstate(b_fock, alpha)
    qubit_init = (spinup(b_spin) - spindown(b_spin))/sqrt(2)
    return tensor(boson_init, qubit_init)
end
