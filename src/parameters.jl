using QuantumOptics

# Bases
N_cutoff = 50 
b_fock = FockBasis(N_cutoff)
b_spin = SpinBasis(1//2)

#Operators
a = destroy(b_fock)
at = create(b_fock)
n = number(b_fock)

sm = sigmam(b_spin)
sp = sigmap(b_spin)
sz = sigmaz(b_spin)
sx = sigmax(b_spin)
sy = sigmay(b_spin)


JUMP_OPERATORS = Dict{String, Operator}([
    ("boson_loss", tensor(a, identityoperator(b_spin))),
    ("boson_gain", tensor(at, identityoperator(b_spin))),
    ("qubit_decay", tensor(identityoperator(b_fock), sm)),
    ("qubit_dephase", tensor(identityoperator(b_fock), sz)),
])

JUMP_ADJOINT_OPERATORS = Dict{String, Operator}([
    ("boson_loss", tensor(at, identityoperator(b_spin))),
    ("boson_gain", tensor(a, identityoperator(b_spin))),
    ("qubit_decay", tensor(identityoperator(b_fock), sp)),
    ("qubit_dephase", tensor(identityoperator(b_fock), sz)),
])