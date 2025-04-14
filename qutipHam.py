from qutip import *

def H_sb(n, c):
    '''
        Defining spin bath model Hamiltonian.
        Inputs:
            - n: total number of qbs making up universe (system and environment)
            - c: determines strength of interaction terms (can be generalized to a vector)
    '''
    H_sb = 0
    for i in range(0, n):
        term = [qeye(2)]*n
        term[0] = 1/2 * sigmaz()
        if i == 0:
            continue
        term[i] = c*sigmaz()
        tensor_term = tensor(term)
        H_sb += tensor_term
    
    return H_sb


def H_ising(n, c):
    '''
        Defining 1D ising model Hamiltonian w/ periodic boundary conditions. 
        The coupling strength is given by c.
        The Hamiltonian is defined as: H = (c/2) * sum_over_i(simgaz_i)*(sigmaz_(i+1))
        No self-Hamiltonian, only interaction Hamiltonian!
    '''
    H_ising = 0
    for i in range(0, n):
        term = [qeye(2)]*n
        if i == n-1:
            term[i] = sigmaz()
            term[0] = sigmaz()
            tensor_term = tensor(term)
            H_ising += tensor_term 
        else:
            term[i] = sigmaz()
            term[i+1] = sigmaz()
            tensor_term = tensor(term)
            H_ising += tensor_term 
    
    return (c/2) * H_ising

def H_cnot(n):
    '''
        CNOT-style interaction Hamiltonian.
        Control is always qubit 0; targets are the rest.
    '''
    H_cnot = 0
    proj1 = basis(2, 1) * basis(2, 1).dag()  # |1><1|
    for i in range(1, n):
        ops = [qeye(2)] * n
        ops[0] = proj1           # control qubit
        ops[i] = sigmax()        # target qubit
        H_cnot += tensor(ops)
    return H_cnot


    