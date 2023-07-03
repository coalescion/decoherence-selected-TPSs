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
    
    return H_ising/2