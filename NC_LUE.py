from qutip import *
import qutipHam
import numpy as np
import scipy as sci
import total_sys_GGMMs as ggmm
import mod_initstates as init
import matplotlib.pyplot as plt
import mod_initstates as init
import mod_timeevol as te
import mod_score as sc
import math


def compressHam(nVar, hamiltonian):
    '''
    Compresses spectrum of hamiltonian into interval [0, 1]. Additionally, 
    makes trace of compressed hamiltonian 1 if possible. 
    
    Output:
    - processedHam: Hermitian, positive-definite, and Tr = 1 hamiltonian with
    eigenvalues in the interval [0, 1], i.e. an imposter density matrix.
    '''
    
    eigvals = hamiltonian.eigenstates()[0]
    eigMax = max(eigvals)
    eigMin = min(eigvals)
    eigRange = eigMax - eigMin
    a = 0
    b = 1
    d = 2**nVar
    alpha = (b-a)/(eigRange)
    alphabeta = a - alpha*eigMin
    tr = (d * -(eigMin/(eigMax-eigMin)))
    
    processedHam = (alpha*hamiltonian + alphabeta*tensor([qeye(2)]*nVar))/tr
    
    #if processedHam.tr() == 1:
    return processedHam
    
    #else:
        #return print(f'Error: Tr(processedHam) = {processedHam.tr()}')

def NC(nVar, hamiltonian1, hamiltonian2):
    '''
    Necessary check for LUE from Martins. 0 == not LUE, 2 := uncertain.
    '''
    diagList1 = []
    diagList2 = []
    nullMatrix = np.zeros((2, 2))
    
    cHam1 = compressHam(nVar, hamiltonian1)
    cHam2 = compressHam(nVar, hamiltonian2)
    
    for i in range(nVar):
        diag1 = np.diag((cHam1.ptrace(i)).eigenenergies())
        diag2 = np.diag((cHam2.ptrace(i)).eigenenergies())
        diff = diag2 - diag1
        if (np.absolute(diff) < 10**(-14)).all():
            continue
        else:
            return 0
    
    return 2

# ################################
# # For testing
# num_env_qbs = 1
# n = num_env_qbs+1
# cc = 1

# dim_tot = 2**n
# el = dim_tot**2-1



# initial_ham = qutipHam.H_ising(n, cc)
# totGGMMs = ggmm.constructSelfGGMMs(n, 0)

# def randomThetas(LUorU):
#     lub = math.pi/4
#     if LUorU == 0:
#         # global unitary (see scram fn)
#         scram_thetas = np.random.uniform(-lub, lub, el)
#         return scram_thetas
#     if LUorU == 1:
#         # local unitary
#         scram_thetas = np.random.uniform(-lub, lub, 3*n)
#         return scram_thetas
#     else:
#         return print('Error: entered LU value that is not 1 or 0.')

# def scram(nVar, hamiltonian, thetasVar, GGMMs, LUorU):
#     '''Inputs:
#     - hamiltonian: Hamiltonian you want to scramble
#     - thetas: Theta parameters corresponding to each GGMMs (unitary, LU = 0) or 
#         to each selfGGMMs set (local unitary, LU = 1)
#     Outputs:
#     - H_scram: Scrambled Hamiltonian '''

#     if LUorU == 0:
#         #U = haar_SU(nVar)
#         U = sc.construct_unitary(nVar, thetasVar, GGMMs)
#         return U*hamiltonian*U.dag()
#     if LUorU == 1:
#         LU = sc.construct_localUnitary(nVar, thetasVar)
#         return LU * hamiltonian * LU.dag()
#     else:
#         return print('Error: entered LU value that is not 1 or 0.')

# # Actual test
# l = 0 
# while l != 100:
#    randHam1 = scram(n, initial_ham, randomThetas(0), totGGMMs, 0)
#    randHam2 = scram(n, initial_ham, randomThetas(0), totGGMMs, 0)
#    print('SUCCESS TO HERE')
#    print(NC(n, randHam1, randHam2))