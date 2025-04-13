# General
from qutip import *
import pennylane as qml
import numpy as np
import numpy.ma as ma
import scipy as sci
import random
import math
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import timeit
import sympy
import itertools

# For Defining Initial Objects
import qutipHam
import total_sys_GGMMs as ggmm
import mod_initstates as init

# For Scrambling Hamiltonians
import matplotlib.pyplot as plt
import mod_timeevol as te
import mod_score as sc




num_env_qbs = 2
n = num_env_qbs+1  # we have one qubit as our system, rest of qubits as environment
dim_sys = 2
dim_env = 2**num_env_qbs
dim_tot = 2**(n)

# coupling constants for Hamiltonian
cc = 1
cx = cy = cz = 1
LU_or_U = 'y' 

initial_hams = [qutipHam.H_ising(n, cc), qutipHam.H_sb(n, cc)]

""" !! currently using spinbath hamiltonian !! (see initial_hams) """
ham_num = 1
initial_ham = initial_hams[ham_num]
print(f'initial ham = {initial_ham}')
tensor_initial_states, initial_states = init.init_states(num_env_qbs)

# # INITIAL STATES OF SYS-ENV (WLOG, SYS in +Z)
# sysRho = np.zeros((2, 2))
# sysRho[0][0] = 1   # specifying pure state where first entry is one
# sysRho = Qobj(sysRho)

# #random dm
# randRho = rand_dm(2**n)
# randRho.dims = [[2]*n]*2 # 2*n input space (n qubits of dim 2), 2*n output space (n qubits of dim 2); 
# randRho = [randRho]*2

# # even superposition environment ready state
# envRho = Qobj(np.ones(shape=(2**(n-1), 2**(n-1))), dims = [[2]*(n-1)]*2)/2**(n-1)
# eSup = tensor(sysRho, envRho)
# eSup = [eSup, eSup]

# Constructing GGMMs (Generalized Gel-Mann Matrices)
totGGMMs = ggmm.constructTotalGGMMs(n)
initial_thetas = [0]*len(totGGMMs)


def scram(nVar, hamiltonian, thetasVar, GGMMs, LUorU):
    '''
   Inputs:
    - hamiltonian: Hamiltonian you want to scramble
    - thetas: Theta parameters corresponding to each GGMMs (LU = 0) or 
      to each selfGGMMs set (LU = 1)
   Outputs:
    - H_scram: Scrambled Hamiltonian 
    '''
    if LUorU == 0:
        U = sc.construct_unitary(nVar, thetasVar, GGMMs)
        return U*hamiltonian*U.dag()
    if LUorU == 1:
        LU = sc.construct_localUnitary(nVar, thetasVar)
        return LU * hamiltonian * LU.dag()
    else:
        return print('Error: entered LU value that is not 1 or 0.')




print(f'SCRAM: {scram(3, initial_ham, initial_thetas, totGGMMs, 1)}')