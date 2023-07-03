from qutip import *
import numpy as np
import mod_GGMM as gg
import itertools
from itertools import permutations

def constructSelfGGMMs(n, j):
    '''
        inputs:
        - n: number of total qubits
        - j: qubit; 0 <-> qubit 1, 1 <-> qubit 2, etc. Implictly, slot 0 <-> qubit 1, etc.
    
        outputs:
        - selfGGMMs
    '''
    A = gg.construct_ggmm_sub(2)
    selfGGMMs = [0]*3
    for i in range(len(A)):
        tensor_list = [qeye(2)]*n
        tensor_list[j] = A[i]
        selfGGMM = tensor(tensor_list)
        selfGGMMs[i] = selfGGMM
    
    return selfGGMMs

def constructAllSelfGGMMs(n):
    allSelfGGMMs = []
    for j in range(n):
        S = constructSelfGGMMs(n, j)
        allSelfGGMMs += S

    return allSelfGGMMs
    
def constructTotalGGMMs(n):
    '''
        inputs:
        - n: number of total qubits
    
        outputs:
        - totGGMMs
    '''

    A = gg.construct_ggmm_sub(2)
    A.insert(0, qeye(2))
    allGGMMs = list(itertools.product(A, repeat=n)) 
    for i in range(len(allGGMMs)):
        allGGMMs[i] = list(allGGMMs[i])
        allGGMMs[i] = tensor(allGGMMs[i])
    allGGMMs.remove(tensor([qeye(2)]*n))
    
    return allGGMMs

def constructTotalMinusSelfGGMMs(n):
    B = constructTotalGGMMs(n)
    C = ggmm.constructAllSelfGGMMs(n) + [tensor([qeye(2)]*n)]
    D = [e for e in B if e not in C]

    return totalMinusAllSelfGGMMs
    
    