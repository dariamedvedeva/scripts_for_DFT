import numpy as np
import math
from scipy.integrate import simps
import sys
import os

print("\n")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("\n")
print(" Script for construction of GF from DFT DOS (Sigma is 0)")
print("\n")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("\n")

path_to_file = "/User/path/"
filename     = path_to_file + "cub.dos"
beta         = 100
idelta       = 1j * math.pi / beta
num_of_orbs  = 10

def main ():
    f = open(filename, "r")
    lines = f.readlines()

    freq    = []
    ro_t2g  = []
    ro_eg   = []
    Delta_t2g  = []
    Delta_eg   = []
    

    for line in lines:
        # ro is density of states from DFT calculation
        freq.append(float(line.split()[0]))
        ro_t2g.append(float(line.split()[1]))
        ro_eg.append(float(line.split()[-1]))
    f.close()

    array_length_freq   = len(freq)
    array_length_dos    = len(ro_t2g)
    
    if (array_length_freq == array_length_dos):
        a = 1
    else:
        os.exit()
    
    # Green's functions
    GF_t2g  = construct_GF_from_DOS(freq, ro_t2g, "t2g")
    GF_eg   = construct_GF_from_DOS(freq, ro_eg,  "eg")

    for i in range(array_length_freq):
        Delta_t2g.append(-(1.0/GF_t2g[i]).imag)
        Delta_eg.append(-(1.0/GF_eg[i]).imag)

    save_Delta_to_file(freq, Delta_t2g, "Delta_t2g")
    save_Delta_to_file(freq, Delta_eg, "Delta_eg")

def save_GF_to_file(freq, func, pointer):
    # function is complex
    f = open("GF_" + pointer + ".dat", "w")
    length = len(freq)
    for i in range(length):
        f.write(str(freq[i]))
        f.write('\t')
        f.write(str(func[i].real))
        f.write('\t')
        f.write(str(func[i].imag))
        if (i != length):
            f.write('\n')
    f.close()

def save_Delta_to_file(freq, func, pointer):
    # function is complex
    f = open(pointer + ".dat", "w")
    length = len(freq)
    for i in range(length):
        f.write(str(freq[i]))
        f.write('\t')
        f.write(str(func[i]))
        if (i != length):
            f.write('\n')
    f.close()

def construct_GF_from_DOS(freq, dos, pointer):
    array_length_freq   = len(freq)
    
    GF = []
    total_occup = simps(dos, freq)
    
#    freqs =
    function_for_this_freq = []
    
    print("total occupation of " + pointer + " is " + str(total_occup))
    for i in range(array_length_freq):
        for j in range(array_length_freq):
            function_for_this_freq.append(dos[j]/(freq[i] + idelta - freq[j]))
        value = simps(function_for_this_freq, freq)
        GF.append(value)
        function_for_this_freq = []
    save_GF_to_file(freq, GF, pointer)
    print("GF for " + pointer + " was written into file.\n")
    return GF

main()
