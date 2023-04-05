#!/usr/bin/env julia 

#%% INITIAL SETUP
include( "../modules/antisymmetry_numeric.jl" )

orbital = "Eg"
N = 2

cg_path = "/home/acalvo/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
asym_path = "/home/acalvo/Bulegoa/AntiSymmetricPart/Oh/$(orbital)_julia/"

compute_asymstates_allN( orbital , cg_path , asym_path ; verbose=true )

##%%
#ISisr = compute_asymstates_N( orbital , N , cg_path , asym_path)
#


