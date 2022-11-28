#!/usr/bin/env julia 

using PartialWaveFunctions 

MAX_S2 = 6
cgs = Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }()

for S1=0:2:MAX_S2,
    S2=0:2:MAX_S2,
    S3=0:2:(S1+S2)

    @show S1, S2, S3 

    dS1 = S1 + 1
    dS2 = S2 + 1
    dS3 = S3 + 1

    cgs[S1,S2,S3] = zeros( ComplexF64 , dS1 , dS2 , dS3 )
    cgsview = @view cgs[S1,S2,S3][:,:,:]

    for (i1,s1) in enumerate((-S1):2:S1),
        (i2,s2) in enumerate((-S2):2:S2),
        (i3,s3) in enumerate((-S3):2:S3)

        cgsview[i1,i2,i3] = CG_doublearg(S1,s1,S2,s2,S3,s3) 

    end

    @show cgs[S1,S2,S3]
end





for (i,s) in enumerate( -2:2:2 ) 
    @show i,s 
end
