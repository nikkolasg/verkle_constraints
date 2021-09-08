#!/usr/bin/env python3

import math

def bench(n,f,N):
    lf = math.ceil(math.log(f,2))
    d = math.ceil(math.log(n,f)) # number of layers in the vtree
    c_h = 505 # constraints to hash 
    c_f = 506 # constraints to perform fixed based multiplication
    c_v = 2249 # constraints to perform variable based multiplication
    c_a = 6 # constraints to perform a simple addition
    c_d = 5 # constraints to perform doubling (squaring in pippenger's)
    c_decomp = 392 # decompression - see section A.3.3 of zcash specs 
    c_check = 20 # check on curve and no small order - see section A.3.3 of zcash specs 
    
    def fo(i):
        return '{:,}'.format(int(i)).replace(',', ' ')
    
    def multiexp(s,decomp=False,check=False):
        ratio = 52 # exp. ratio found, between the size and the number of additions
        order = 256
        cost = s * ratio * c_a # number of additions
        cost += 2 * order  * c_d # number of point doubling roughly
        if check == True:
            cost = s * c_check # check points given are correct
        if decomp == True:
            cost += c_decomp * s # (de)compress the points
        return cost
    
    # Verkle tree part
    # verkle tree encoding checks - we only take the x coordinate (we check points are correct)
    verkle = N * (d)  * c_check
    
    # Multi point opening part
    multi = 3 * c_h * N / 8 # hashes
    multi += 3 * d * N # compute g2(t)
    nmexp = d* N # Individual verification, no batching
    # knmexp = 
    multi += multiexp(nmexp) # commitment to g1 - points are already checked at previous part
    
    # IPA verification
    ipa = lf * c_h # intermediate challenges computation
    ipa += multiexp(lf * 2) # intermediate commitments computation
    ipa += lf * 3 # final polynomial for bases
    ipa += multiexp(f,check=True) # Compute g'
    ipa += f * 3 # Compute b'
    ipa += 2 * c_v # Final check
    
    print(f"Number of nodes (n): 2^{int(math.log(n,2))}")
    print("Branching factor (f): ",f)
    print("Depth of the tree (d): ",d)
    print(f"Number of openings to compress: d*N = {d}*{N} = {d*N}")
    print("Prover storage overhead (all commitments): ", int(sum([f**i for i in range(0,d)]) / 1024 )," kB")
    print()
    
    
    print("TOTAL constraints: ",fo(ipa + multi + verkle))
    print("     Verkle Tree part constraints:   ", fo(verkle))
    print("     Multi points constraints:       ", fo(multi))
    print("     IPA constraints:                ", fo(ipa))
    print("")
    print("--------------------------------------------------")
    print("")

tests = [
    [
        2**30, # number of nodes in total
        34000, # width of each node in the vtree
        2350*10, # total number of verkle tree openings - 10 challenge per sector, 2350 sectors per partitions
    ],
    [
        2**30, # number of nodes in total
        4096, # width of each node in the vtree
        2350*10, # total number of verkle tree openings - 10 challenge per sector, 2350 sectors per partitions
    ],
    [
        2**31, # number of nodes in total
        4096, # width of each node in the vtree
        2350*10, # total number of verkle tree openings - 10 challenge per sector, 2350 sectors per partitions
    ]

]

[bench(t[0],t[1],t[2]) for t in tests]