#!/usr/bin/env python3

import math
from random import randint
import matplotlib.pyplot as plt
import numpy as np
curveorder = 2^2 * 13108968793781547619861935127046491459309155893440570251786403306729687672801

c_h = 505 # constraints to hash 
c_f = 506 # constraints to perform fixed based multiplication
c_v = 2249 # constraints to perform variable based multiplication
c_a = 6 # constraints to perform a simple addition
c_d = 5 # constraints to perform doubling (squaring in pippenger's)
c_decomp = 392 # decompression - see section A.3.3 of zcash specs 
c_check = 20 # check on curve and no small order - see section A.3.3 of zcash specs 
 
def pip(s):
    bases = s
    scalars = [randint(1,curveorder) for _ in range(s)]
    if bases < 32:
        c = 3
    else:
        c = math.ceil(math.log(bases))

    mask = 1 << c
    cur = 0
    num_bits = 254
    windows = 0
    buckets = 0
    adds = 0
    squares = 0
    added = 0

    while cur <= num_bits:
        buckets = ((1 << c) - 1)

        for s in scalars:
            index = s & mask;
            if index != 0:
                # buckets[index - 1].add_assign_mixed(&g);
                adds += 1
                added = 1

        adds += 2 * buckets

        windows += 1

        cur += c;

    if added == 0:
        # try again
        return pip(s)

    squares += c * windows
    adds += 1 * windows

    return adds,squares

def fo(i):
    return '{:,}'.format(int(i)).replace(',', ' ')
    
def poseidon(s):
    return s / 8 * 505

def multiexp(s,decomp=False,check=False):
    adds,doubles = pip(s)
    return c_a * adds + c_d * doubles
    

def merkle(n,f,sectors,challenges):
    depth = math.log(n,8)
    return {
        # *8 because at each lvl we hash 8 elements, not 1
        "total": poseidon(sectors*challenges*depth*8), 
    }

def verkle(n,f,sectors,challenges,log=True):
    N = sectors * challenges
    lf = math.ceil(math.log(f,2))
    d = math.ceil(math.log(n,f)) # number of layers in the vtree excluding leaves
   
    # Verkle tree part
    # verkle tree encoding checks - we only take the x coordinate (we check points are correct)
    # d-1 because we only do this for internal nodes
    verkle = N * (d-1)  * c_check
    
    # Multi point opening part
    # number of tuples (C_,z_i,y_i):
    # N*d since we have N*d (C_,z_i,y_i) tuples since for each challenge we have
    # to open  - the leaves are data, so we dont count them. However, we only
    # count the root once so we can d-1 and add the number of sectors to count
    # the root only once.
    nmexp = (d-1) * N + sectors # Individual verification, no batching
    # need to hash all (C,z_i,y_i) tuples for Fiat Shamir randomness
    # * 3 since we hash each individual element and a point must be hashed via
    # its two coordinates
    multi = poseidon(nmexp*4)
    multi += 3 * (d-1) * N  + sectors# compute g2(t) only one division 
    multi += multiexp(nmexp) # commitment to g1 - points are already checked at previous part
    
    # IPA verification
    ipa = poseidon(lf) # intermediate challenges computation
    ipa += multiexp(lf * 2) # intermediate commitments computation
    ipa += lf * 3 # final polynomial for bases
    ipa += multiexp(f,check=True) # Compute g'
    ipa += f * 3 # Compute b'
    ipa += 2 * c_v # Final check
    
    if log is True:
        print(f"Number of nodes (n): 2^{int(math.log(n,2))}")
        print("Branching factor (f): ",f)
        print("Depth of the tree (d): ",d)
        print(f"Number of openings to msm: {d-1} * {N} + {sectors} = {nmexp}")
        print("Prover storage overhead (all commitments): ", int(sum([f**i for i in range(0,d)]) / 1024 )," kB")
        print()
        
        
        print("TOTAL constraints: ",fo(ipa + multi + verkle))
        print("     Verkle Tree part constraints:   ", fo(verkle))
        print("     Multi points constraints:       ", fo(multi))
        print("     IPA constraints:                ", fo(ipa))
        print("")
        print("--------------------------------------------------")
        print("")
    return {
            "branch":f,
            "total": ipa+multi+verkle,
            "verkle_check": verkle,
            "multipoints":multi,
            "ipa":ipa,
            "msm_size":nmexp,
    }


if __name__ == '__main__':
    n = 2**30
    sectors = 2350
    challenges = 10
    step = 10000
    tests = [ [n,f,sectors,challenges] for f in [2**i for i in range(8,17)] ]
    tests_merkle = [ [n,f,sectors,challenges] for f in [2**i for i in range(8,17)] ]
    results = [verkle(*t) for t in tests]
    results_merkle = [merkle(*t) for t in tests]
    keys = results[0].keys()
    import csv
    import json 
    with open('results_wpost.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
        for row in results_merkle:
            writer.writerow(row)

    # print(json.dumps(results, indent=4, sort_keys=True))

    branch = [d["branch"] for d in results]
    totals = [d["total"] for d in results]
    totals_merkle = [d["total"] for d in results_merkle]
    mpos = [d["multipoints"] for d in results]
    ipas = [d["ipa"] for d in results]
    plt.plot(branch,totals,label="total constraints")
    plt.plot(branch,totals_merkle,label="total constraints MERKLE")
    plt.plot(branch,mpos,label="multi points opening")
    plt.plot(branch,ipas,label="inner product opening")
    plt.xlabel("branching factor")
    plt.ylabel("number of constraints")
    plt.yscale("log")
    plt.legend()
    plt.show()
