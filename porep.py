#!/usr/bin/env python3
from wpost import verkle,poseidon,fo,merkle
import matplotlib.pyplot as plt
import math
import json
# Porep:
# Challenges: 176 challenges
#     for each challenges: 
#       * prove inclusion proof of challenged node to commD (sha256) ??
#       * prove inclusion proof of the node in the stacked DRG
#           * Reveal labels
#       * prove inclusion proof of the 14 parents
#           * Reveal labels
#       * Prove knowledge of key used to encode replica ??
#       * Prove inclusion proof of nodes in Replica

layers = 11
hash_labels = 30
def porep(n,challenges,parents,branch):
    ## Opening challenged nodes
    ### We open "challenges" nodes in the big column wise verkle tree
    inclusion_drg = challenges 
    ## Opening parent nodes
    inclusion_parents = challenges * parents
    c1 = verkle(n,branch,1,inclusion_drg + inclusion_parents,log=False)
    ## Opening replica
    inclusion_replica = challenges
    c2 = verkle(n,branch,1,inclusion_replica,log=False)
    total = c1['total'] + c2['total']
    print()
    print(f"Branching factor: {branch}")
    print(json.dumps(c1,indent=4))
    print(f"Constraints for inclusion DRG: {fo(c1['total'])}")
    print(f"Constraints for inclusion Replica: {fo(c2['total'])}")
    print(f"Total constraints for inclusion: {fo(total)}")
    print()
    print(f"------------------------------------------")
    return {
        "total": total,
        "drg":c1['total'],
        "replica":c2['total'],
        "branch": branch,
    }

if __name__ == '__main__':
    n = 2**30
    challenges = 176
    parents = 14
    tests_verkle = [ [n,challenges,parents,f] for f in [2**i for i in range(8,17)] ]
    tests_merkle = [ [n,f,challenges,parents] for f in [2**i for i in range(8,17)] ]
    results_verkle = [porep(*t) for t in tests_verkle]
    # challenges * (1 -for challenged nodes- + parents - par. incl. proof.- + 1
    # for replica
    results_merkle = [merkle(t[0],t[1],1,t[2]*(1 + t[3] + 1)) for t in tests_merkle]
    keys = results_verkle[0].keys()
    import csv
    import json
    with open('results_porep.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        for row in results_verkle:
            writer.writerow(row)
        for row in results_merkle:
            writer.writerow(row)

    print(results_merkle)
    branch = [d["branch"] for d in results_verkle]
    inclusion_verkle = [d["total"] for d in results_verkle]
    inclusion_merkle = [d["total"] for d in results_merkle]
    plt.plot(branch,inclusion_verkle,label="inclusion_verkle")
    plt.plot(branch,inclusion_merkle,label="inclusion_merkle")
    plt.xlabel("branching factor")
    plt.ylabel("number of constraints")
    plt.yscale("log")
    plt.legend()
    plt.show()
