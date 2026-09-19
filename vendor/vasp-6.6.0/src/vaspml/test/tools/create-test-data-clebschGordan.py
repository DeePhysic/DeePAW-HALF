#!/usr/bin/env python3

from sympy.physics.quantum.cg import CG
from sympy import S, N

def cg_data_full(j1max, j2max):
    print("    // clang-format off")
    count = 0
    for j1 in range(0, j1max + 1):
        for m1 in range(-j1, j1 + 1):
            for j2 in range(0, j2max + 1):
                for m2 in range(-j2, j2 + 1):
                    for J in range(0, j1 + j2 + 1):
                        for M in range(-J, J + 1):
                            cg = CG(S(j1), S(m1), S(j2), S(m2), J, M).doit()
                            if cg != 0:
                                #print(f"{count:5d} {j1:3d} {m1:3d} {j2:3d} {m2:3d} {J:3d} {M:3d}, {cg}, {N(cg)}")
                                print(f"    tc->cgs[{count:5d}].coeff = {N(cg):23.16E}; // j1 m1 j2 m2 J M coeff: {j1:3d} {m1:3d} {j2:3d} {m2:3d} {J:3d} {M:3d} {cg}")
                            count += 1
    print("    // clang-format on")

def cg_data_m0(j1max, j2max):
    print("    // clang-format off")
    count = 0
    for j1 in range(0, j1max + 1):
        for j2 in range(0, j2max + 1):
            for J in range(0, j1 + j2 + 1):
                cg = CG(S(j1), 0, S(j2), 0, J, 0).doit()
                if cg != 0:
                    #print(f"{count:5d} {j1:3d} {j2:3d} {J:3d}, {cg}, {N(cg)}")
                    print(f"    tc->cgs[{count:5d}].coeff = {N(cg):23.16E}; // j1 j2 J coeff: {j1:3d} {j2:3d} {J:3d} {cg}")
                count += 1
    print("    // clang-format on")

if __name__ == '__main__':
    #cg_data_full(3, 3)
    cg_data_m0(6, 6)
