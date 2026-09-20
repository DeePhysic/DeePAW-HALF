# DeePAW-HALF acceleration of VASP SCF convergence

## Results

For a VASP 6.6.0 SCF calculation of HfO₂, the benchmark used
`KSPACING=0.35`, 36 irreducible k points, `ALGO=All`, and `EDIFF=1E-4`.
DeePAW-HALF initialization was compared with native VASP SAD initialization.

### `LMAXPAW=-1`

| Metric | DeePAW-HALF | VASP SAD | HALF acceleration |
|---|---:|---:|---:|
| Electronic iterations | **4** | 16 | **75% fewer; 4.0× iteration-convergence speedup** |
| End-to-end wall time | **193.44 s** | 266.96 s | **1.38× speedup; 27.54% less time** |
| Final energy | -121.10188267 eV | -121.10187886 eV | Difference: 3.813×10⁻⁶ eV/cell |

### Default VASP LMAX

| Metric | DeePAW-HALF | VASP SAD | HALF acceleration |
|---|---:|---:|---:|
| Electronic iterations | **5** | 16 | **68.75% fewer; 3.2× iteration-convergence speedup** |
| End-to-end wall time | **145.76 s** | 262.70 s | **1.80× speedup; 44.51% less time** |
| Final energy | -121.04035930 eV | -121.04037650 eV | Difference: 1.7194×10⁻⁵ eV/cell |

## Conclusion

DeePAW-HALF reduced the VASP electronic iteration count from 16 to 4–5, a
reduction of 68.75%–75%. Including HALF initial-wavefunction generation, the
measured end-to-end speedup was **1.38×–1.80×**, reducing wall time by
**27.54%–44.51%**.

For both settings, the final HALF and SAD energies differed by less than
`2×10⁻⁵ eV/cell`. DeePAW-HALF therefore preserved the converged VASP result
while substantially accelerating convergence.

Each path currently has one timing sample. The wall times report the observed
end-to-end effect, while the electronic iteration count is the more stable
acceleration metric.

中文版：
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md`](HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md).
