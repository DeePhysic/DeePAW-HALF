# Si primitive-cell CLI bands: EVD versus matrix-free ACC

The single-GPU RTX PRO 6000 test used the two-atom Si primitive cell, a
DeepAW smooth density, its PAW POTCAR, `ENCUT=520 eV`, 12 bands, and 60 points
along `GXWKGLUWLK,UX`.

![Si primitive-cell bands and per-band ACC errors](assets/si_primitive_bandstructure_evd_acc.png)

| solver | wall time | max error, occupied 4 | max error, lowest 8 | max error, all 12 |
|---|---:|---:|---:|---:|
| dense EVD | 7.37 s | reference | reference | reference |
| ACC, at most 40 iterations | 8.90 s | `2.49e-13 eV` | `3.79e-12 eV` | `1.820 eV` |
| ACC, at most 200 iterations | 31.69 s | `2.96e-13 eV` | `2.96e-13 eV` | `0.464 eV` |

Default ACC and EVD give indirect gaps of `0.5635976766544735 eV` and
`0.5635976766545410 eV`, respectively. ACC therefore reproduces the occupied
subspace, low conduction bands, and gap needed by HALF--VASP initialization.
The visible error is confined to the highest requested empty bands. This small
case has only about 1100 plane waves per k point, so its iterative FFT and
orthogonalization overhead makes default ACC 21% slower than dense EVD. ACC is
intended to avoid the `NPL x NPL` dense matrices of larger systems, not to win
this small-matrix timing.

The `half bands` command now accepts `--acc-tol`, `--acc-max-iter`, and
`--acc-block-size`. Raw records are available for
[`EVD`](assets/si_primitive_bands_evd.json),
[`ACC-40`](assets/si_primitive_bands_acc.json), and
[`ACC-200`](assets/si_primitive_bands_acc_strict.json).
