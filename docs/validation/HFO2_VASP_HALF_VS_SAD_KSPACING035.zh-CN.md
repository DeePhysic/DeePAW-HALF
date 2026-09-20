# DeePAW-HALF 加速 VASP 自洽收敛

## 测试结果

在 HfO₂ 的 VASP 6.6.0 自洽计算中，测试采用 `KSPACING=0.35`、36 个不可约
k 点、`ALGO=All` 和 `EDIFF=1E-4`，比较 DeePAW-HALF 初始化与 VASP 原生
SAD 初始化的收敛速度。

### `LMAXPAW=-1`

| 指标 | DeePAW-HALF | VASP SAD | HALF 加速效果 |
|---|---:|---:|---:|
| 电子迭代数 | **4** | 16 | **减少 75%，迭代收敛加速 4.0×** |
| 端到端 wall time | **193.44 s** | 266.96 s | **1.38× 加速，节省 27.54%** |
| 最终能量 | -121.10188267 eV | -121.10187886 eV | 相差 3.813×10⁻⁶ eV/cell |

### VASP 默认 LMAX

| 指标 | DeePAW-HALF | VASP SAD | HALF 加速效果 |
|---|---:|---:|---:|
| 电子迭代数 | **5** | 16 | **减少 68.75%，迭代收敛加速 3.2×** |
| 端到端 wall time | **145.76 s** | 262.70 s | **1.80× 加速，节省 44.51%** |
| 最终能量 | -121.04035930 eV | -121.04037650 eV | 相差 1.7194×10⁻⁵ eV/cell |

## HALF 独立计算的 HfO₂ 能带

下图由 DeePAW-HALF 独立完成，未启动或进入 VASP 流程。计算采用 PBE、
500 eV 截断能、MIMIC_US、60 条能带和 100 个高对称路径采样点；能量以
价带顶（VBM）为零点。

![DeePAW-HALF 独立计算的 HfO2 能带](assets/hfo2_half_direct_bs.png)

HALF 得到的采样带隙为 **4.576 eV**，Γ 点直接带隙为 **4.629 eV**。

### 与 VASP 的能带对比

VASP 6.6.0 与独立 HALF 计算使用相同的 100 个 k 点和 60 条能带；两组能量
分别以各自的 VBM 对齐。

![HfO2 的 HALF 与 VASP 能带对比](assets/hfo2_half_vasp_band_comparison.png)

| 指标 | DeePAW-HALF | VASP 6.6.0 | 差值 |
|---|---:|---:|---:|
| 采样带隙 | 4.576341 eV | 4.575784 eV | 0.557 meV |
| Γ 点直接带隙 | 4.629091 eV | 4.631752 eV | 2.661 meV |

HALF 在未进入 VASP 流程的情况下复现了 VASP 的 HfO₂ 能带色散和带隙。

## 结论

DeePAW-HALF 将 VASP 的电子迭代数从 16 次降低到 4–5 次，即减少
68.75%–75%。包含 HALF 初始波函数生成在内，实测端到端加速为
**1.38×–1.80×**，wall time 减少 **27.54%–44.51%**。

两种设置下，HALF 与 SAD 的最终能量差均小于 `2×10⁻⁵ eV/cell`，说明 HALF
在显著加快收敛的同时保持了最终 VASP 结果的一致性。

每条路径目前只有一次计时样本，因此端到端时间用于展示本次实测效果；电子迭代
数是更稳定的加速指标。

English version:
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.md`](HFO2_VASP_HALF_VS_SAD_KSPACING035.md)。
