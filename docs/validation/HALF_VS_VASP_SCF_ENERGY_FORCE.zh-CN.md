# HALF 能量和力与完全自洽 VASP 的对比

本文是项目对外采用的能量/力精度对比。参照计算是 `LMAXMIX=-1` 的完全自洽
VASP；固定密度 `ICHARG=11` 结果只保留为实现回归，不进入下面的精度表。

两侧均采用 PBE、`KSPACING=0.35`、`ISPIN=1`、`ISMEAR=0` 和
`SIGMA=0.02 eV`。VASP 使用 `EDIFF=1E-4`、`ALGO=All` 并完成电子自洽。
Si 使用 520 eV，HfO2 使用 521 eV。

## 总能量

| 体系 | 原子数 | HALF Harris (eV/cell) | 完全自洽 VASP (eV/cell) | HALF - VASP (eV/cell) | 每原子差 (meV/atom) |
|---|---:|---:|---:|---:|---:|
| Si | 8 | `-43.3066227505` | `-43.2957423300` | `-0.0108804205` | `-1.360` |
| HfO2 | 12 | `-121.8844407617` | `-122.0024819000` | `+0.1180411383` | `+9.837` |

HfO2 的绝对能量差主要来自与位置无关的 POTCAR 原子参考能 `EATOM` 重构。
它会平移总能量，但不进入 Hamiltonian 本征问题，也不贡献原子力。

## 原子力

| 体系 | 分量 MAE (meV/Angstrom) | 分量 RMSE (meV/Angstrom) | 最大分量差 (meV/Angstrom) | HALF 最大力模 (meV/Angstrom) | VASP 最大力模 (meV/Angstrom) |
|---|---:|---:|---:|---:|---:|
| Si | `0.414` | `0.520` | `0.913` | `0.0158` | `1.096` |
| HfO2 | `10.988` | `14.566` | `28.154` | `95.31` | `122.19` |

这组差异同时包含 HALF 泛函的剩余误差，以及 DeepAW 输入密度与最终 VASP
自洽密度的差异，因此是端到端应用应采用的精度结果。

## 对 VASP 初始波函数和 SCF LOOP 的影响

`EATOM`、Ewald 记账、atomic PAW double counting 和解析力响应都在本征问题
之后计算。优化这些量会改变报告的能量或力，但不会改变 $H$、$S$、本征矢、
传给 VASP 的波函数或 VASP 的 SCF LOOP 数。

只有位于波函数生成路径上的改动才可能减少 LOOP：

1. 让 DeepAW 密度更接近最终自洽密度；
2. 继续提高局域/NLCC 势和 PAW SETDIJ 与 VASP 的一致性；
3. 降低占据子空间残差，提高全带优化收敛质量；
4. 在直接内存传递时严格匹配平面波顺序、k 点、占据和基组截断。

新的 POTCAR 插值进入 $V_{\mathrm{eff}}$，因此可能小幅改善初始子空间。已经完成
的 Harris 力响应以及后续 `EATOM` 修正不在这条路径上，单独应用不会减少 SCF
LOOP。是否真正减少 LOOP 必须用完全相同的 VASP 输入重新做 MP-85 受控 A/B
测试，不能从能量/力改善直接推断。

## 数据位置

- HALF：`/data/limusen/deepaw_half_force_compare/final_vasp_interp_gpu`
- Si 完全自洽：`/data/limusen/deepawchgs/020_mp-149_Si`
- HfO2 完全自洽：
  `/data/limusen/deepawchgs_retry_encut521/054_mp-352_HfO2`
- 可机读汇总：
  [`half_vs_vasp_scf_energy_force.json`](half_vs_vasp_scf_energy_force.json)
