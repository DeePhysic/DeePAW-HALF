# Harris 力的 L 通道与密度响应优化

本文记录 `force-accuracy-experiments` 分支对 PAW augmentation 角动量通道和
Harris 固定密度力的进一步实现、诊断与验证。测试均在 Pro 6000 裸机上完成，
使用 CUDA 13.2、NVHPC 26.5、PBE、VASP `IBZKPT`；Si 为 520 eV/24 bands，
HfO2 为 521 eV/80 bands。

## 1. L 通道必须一致截断

augmentation 核展开为

$$
Q_{ij}^{I}(\mathbf r)=\sum_{LM}C_{ij}^{LM}Q_{ij}^{L}
g_L(r_I)Y_{LM}(\widehat{\mathbf r_I}),
\qquad \mathbf r_I=\mathbf r-\mathbf R_I.
$$

同一个展开同时进入三个量：

$$
D_{ij}^{I}=D_{ij}^{\mathrm{ion}}+
\int V_{\mathrm{eff}}(\mathbf r)Q_{ij}^{I}(\mathbf r)\,d\mathbf r,
$$

$$
\rho_{\mathrm{aug}}(\mathbf r)=
\sum_{Iij}P_{ij}^{I}Q_{ij}^{I}(\mathbf r),
$$

$$
\mathbf F_I^{\mathrm{aug}}=-\sum_{ij}P_{ij}^{I}
\frac{\partial D_{ij}^{I}}{\partial\mathbf R_I}.
$$

因此实验性的 `HALF_QDEP_LMAX` 同时作用于 SETDIJ、augmentation density 和
解析 `dD/dR`；只截断其中一项会破坏能量—力一致性。默认值 `-1` 使用 POTCAR
允许的全部通道：

```bash
cmake -DHALF_QDEP_LMAX=-1 ...  # 正式默认：全部 L
cmake -DHALF_QDEP_LMAX=2  ...  # 仅用于通道诊断
```

HfO2 的 `Lmax=0,2,4` 扫描得到逐位相同的能量与力。这表明该高对称结构在完整
k 点加权后的非球形 onsite multipole 消失；原有 HfO2 力误差并非高 L 求积精度
导致，继续增加角向积分点或保留更高 L 都不会改善它。

## 2. 真正缺失的 Harris 响应

令输入固定密度为 $n_{\mathrm{in}}$，由固定势本征态重建的输出密度为
$n_{\mathrm{out}}$，并定义

$$
\Delta n=n_{\mathrm{out}}-n_{\mathrm{in}}.
$$

Harris 泛函对密度的线性响应场为

$$
g_{H\!xc}(\mathbf r)=
v_H[\Delta n](\mathbf r)+
f_{xc}[n_{\mathrm{in}}+n_{\mathrm{core}}]\,\Delta n(\mathbf r).
$$

相应的离子力修正是

$$
\mathbf F_I^{\mathrm{Harris}}=-\int g_{H\!xc}(\mathbf r)
\frac{\partial n_{\mathrm{atom}}^I(\mathbf r-\mathbf R_I)}
{\partial\mathbf R_I}\,d\mathbf r,
$$

其中 $n_{\mathrm{atom}}$ 来自 POTCAR 的 `PSPRHO`。`PSPCOR` 只属于独立的
NLCC/`FORCOR` 项，不能在这里代替 `PSPRHO`。

旧实现只计算了 XC kernel，并错误地与 `PSPCOR` 收缩，遗漏了通常占主导的
Hartree 响应。新实现对完整有效势做密度方向的中心导数：

$$
g_{H\!xc}\approx
\frac{V_{\mathrm{eff}}[n_{\mathrm{in}}+t\Delta n]
-V_{\mathrm{eff}}[n_{\mathrm{in}}-t\Delta n]}{2t},
\qquad t=10^{-4},
$$

离子位置没有被移动；这仍是解析力中的密度响应计算，不是总能量有限差分。

## 3. POTCAR 表插值

CPU 与 CUDA 现在共用与 POTCAR 构造一致的规则：

- 局域势表使用 $G=0$ 一阶导数为零、末端二阶导数为零的三次样条；
- `PSPCOR` 和 `PSPRHO` 使用四点三次插值；
- CUDA device 端直接执行相同四点公式，不再为 core 表建立无用的二阶导数数组。

`half-energy-check` 增加了四点三次多项式精确性和局域势边界条件回归。

## 4. 正确的固定密度对照

普通 `ICHARG=1, NELM=1, ALGO=All` 在一次直接对角化后本征残差为零，VASP
会把该电子步标记为收敛，并不把 convergence correction 加入最终力。因此它
不是 Harris 响应项的正确对照。验证采用相同 CHGCAR、结构、k 点和 band 数，
但使用 `ICHARG=11` 明确固定密度，使 VASP 输出并加入该修正。

| 体系 | 旧 HALF MAE | 新 HALF MAE | 新最大分量差 | MAE 改善 |
|---|---:|---:|---:|---:|
| Si | `2.01117e-4` | `6.99076e-6` | `1.51799e-5` | `28.8x` |
| HfO2 | `3.64625e-2` | `1.00706e-3` | `2.57082e-3` | `36.2x` |

单位均为 eV/Angstrom。Hf 第一原子的 VASP convergence correction 为约
`(0.192, 0.0448, 0.0307)` eV/Angstrom；HALF 为
`(0.1906, 0.0437, 0.0299)` eV/Angstrom，说明 Hartree+XC 响应及
`PSPRHO` 收缩已逐项对齐。

最终 CUDA 结果：

| 体系 | HALF free energy (eV) | HALF–VASP Harris 力 MAE | 最大分量差 |
|---|---:|---:|---:|
| Si | `-43.3066227505` | `6.99076e-6` | `1.51799e-5` |
| HfO2 | `-121.8844407617` | `1.00706e-3` | `2.57082e-3` |

## 5. 回归与数据位置

CUDA direct-grid、CPU direct-grid 和 legacy shell 三种构建均通过 `9/9` CTest。

- 最终 HALF：`/data/limusen/deepaw_half_force_compare/final_vasp_interp_gpu`
- VASP 固定密度基准：
  `/data/limusen/mimic_us_energy_force_compare/{Si,HfO2}_harris_icharg11`
- L 通道扫描：`/data/limusen/deepaw_half_force_compare/l_channel_scan`
- 可机读结果：[`harris_force_l_response.json`](harris_force_l_response.json)
