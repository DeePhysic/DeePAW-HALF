# onsite PAW / one-centre 构造对齐与优化

本文记录 `force-accuracy-experiments` 分支对 VASP 6.6.0 `MIMIC_US`
构造路径的逐项核对和实现。这里只描述算法与验证结果，不包含专有源码。

## 1. VASP 路径与原 HALF 路径的差异

在 `LMAXMIX=-1` 的固定密度计算中，密度相关的完整 PAW one-centre
泛函被关闭；保留下来的 onsite 广义本征问题为

$$
H=H_{\mathrm{smooth}}+
\sum_{Iij}|\widetilde p_i^I\rangle D_{ij}^I
\langle\widetilde p_j^I|,
\qquad
S=1+\sum_{Iij}|\widetilde p_i^I\rangle q_{ij}^I
\langle\widetilde p_j^I|.
$$

其中

$$
D_{ij}^I=D_{ij}^{\mathrm{ion}}+
\sum_{LM}C_{ij}^{LM}Q_{ij}^{L}D_{LM}^{I},
$$

$$
D_{LM}^{I}=\int V_{\mathrm{eff}}(\mathbf r)
g_L(|\mathbf r-\mathbf R_I|)
Y_{LM}(\widehat{\mathbf r-\mathbf R_I})\,d\mathbf r.
$$

VASP 把归一化补偿函数 $g_LY_{LM}$ 直接放到离子附近的 FFT 网格点上。
原 HALF 则先用三次样条把三维势插值到径向球壳，再做径向和角向求积。
后者增加了与 FFT 网格无关的第二套离散，也使 $D$ 与其位置导数更难严格一致。

## 2. 直接 FFT 网格 SETDIJ

新实现使用与 VASP 相同的离散目标：

$$
D_{LM}^{I}\approx\frac{\Omega}{N_{\mathrm{grid}}}
\sum_{\mathbf r_g}V_{\mathrm{eff}}(\mathbf r_g)
g_L(r_{gI})Y_{LM}(\widehat{\mathbf r_{gI}}),
\qquad r_{gI}\le r_c.
$$

补偿函数仍是两个球 Bessel 函数的线性组合：

$$
g_L(r)=A_1j_L(q_1r)+A_2j_L(q_2r),
$$

其中 $j_L(q_nr_c)=0$，并同时满足边界一阶导数为零与

$$
\int_0^{r_c}g_L(r)r^{L+2}\,dr=1.
$$

CPU 和 CUDA Fortran 后端均实现了直接网格路径。CUDA 使用一个
`(ion,LM)` 对应一个 thread block 的归约，势、补偿核和 $D_{LM}$ 全程留在 GPU。
旧球壳算法保留为 CMake A/B 开关：

```bash
-DHALF_QDEP_DIRECT_GRID=ON   # 默认
-DHALF_QDEP_DIRECT_GRID=OFF  # 旧球壳积分，仅用于回归比较
```

## 3. one-centre 径向矩与 overlap

VASP 在读入 partial waves 后，从 POTCAR 的对数径向网格重新构造

$$
Q_{ij}^{L}=\int_0^{r_c}
\left[\phi_i(r)\phi_j(r)-\widetilde\phi_i(r)\widetilde\phi_j(r)\right]
r^L\,dr.
$$

新实现不再在不同模块中分别重复积分，而是在 POTCAR 解析完成时建立
`qpaw_l(i,j,L)`。积分权重严格采用对数坐标 $x=\ln r$ 上的 Simpson 规则，
$dr=r\,dx$。随后：

- overlap 使用 $q_{ij}=Q_{ij}^{0}$；
- `SETDIJ` 使用相同的 $Q_{ij}^{L}$；
- augmentation density 使用相同的 $Q_{ij}^{L}$；
- 原子 PAW double-counting 使用相同的 $Q_{ij}^{0}$。

因此 $H$、$S$、能量和力不再各自拥有略有差异的径向矩定义。Si 和 HfO2
POTCAR 中，重构后的有效同角动量 $Q^0$ 与文件表值最大差分别仅为
`2.48e-7` 和 `3.71e-7`，说明输入一致，但显式重构消除了后续漂移。

## 4. 解析导数

位置导数作用在同一个离散补偿核上：

$$
\frac{\partial D_{LM}^{I}}{\partial R_{I\alpha}}=
\frac{\Omega}{N_{\mathrm{grid}}}
\sum_{\mathbf r_g}V_{\mathrm{eff}}(\mathbf r_g)
\frac{\partial}{\partial R_{I\alpha}}
\left[g_L(r_{gI})Y_{LM}(\widehat{\mathbf r_{gI}})\right].
$$

这不是移动原子后重新计算总能量；它是单次 Hamiltonian 中紧支撑 QDEP 核的
局部导数。`half-force-check` 再用原子位移中心差分作为独立 oracle。结果为：

| 体系 | `max |analytic dD/dR - FD dD/dR|` (eV/Angstrom) |
|---|---:|
| Si | `2.75e-9` |
| HfO2 | `2.43e-9` |

相较 `24×48` 球壳积分的约 `4.36e-5`（Si）与 `2.56e-4`
（HfO2），离散一致性提高约四到五个数量级。

## 5. 完整能量/力结果

Pro 6000、CUDA 13.2、NVHPC 26.5、PBE、`ISPIN=1`、`SIGMA=0.02 eV`，
使用 VASP `IBZKPT`；Si 为 520 eV/24 bands，HfO2 为 521 eV/80 bands。

| 体系 | HALF free energy (eV) | HALF-VASP (eV/cell) | 力分量 MAE (eV/Angstrom) | 最大分量差 |
|---|---:|---:|---:|---:|
| Si | `-43.306625404` | `-0.013307174` | `0.000739611` | `0.001308693` |
| HfO2 | `-121.884427193` | `+0.088155727` | `0.068379790` | `0.319947241` |

直接网格使 Si 相对 VASP 的能量差由约 `20.7 meV/cell` 降到
`13.3 meV/cell`，并把 HALF 最大残余力降到 `0.457 meV/Angstrom`。
HfO2 的内部 $D$ 导数已经达到数值精度，但相对 VASP 的力差基本不变。
这证明 HfO2 的剩余误差不再来自 QDEP 网格、角向积分或 $Q^L$ 径向矩；下一步
应逐项核对原子 AE/PS double-counting 和输出 augmentation density 对局域力的
贡献，不能再靠增加球面积分点解决。

解析力阶段还缓存了与 k 点无关的 `D_ij` 和 `dD_ij/dR`。HfO2 完整测试
wall time 从 `121.87 s` 降到 `97.78 s`；该数字包含本征求解且测试时 GPU
存在其他进程，只说明消除了重复构造，不作为独占 GPU 性能结论。

## 6. 验证位置

- 直接网格结果：`/data/limusen/deepaw_half_force_compare/direct_grid_setdij`
- VASP 对数权重结果：`/data/limusen/deepaw_half_force_compare/direct_grid_vasp_weights`
- 构建：`/data/limusen/deepaw-half-acc-src/build/direct-grid-cpu` 与
  `/data/limusen/deepaw-half-acc-src/build/direct-grid-cuda-nvhpc`

