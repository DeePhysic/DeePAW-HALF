# 解析力的不可约 k 点空间群展开

本文记录 `force-kpoint-symmetry` 分支实现的解析力空间群展开。它使独立 HALF
CLI 能够在不破坏晶体对称性的情况下，只求解不可约 k 点。

## 方法

对每个不可约代表点 $\mathbf k$，HALF 只对角化一次 PAW Hamiltonian。该点的
权重已经包含完整 k-star；对无磁标量计算还包含时间反演伙伴。标量密度按空间群
操作 $g=(R_g,\tau_g)$ 展开：

$$
\overline n_{\mathbf k}(\mathbf r)=
\frac{1}{|G|}\sum_{g\in G}
n_{\mathbf k}\!\left(R_g^{-1}(\mathbf r-\tau_g)\right).
$$

HALF 先从代表点波函数及 onsite occupancy 重建平滑密度和 POTCAR PAW
augmentation density，再对两者之和做空间群展开。这与逐个旋转 star 中的
波函数和 PAW occupancy 等价，但不需要保存完整网格波函数或显式 Wigner 矩阵。

当操作 $g$ 把原子 $I$ 映射到 $gI$ 时，Cartesian 力按下式旋转：

$$
\overline{\mathbf F}_{gI}\mathrel{+}=
\frac{1}{|G|}\,C_g\mathbf F_I,
\qquad
C_g=L^{T}R_g(L^{T})^{-1}.
$$

Ewald、局域势、projector、PAW augmentation、NLCC 与 Harris 响应力分别使用
相同的原子置换和 Cartesian 旋转。密度和无磁原子力在时间反演下为偶，因此
时间反演只由 k 点 multiplicity 处理，不需要额外变换。

空间群旋转和平移由 spglib 提供。选择不可约网格时，输入 DeepAW 密度也会在
构造 $V_{\mathrm{eff}}$ 前使用同一空间群对称化，保证 Hamiltonian 与 k 点约化
采用一致的对称性。

## Si 原胞验证

CLI 测试采用二原子 diamond-Si 原胞、裸机 RTX PRO 6000、PBE、520 eV、
`KSPACING=0.35`、12 bands、CUDA EVD 和单 GPU。

| 指标 | 完整网格 | 不可约点 + 空间群展开 |
|---|---:|---:|
| 实际求解 k 点 | 216 | 16 |
| 空间群操作数 | 1 | 48 |
| free energy (eV/cell) | `-10.8288461353` | `-10.8288461451` |
| 能量差 | 参照 | `-9.78e-9 eV` |
| 最大力模 | `1.65e-7 eV/Angstrom` | `2.45e-18 eV/Angstrom` |
| 最大力分量差 | 参照 | `1.06e-7 eV/Angstrom` |
| wall time | `35.41 s` | `4.42 s` |

力的差异处于完整网格自身的数值噪声范围。不可约实现加速 **8.01 倍**，需要的
对角化次数减少 13.5 倍。

JSON 输出增加 `space_group_operation_count`；启用该路径时，
`force_kpoint_expansion` 为 `observable-space-group`。

## 回归

CUDA direct-grid、CPU direct-grid 和 legacy-shell 三种构建均通过 `9/9` CTest。
可机读结果见
[`kpoint_symmetry_force_expansion.json`](kpoint_symmetry_force_expansion.json)，
原始运行位于
`/data/limusen/deepaw_half_si_primitive_cli_timing_20260921`。
