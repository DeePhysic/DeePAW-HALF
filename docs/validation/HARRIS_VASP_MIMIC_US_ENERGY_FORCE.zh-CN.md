# HALF Harris 与 VASP MIMIC_US 能量、力对比

> 历史参照说明：本文的普通 `ICHARG=1`、`NELM=1` VASP 力在一次精确
> 对角化后未包含固定密度 Harris 收敛修正，因而不能作为最终解析力 oracle。
> 当前实现与严格 `ICHARG=11` 参照、完整 Hartree+XC 响应及 $L$ 通道诊断见
> [`HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.zh-CN.md`](HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.zh-CN.md)。

## 1. 对比定义

本报告中的 `MIMIC_US` 严格指 VASP 的

```text
LMAXMIX = -1
```

模式，而不是某个电子迭代步的名称。为了避免拿收敛 SCF 结果与固定密度
Harris 结果混比，VASP 计算使用 `NELM=1`、`NELMIN=1`，只截取从同一 DeepAW
平滑密度及 HALF 初始波函数出发的未收敛 MIMIC_US 状态。

共同设置为 PBE、`ISPIN=1`、`ISMEAR=0`、`SIGMA=0.02`、
`KSPACING=0.35`、`KGAMMA=.TRUE.`、`LHALF_INIT=.TRUE.`、
`LHALF_API=.FALSE.`、`LMAXMIX=-1`。Si 使用 `ENCUT=520 eV`，HfO2 使用
`ENCUT=521 eV`。

CHGCAR 只提供 `NGX*NGY*NGZ` 个平滑密度值。HALF 不读取 CHGCAR 文件尾部的
augmentation occupancies。所有 onsite occupancy、AE-PS multipoles、
compensation functions 与 augmentation density 均由当前波函数和匹配 POTCAR
重建。

## 2. 能量结果

### 2.1 POTCAR 自动重建后的绝对能量

| 体系 | HALF Harris free energy (eV) | VASP MIMIC_US `NELM=1` TOTEN (eV) | HALF - VASP (eV) |
|---|---:|---:|---:|
| Si | -43.314024964 | -43.293318230 | -0.020706734 |
| HfO2 | -121.886526811 | -121.972582920 | 0.086056109 |

HALF 现在直接用 POTCAR 的原子占据矩阵、AE/PS partial waves、core density、
`DEXC` 和补偿电荷重建球形原子 PAW double counting，不再要求调用者传入
VASP 输出中的数值，也不读取 CHGCAR augmentation 尾部。命令行中的
`--paw-atomic-double-counting` 只保留为显式覆盖/诊断入口。

自动重建项与 VASP 当前 MIMIC_US 状态为：

| 体系 | HALF POTCAR 重建 (eV) | VASP (eV) | HALF - VASP (eV) |
|---|---:|---:|---:|
| Si | 140.758772710 | 139.835990310 | 0.922782400 |
| HfO2 | 36.089365550 | 27.477902270 | 8.611463280 |

AE 与 PS 原子项各自很大并强烈抵消；当前径向 PBE 实现尚未逐分项达到 VASP
精度。不过，将该 POTCAR 重建项纳入完整 Harris 泛函后，最终总能量误差已经
降到 Si `0.021 eV/cell`、HfO2 `0.086 eV/cell`。这不是用 VASP 数值标定得到的。

重建从 POTCAR 的球形原子占据 `QATO` 开始。对角角动量通道的价电子径向密度为

$$
\widetilde n^{\mathrm{AE/PS}}_{00}(r)=
Y_{00}\sum_{ab}\delta_{l_a l_b}(2l_a+1)
Q^{\mathrm{atom}}_{ab}
W^{\mathrm{AE/PS}}_a(r)W^{\mathrm{AE/PS}}_b(r).
$$

PS 侧再加入由 AE-PS 多极矩确定的球形补偿电荷 $\widehat n_{00}(r)$，并检查

$$
\int n^{\mathrm{AE}}(\mathbf r)\,d\mathbf r
=\int\left[n^{\mathrm{PS}}(\mathbf r)+\widehat n(\mathbf r)\right]d\mathbf r.
$$

当前实现计算的两个位置无关原子项为

$$
E_{\mathrm{dc}}^{\mathrm{AE}}=
-E_H[n^{\mathrm{AE}}]
+E_{\mathrm{xc}}[n^{\mathrm{AE}}+n_c^{\mathrm{AE}}]
-\int n^{\mathrm{AE}}v_{\mathrm{xc}}[n^{\mathrm{AE}}+n_c^{\mathrm{AE}}]d\mathbf r
-E_{\mathrm{xc}}[n_c^{\mathrm{AE}}]-E_{\mathrm{DEXC}},
$$

$$
E_{\mathrm{dc}}^{\mathrm{PS}}=
+E_H[n^{\mathrm{PS}}+\widehat n]
-E_{\mathrm{xc}}[n^{\mathrm{PS}}+\widehat n+n_c^{\mathrm{PS}}]
+\int(n^{\mathrm{PS}}+\widehat n)
v_{\mathrm{xc}}[n^{\mathrm{PS}}+\widehat n+n_c^{\mathrm{PS}}]d\mathbf r.
$$

整胞的 `paw_atomic_double_counting_eV` 是各元素
$E_{\mathrm{dc}}^{\mathrm{AE}}+E_{\mathrm{dc}}^{\mathrm{PS}}$ 乘原子数后的和。

### 2.2 分量对比

| 体系/分量 | HALF (eV) | VASP 未收敛 MIMIC_US (eV) | HALF - VASP (eV) |
|---|---:|---:|---:|
| Si band | 22.223275800 | 23.492086290 | -1.268810490 |
| Si Hartree double counting | -60.273397505 | -60.602694260 | 0.329296755 |
| Si combined XC | -71.960292082 | -71.951882680 | -0.008409402 |
| Si Ewald | -911.952813534 | -911.952858480 | 0.000044946 |
| Si local G=0 / PSCENC | 13.355229652 | 13.355229650 | 0.000000002 |
| Si atomic reference / EATOM | 824.535200000 | 824.530810940 | 0.004389060 |
| HfO2 band | -816.594587937 | -801.166613980 | -15.427973957 |
| HfO2 Hartree double counting | -1369.342375420 | -1376.560388480 | 7.218013060 |
| HfO2 combined XC | 314.708708971 | 315.196230200 | -0.487521229 |
| HfO2 Ewald | -5797.173198768 | -5797.174463800 | 0.001265032 |
| HfO2 local G=0 / PSCENC | 607.712760793 | 607.712832880 | -0.000072087 |
| HfO2 atomic reference / EATOM | 6902.712800000 | 6902.541917990 | 0.170882010 |

在实现 POTCAR 自动重建之前，把 VASP 当前步打印的 PAW double counting 与
VASP `EATOM` 临时补到旧 HALF 结果后，残差曾为：

| 体系 | 临时补项后的 HALF - VASP (eV) |
|---|---:|
| Si | -0.947878191 |
| HfO2 | -8.696289181 |

该旧诊断说明问题不只是一个位置无关常数。当前自动重建显著改善了最终总能量，
但 band、Hartree、径向 PBE 与 onsite 装配仍未逐项完全等同于 VASP。

## 3. 力结果

VASP 力取上述 `LMAXMIX=-1`、`NELM=1` 计算结束时写出的
`TOTAL-FORCE`，不是收敛 SCF 力。

| 体系 | HALF 最大绝对分量 (eV/Angstrom) | VASP 最大绝对分量 (eV/Angstrom) | 分量 MAE | 分量 RMSE | 最大分量差 |
|---|---:|---:|---:|---:|---:|
| Si | 0.006975088 | 0.000858000 | 0.005320861 | 0.005455295 | 0.007808088 |
| HfO2 | 0.104255268 | 0.215302000 | 0.068045473 | 0.115119451 | 0.319542268 |

这里已修正固定密度 Harris 响应：DeepAW 平滑输入密度是外部固定场，不能用
POTCAR `PSPRHO` 近似成随原子移动的输入密度。新实现只保留 POTCAR core density
移动引起的 XC-kernel 响应，不加入 Hartree 响应。该修正使 Si 最大力从
`0.0499` 降至 `0.0070 eV/Angstrom`，HfO2 从 `0.2940` 降至
`0.1043 eV/Angstrom`。

为了区分“力公式错误”和“HALF/VASP 泛函差异”，又对同一个 HALF Harris 能量
做了不进入生产路径的中心差分：

| 体系/分量 | HALF 解析力 | HALF 能量中心差分 | 绝对差 (eV/Angstrom) |
|---|---:|---:|---:|
| Si atom 1, x | 0.006100663 | 0.004351093 | 0.001749570 |
| HfO2 atom 1, x | -0.104255268 | -0.107329684 | 0.003074416 |

HfO2 的解析导数与自身能量同号且接近，而 VASP 同一分量为
`+0.215287 eV/Angstrom`。这表明当前最大的 VASP 差异来自 HALF 与 VASP 的
onsite PAW/QDEP/one-centre 泛函尚未完全同构，而不是简单漏掉总能量负梯度的
符号。解析力内部仍有 `1e-3 eV/Angstrom` 量级的离散不一致需要继续收敛。

augmentation 电子数也存在可测差异：

| 体系 | HALF POTCAR 重建 augmentation (e) | VASP MIMIC_US (e) | 差值 (e) |
|---|---:|---:|---:|
| Si | -1.84190385 | -1.83255990 | -0.00934395 |
| HfO2 | 16.06675370 | 16.01618910 | 0.05056460 |

## 4. 结论

1. Ewald 和局域势零频修正已经与 VASP 对齐。
2. PAW atomic double counting 已由 POTCAR 自动重建；CHGCAR augmentation
   尾部不参与计算。原子 AE/PS 大数抵消的逐分项精度仍需继续提高。
3. POTCAR augmentation density 已进入 HALF 力，固定密度 Harris 响应也已改为
   core-only XC 响应，但 onsite、band/Hartree/XC 与 VASP MIMIC_US 尚未闭合。
4. 当前 HALF 解析力不能宣称复现 VASP MIMIC_US 力；Si 的对称性允许力应接近
   零，而 HALF 仍出现约 `0.007 eV/Angstrom`，这是明确的待修问题。
5. 下一步应以同一 HALF 泛函的中心差分为第一层 oracle，再用 VASP
   `SETDIJ`、augmentation occupancy、PAW double counting 和逐项力闭合第二层
   VASP parity。

## 5. 如何继续提高力精度

1. **保持变分一致。** 当前不含 augmentation 的 Si 解析力与相同 HALF 能量的
   中心差分只差约 `1.5e-4 eV/Angstrom`；主要剩余误差集中在
   `Qij*Veff`。下一步要让 augmentation density 对局域势位移的积分与构造
   `Dij` 时使用完全相同的径向、角向和周期三次 B-spline 求积，避免“Dij 用球面
   求积、局域力用 FFT 网格采样”的离散不一致。
2. **闭合 onsite 能量与导数。** 用同一个 onsite density matrix 同时产生
   augmentation charge、`Dij`、AE/PS one-centre energy 及其导数，避免多个路径
   独立重建同一物理量。
3. **提高径向/角向求积。** 对 Hf 等高角动量 PAW 数据提高 Gauss-Legendre 与
   方位角阶数，并以阶数加倍后的力变化作为收敛判据；径向 PBE 梯度使用 POTCAR
   对数网格上的解析/高阶离散，而不是低阶差分。
4. **先验证原始力，再施加对称性。** 对每个对称不等价原子和三个笛卡尔方向做
   `h=0.002, 0.001, 0.0005 Angstrom` 三点收敛；原始解析导数闭合后，再使用空间群
   投影和 acoustic sum rule 去除纯数值噪声，不能用对称化掩盖公式误差。
5. **全程 complex128。** projector、occupancy、Rayleigh quotient 和力归约保持
   complex128/float64，并使用确定性归约检查 GPU 并行求和噪声。

## 6. 原始文件

- HALF Si：`/data/limusen/deepaw_half_force_compare/Si/half_forces_core_response.json`
- HALF HfO2：`/data/limusen/deepaw_half_force_compare/HfO2/half_forces_core_response.json`
- VASP Si：`/data/limusen/mimic_us_energy_force_compare/Si/OUTCAR`
- VASP HfO2：`/data/limusen/mimic_us_energy_force_compare/HfO2/OUTCAR`
