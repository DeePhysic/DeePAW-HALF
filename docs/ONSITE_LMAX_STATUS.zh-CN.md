# HALF 运行时 onsite 通道：实现状态与物理契约

`--onsite-lmax` 与 VASP 的 `LMAXMIX` 是不同控制量。HALF 的默认值 `-1`
仍执行已有的固定输入平滑密度、MIMIC_US 型 Harris 路径。非负值必须启用
one-centre PAW 泛函，不能用 `HALF_QDEP_LMAX` 截断补偿电荷来冒充。

目前 `half bands` 和 `half energy` 均解析该参数。默认构建对非负值报错；
仅在显式设置 `-DHALF_EXPERIMENTAL_ONSITE=ON` 时开放实验路径。它不是
已经完成验证的生产功能，不能据此宣称 BN 高 $l$ EOS 已解决。

`half_onsite_density` 从 POTCAR 的 `QATO` 构建原子参考
onsite 占据矩阵，并对任意指定非负 $L_{\max}$，由当前占据矩阵 $P^I_{ab}$
构建 AE、PS 的径向多极密度与补偿矩：

$$
\rho^{I,\mathrm{AE/PS}}_{LM}(r)
=\rho^{I,\mathrm{core}}_{00}(r)\,\delta_{L0}\delta_{M0}
+\sum_{ab}P^I_{ab}\,G^{LM}_{ab}
  w^{I,\mathrm{AE/PS}}_a(r)w^{I,\mathrm{AE/PS}}_b(r),
$$

$$
Q^I_{LM}=\sum_{ab}P^I_{ab}\,G^{LM}_{ab}\,Q^L_{ab},
\qquad
G^{LM}_{ab}=\int Y_a(\hat r)Y_b(\hat r)Y_{LM}(\hat r)\,d\Omega.
$$

这里的径向量采用 POTCAR 的 $r^2\rho_{LM}(r)$ 约定。
令 $P_0$ 是 POTCAR `QATO` 给出的原子参考占据，$F_L(P)$ 是截至
$L$ 的 AE 减 PS one-centre Hartree/PBE 泛函（PS 侧含补偿电荷，XC
侧含各自的芯密度）。为使 POTCAR 已包含的参考一阶势不重复进入哈密顿量，
实验路径使用

$$
\Delta F_L(P)=F_L(P)-F_L(P_0)
-\left.\frac{\partial F_L}{\partial P}\right|_{P_0}:(P-P_0),
$$

$$
\Delta D_L(P)=\frac{\partial F_L}{\partial P}(P)
-\frac{\partial F_L}{\partial P}(P_0),
\qquad
\Delta E_{\mathrm{dc}}=\Delta F_L(P)-P:\Delta D_L(P).
$$

这里冒号是占据矩阵的全指标收缩。$\Delta D_L$ 加到每个原子的
非局域 PAW 投影算子；$\Delta E_{\mathrm{dc}}$ 加到本征值求和后的
总能量表达式。这样在 $P=P_0$ 时，两项都精确归零。

单元测试覆盖 s 通道球形密度、p 通道四极矩、`Lmax=0` 截断、原子参考占据，
以及 one-centre 能量对占据的数值导数。

实验路径已从本征波函数计算 onsite 占据，自洽迭代占据依赖的 `Dij`，
并将 one-centre Hartree/PBE 非线性残差及双计数修正接入能量。CPU 与 CUDA
求解器都接受该 `Dij`；CPU MPI 的占据矩阵也执行全局求和。CUDA Fortran
已在 GPU 上执行投影子–波函数收缩、占据矩阵构造、径向 AE/PS 多极密度、
Hartree/PBE 求积及占据导数的多次能量求值。POTCAR 的静态 Gaunt 系数、
补偿函数系数准备，以及最终解析力组合仍在主机侧，因此**尚不能称整个
HALF 工作流全 GPU**。

最终的大核稳定差分步长下，BN 粗网格单进程与两进程 MPI 的
`l=0` 自由能分别为 `14.657708959716214` 和
`14.657708946668853 eV`，差约 `13.0 neV`。
先前将 CPU one-centre 泛函接入 CUDA 求解器的过渡版本中，同一输入的
`l=0` 自由能 CPU/GPU 差 `2.79e-8 eV`、解析力最大分量差
`2.46e-7 eV/Å`；`l=2` 能量差 `5.93e-8 eV`，三点显式能带路径的
最大本征值差为 `2.74e-9 eV`。当前径向求积也搬上 GPU 后，BN
`l=0/2` 自由能差约 `3e-7 eV`；重测后的 `l=0` 解析力最大分量差
在最终大核稳定差分步长下为 `5.93e-7 eV/Å`，三点显式能带最大本征值差
`5.70e-8 eV`。
这些 GPU0 实测在 EOS 批次旁进行，没有使用 GPU1；不应用其运行时间
比较性能。
GPU 径向泛函的独立 CPU 对照（相同预设占据矩阵）在 TiO₂ `l=2/4`
下通过；增大中心差分占据步长以压低重元素大能量相消噪声后，TiO₂
`l=4` 的能量、双计数、最大 `Dij` 差分别为 `4.42e-8 eV`、
`5.34e-8 eV`、`3.37e-8 eV`，CeO₂ `l=6` 对应为
`2.13e-7 eV`、`1.04e-6 eV`、`2.56e-7 eV`。这比仅比较完整本征
求解更直接地检查了含 d/f 项的 CUDA one-centre 泛函。此前步长较小
的 TiO₂ 完整固定点 CPU/GPU 自由能差约 `3.58e-5 eV`；新步长的
固定点对照为 `5.69e-6 eV`，两端均在第 17 次达到 `Dij` 变化
小于 `1e-6 eV`。仍需在生产级网格上验证。
TiO₂ 的 120 eV、2.0 Å⁻¹、28 bands 粗网格最终代码试算在 `-1/0/2/4`
下分别为 `73.733226638/76.843113991/77.016308956/77.095052053 eV`；
`l=2` 的 `Dij`
自洽在旧的固定 0.5 混合下第 45 次达到 `9.48e-8 eV`，改用带增残回退的
自适应更新后第 20 次达到 `5.46e-8 eV`；同设置的 CUDA/CPU 自由能
相差 `2.98e-7 eV`，两者均在第 20 次收敛。这些只证明通道有实质响应和流程可运行，
`l=4` 比 `l=2` 高约 `0.0787431 eV`。这些结果不构成 VASP 数值一致性或
EOS 准确性验证。

生产模式开放前还必须完成：

1. 对照 VASP 的完整 AE/PS one-centre 路径，逐项确认当前相对 QATO 的
   Hartree/PBE 非线性残差与 POTCAR `DION`、原子参考能量严格相容；特别是
   离子、动能与补偿项不能重复计算或漏算。
2. 对 Si、BN、含 d/f 投影子体系核验 onsite 占据、`Dij`、能量与能带，
   并确认不同 $l$ 的物理效应和 `-1` 的严格回归。
3. 将 POTCAR 静态准备和力组合也迁入 GPU，并减少当前 CUDA 径向泛函
   的多次有限差分求值；与 CPU 在实际自洽占据固定点上达到更严格的一致性。
4. 与该算子一致的解析力，包括投影子位置导数与 one-centre 响应；
   用有限差分及 VASP 的对应 PAW 路径分别验证。在 BN 的 200 eV、
   1.5 Å⁻¹ 粗网格上，`l=0` 减 `-1` 的解析力增量为
   `-0.006258032 eV/Å`，相同两套能量的中心差分增量为
   `-0.006271455 eV/Å`，相差 `1.34e-5 eV/Å`。这是早期差分步长下 onsite *增量*
   的内部导数检验；两套绝对力各自对其有限差分仍相差约
   `0.164 eV/Å`，主要是原有 Harris 路径的共同差距，不能把增量检验
   表述为绝对力已通过。

最终差分步长 `1e-3`、同一 200 eV/1.5 Å⁻¹ 输入下，CUDA 解析力
onsite 增量为 `-0.006256910 eV/Å`，能量中心差分增量为
`-0.006279059 eV/Å`，相差 `2.21e-5 eV/Å`。这仍是增量验证，
并未改变绝对力尚待核验的结论。
同一脚本的 NVHPC CPU 和 ifx CPU 增量误差分别为
`3.21e-5`、`3.41e-5 eV/Å`；逐项数据见
`docs/validation/onsite_force_increment_bn_final_step.json`。

只有这些项在 CPU 和 CUDA 两套路径都完成且通过数值验证后，才移除
非负值的安全报错。
