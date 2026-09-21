# Harris 无矩阵全带加速原理与推导

本文说明 DeePAW-HALF CUDA `acc` 求解器及其 VASP 6.6 直接接入的理论、算法和
实现。目标是让数值过程可以审计：下述每个大规模运算都能在
`half_cuda_assembly.cuf` 或 `half_cuda_iterative_solver.cuf` 中找到对应实现。

## 1. 固定密度 Harris 问题

DeepAW 给出平滑价电子密度 $\widetilde{\rho}_0(\mathbf r)$。HALF 不更新这个密度，
而是只构造一次固定哈密顿量：

$$
\widehat H_0 = \widehat H[\widetilde{\rho}_0].
$$

初始轨道来自 PAW 广义本征问题：

$$
\widehat H_0 |\widetilde\psi_n\rangle
= \varepsilon_n \widehat S |\widetilde\psi_n\rangle.
$$

这与 VASP 的 SCF 过程不同。HALF 求解的是一次性 Harris 问题；VASP 把返回轨道
作为初始子空间，然后照常更新密度、优化轨道并完成自洽。

固定局域有效势为：

$$
V_{\mathrm{eff}}(\mathbf r)
= V_{\mathrm{ion}}^{\mathrm{loc}}(\mathbf r)
+ V_H[\widetilde\rho_0](\mathbf r)
+ V_{\mathrm{xc}}[\widetilde\rho_0+\rho_{\mathrm{core}}](\mathbf r).
$$

部分芯密度只进入 XC 的非线性芯修正，不会加入 Hartree 密度。

## 2. PAW 广义本征方程

从平滑赝波函数到全电子波函数的 PAW 变换为：

$$
\widehat{\mathcal T}
= 1 + \sum_{Ii}
\left(|\phi_i^I\rangle-|\widetilde\phi_i^I\rangle\right)
\langle\widetilde p_i^I|.
$$

因此赝空间的重叠算符不是单位算符：

$$
\widehat S
= \widehat{\mathcal T}^{\dagger}\widehat{\mathcal T}
= 1 + \sum_{Iij}
|\widetilde p_i^I\rangle Q_{ij}^I
\langle\widetilde p_j^I|.
$$

在平面波基 $|\mathbf k+\mathbf G\rangle$ 中定义投影矩阵：

$$
B_{\mathbf G,Ii}
= \langle\mathbf k+\mathbf G|\widetilde p_i^I\rangle.
$$

令 $T$ 为对角动能矩阵，$V$ 为局域势卷积，$D$ 和 $Q$ 为按原子分块的 PAW
矩阵。如果显式构造稠密矩阵，则有：

$$
H = T + V + B D B^{\dagger},
$$

$$
S = I + B Q B^{\dagger}.
$$

在 MIMIC_US 路径中，onsite 哈密顿量与原子有关：

$$
D_{ij}^{I}
= D_{ij}^{I,\mathrm{ION}}
+ \int V_{\mathrm{eff}}(\mathbf r)
Q_{ij}^{I,\mathrm{DEP}}(\mathbf r-\mathbf R_I)\,d^3r.
$$

`evd`、`evj` 和 `evx` 会显式构造 $H$ 与 $S$；`acc` 不构造任何
$N_{\mathrm{PL}}\times N_{\mathrm{PL}}$ 稠密矩阵。

## 3. 为什么必须采用全带 Rayleigh-Ritz

把所需的 $m=N_{\mathrm{BANDS}}$ 条轨道作为矩阵 $X$ 的列。固定密度低能子空间
可写为带约束的迹最小化：

$$
\min_X \operatorname{Tr}(X^{\dagger} H X)
\quad\text{且满足}\quad
X^{\dagger} S X = I_m.
$$

引入 Hermitian 拉格朗日乘子矩阵 $\Lambda$：

$$
\mathcal L(X,\Lambda)
= \operatorname{Tr}(X^{\dagger} H X)
- \operatorname{Tr}\left[\Lambda(X^{\dagger}S X-I_m)\right].
$$

对 $X^{\dagger}$ 求驻值得到：

$$
H X = S X \Lambda.
$$

所以所有能带通过共同变分子空间相互耦合。能带分块只是 FFT 和矩阵收缩的执行
单位，不能被当作彼此独立的小本征问题。这一点与 VASP `ALGO=All` 的本质一致：
执行过程可以分块，但优化和正交性必须覆盖全部目标能带。

## 4. 无矩阵 H 和 S 应用

对包含若干能带的输入块 $X_b$，HALF 计算：

$$
Y_b = B^{\dagger}X_b,
$$

$$
H X_b
= T X_b
+ \mathcal F\!\left[
V_{\mathrm{eff}}(\mathbf r)
\mathcal F^{-1}[X_b]
\right]
+ B D Y_b,
$$

$$
S X_b = X_b + B Q Y_b.
$$

$\mathcal F$ 与 $\mathcal F^{-1}$ 包含平面波表示采用的 FFT 归一化。实现中，
波函数规模的数据不会返回 CPU，步骤为：

1. 把平面波系数 scatter 到 FFT 网格；
2. 用 cuFFT 逆变换到实空间；
3. 用 CUDA kernel 乘 $V_{\mathrm{eff}}(\mathbf r)$；
4. 用 cuFFT 正变换并 gather 回平面波顺序；
5. 加入对角动能项；
6. 用 cuBLAS 和 CUDA 投影 kernel 完成 $B^{\dagger}X_b$、$BDY_b$、$BQY_b$。

`--acc-block-size` 只控制一次 H/S 应用处理多少条带，并不会切断不同执行块之间
的全带耦合。

## 5. 重启式分块迭代

### 5.1 初始子空间

HALF 取动能最低的 $m$ 个平面波单位向量作为初始态。VASP 任意 k 点上的索引通过
$O(N_{\mathrm{PL}}\log N_{\mathrm{PL}})$ 堆排序得到。原来的插入排序在大规模偏移
k 点上会退化到 $O(N_{\mathrm{PL}}^2)$，已经移除。

首次应用 H 和 S 后，通过一次 Rayleigh-Ritz 得到 S 正交的 Ritz 基。

### 5.2 广义残差

对 Ritz 向量 $X$ 和对角 Ritz 值 $\Lambda$，残差为：

$$
R = H X - S X\Lambda.
$$

程序报告的收敛量是：

$$
r_{\max} = \max_n \|R_n\|_2.
$$

由于 Ritz 向量已在 S 度量下归一化，$r_{\max}$ 具有能量单位，可以直接衡量
本征方程误差。`--acc-tol` 是目标值，`--acc-max-iter` 是硬上限；可能先达到
迭代上限。

### 5.3 动能预条件

平面波动能对角元可作为 $H-\varepsilon_nS$ 逆算符的低成本近似。HALF 使用：

$$
P_{\mathbf G n}
= -\frac{R_{\mathbf G n}}
{\left|T_{\mathbf G}-\varepsilon_n\right|+\delta},
$$

当前实现取 $\delta=1\ \mathrm{eV}$。它抑制残差中的高动能分量，而无需构造或
分解 H。

### 5.4 S 度量正交化

修正块必须在 PAW 的 S 度量下与当前 Ritz 空间正交。先计算：

$$
C = X^{\dagger} S P,
$$

然后同时更新修正块及已经计算好的算符像：

$$
P \leftarrow P-XC,
$$

$$
HP \leftarrow HP-(HX)C,
$$

$$
SP \leftarrow SP-(SX)C.
$$

再构造修正块度量：

$$
M=P^{\dagger}SP.
$$

若

$$
M=U\,\operatorname{diag}(\mu)\,U^{\dagger},
$$

则归一化修正块为：

$$
P \leftarrow P U\,\operatorname{diag}(\mu^{-1/2}).
$$

$HP$ 和 $SP$ 也乘同一个右变换。只有这个 $m\times m$ 度量本征求解使用
cuSOLVER。

### 5.5 全带 Rayleigh-Ritz 更新

构造重启试探空间及其算符像：

$$
Z=[X\;P],\qquad HZ=[HX\;HP],\qquad SZ=[SX\;SP].
$$

投影算符为：

$$
H_Z=Z^{\dagger}HZ,
$$

$$
S_Z=Z^{\dagger}SZ.
$$

用 cuSOLVER 求解小型广义本征问题：

$$
H_ZY=S_ZY\Theta.
$$

保留最低的 $m$ 列 $Y_m$，并对全部能带一起旋转：

$$
X\leftarrow ZY_m,\qquad
HX\leftarrow HZY_m,\qquad
SX\leftarrow SZY_m.
$$

子空间维度不超过 $2m$。因此 `acc` 模式下 cuSOLVER 只处理 $m\times m$ 或
$2m\times2m$ 小矩阵，从不会接收 $N_{\mathrm{PL}}\times N_{\mathrm{PL}}$ 稠密
哈密顿量。

## 6. 面向 VASP 初始波函数的提前停止

独立能带计算需要严格收敛的本征对；VASP 初始化只需要一个与低能本征空间具有
足够重叠的子空间。VASP 随后还会自行正交化、进行全带优化、更新密度并完成 SCF。

因此 HALF 允许采用较宽松的残差阈值，也允许达到迭代上限后返回。这是明确的
精度—速度折中，并不表示返回态已经是最终 SCF 轨道。在已验证的 20 原子
CsPbBr3 体系中，40 次 ACC 迭代后各 k 点最大残差为 0.022--0.410 eV；但 ACC
与稠密 VASP_LIKE 都只需要 5 个 VASP SCF LOOP，并得到相同的最终能量。

## 7. 复杂度和显存

记 $N=N_{\mathrm{PL}}$、$m=N_{\mathrm{BANDS}}$、$N_r$ 为实空间 FFT 网格规模、
$b$ 为 H/S 执行块大小、$N_p$ 为 PAW 投影子总数。

稠密路径至少需要保存 H 和 S：

$$
M_{\mathrm{dense}} \ge 2\times16N^2\ \mathrm{bytes}
=32N^2\ \mathrm{bytes},
$$

这里还没有包含求解器工作区和副本。稠密组装为 $O(N^2)$，完整广义本征求解为
$O(N^3)$；EVX 区间求解能减少本征矢工作量，但不能消除稠密 H/S 的存储和组装。

ACC 的态数据存储随平面波数线性增长：

$$
M_{\mathrm{ACC}}=O(Nm+N_rb+NN_p+m^2).
$$

一次无矩阵 H/S 应用的近似复杂度为：

$$
O\!\left(mN_r\log N_r + NmN_p + Nm\right).
$$

投影 cuBLAS 乘法为 $O(Nm^2)$，小型 Rayleigh-Ritz 为 $O(m^3)$。因此主要标度
从完整基组维度转移到了实际请求的能带数。

CsPbBr3 验证中 $N=20640$--$20780$，$m=105$。在 $N=20767$ 时，仅两份稠密
complex128 H/S 就需要约 13.8 GB，而 ACC 完全不分配它们。

## 8. VASP 直接内存接入

适配后的 VASP 源码保存在 `vendor/vasp-6.6.0`。调用流程为：

1. VASP 读取 INCAR、POSCAR、POTCAR 和 DeepAW 平滑 CHGCAR；
2. `ALLOCW` 后，VASP 把晶格、分数坐标、物种、FFT 网格、ENCUT 和 NBANDS 直接
   传给 libhalf；
3. 对每个 VASP k 点，adapter 从 `wavedes1` 构造严格一致的串行
   $(G_x,G_y,G_z)$ 表；
4. `half_solve_kpoint_mapped()` 验证两个基组完全相同，执行 ACC，并按 VASP 顺序
   返回平面波系数；
5. VASP 用 `DIS_PW_BAND` 分发每条全局能带，把 Ritz 值写入 `W%CELTOT`，跳过随机
   `WFINIT`，然后继续正常的 `PROALL`、`ORTHCH` 和 SCF 流程。

对应 INCAR 设置为：

```text
LHALF_INIT = .TRUE.
LHALF_API  = .FALSE.
HALF_MODE  = ACC
```

`LHALF_API=.FALSE.` 使用现有 DeepAW CHGCAR；设为 `.TRUE.` 时通过
`HALF_ESCN_URL` 选择 DeePAW-eSCN 密度 API，后续波函数内存传递过程不变。

当前 adapter 要求 ISPIN=1、无 SOC/非共线自旋、使用 full-complex VASP，并设置
KPAR=1。这些是 adapter 接线限制，不是无矩阵方程本身的限制。

## 9. 实现位置

| 运算 | 源码 |
|---|---|
| GPU 固定密度有效势 | `src/half_cuda_potential.cuf` |
| PAW/MIMIC_US 数据与无矩阵 H/S | `src/half_cuda_assembly.cuf` |
| 残差、预条件、S 正交化、重启 Ritz 求解 | `src/half_cuda_iterative_solver.cuf` |
| 稠密/ACC 分发 | `src/half_cuda_solver.cuf` |
| C 和 Fortran 求解器选择器 | `include/half.h`、`include/half_api.f90` |
| CLI 控制和 JSON 诊断 | `app/half_cli.F90` |
| VASP 内存 adapter | `vendor/vasp-6.6.0/src/half_vasp_init.F` |
| VASP INCAR 验证 | `vendor/vasp-6.6.0/src/reader.F` |

Pro 6000 实测结果和完整计算条件记录在
[`validation/cspbbr3_vasp_acc_pro6000.json`](validation/cspbbr3_vasp_acc_pro6000.json)。
