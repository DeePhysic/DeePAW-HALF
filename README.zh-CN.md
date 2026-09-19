# DeePAW–HALF

**HALF**（**H**arris **A**ssociative **L**inearized **A**ugmented **P**lane
**W**ave **F**ortran）是 DeePAW–HAPPY 的 CPU/CUDA Fortran 固定密度电子结构
重建实现。它不是新的机器学习哈密顿量：输入是 DeepAW 或 VASP 产生的平滑
`CHGCAR` 密度以及匹配的 PAW `POTCAR`，输出为固定密度下的本征值与验证报告。

```text
结构 -> DeepAW -> 平滑 CHGCAR -> HALF -> H(k), S(k), 能带/本征值
```

HAPPY 是移植过程中的数值 oracle。当前 CUDA Gamma/MIMIC_US 固定密度路径已在
Si 与 HfO2 上通过逐本征值验证；CPU QDEP、任意 k 点、多 k 点能带、spglib
不可约 k 网格、占据数、Ewald 和 Harris 固定密度总能量也已闭合。Si 总能量及
各分量与 HAPPY 的差小于 `4e-12 eV`；中心差分力在 Si 验证中与 HAPPY 相差小于
`2.3e-9 eV/Angstrom`。原生 `vaspwave.h5` 输出以及 HDF5
电荷/结构/内嵌 POTCAR 输入已经支持；力仍在移植中。

## 数值模型：从固定密度到本征值

HAPPY 的目标模型、也是 HALF 的长期移植目标，是**固定密度**重建：它读取
`CHGCAR` 中的平滑价电子密度 `rho~(r)`，不做自洽场（SCF）迭代。在给定 ENCUT
下，波函数用平面波 `|k+G>` 展开，并求解：

$$
H(\mathbf{k})\,\mathbf{c}_n
=\varepsilon_n S(\mathbf{k})\,\mathbf{c}_n .
$$

对于完整的 HAPPY `MIMIC_US` 算符，稠密矩阵为：

$$
H_{\mathbf G\mathbf G'}(\mathbf k)=
\frac{\hbar^2|\mathbf k+\mathbf G|^2}{2m_e}\delta_{\mathbf G\mathbf G'}
+V_{\mathrm{eff}}(\mathbf G-\mathbf G')
+\sum_{Iij}\beta_i^I(\mathbf G)D_{ij}^I\beta_j^{I*}(\mathbf G'),
$$

$$
S_{\mathbf G\mathbf G'}(\mathbf k)=\delta_{\mathbf G\mathbf G'}
+\sum_{Iij}\beta_i^I(\mathbf G)Q_{ij}^I\beta_j^{I*}(\mathbf G').
$$

第一项是动能；`Veff(G-G')` 是固定密度有效局域势的傅里叶系数。最后两项是
PAW/USPP 增强：`beta` 是倒空间 PAW 投影子，`Q` 修正重叠矩阵，`D` 是 onsite
哈密顿量矩阵。

有效势由固定平滑密度推出：

$$
V_{\mathrm{eff}}(\mathbf r)=V_{\mathrm{ion}}^{\mathrm{local}}(\mathbf r)
+V_H[\widetilde\rho](\mathbf r)
+V_{\mathrm{xc}}[\widetilde\rho+\rho_{\mathrm{core}}](\mathbf r),
$$

$$
V_H(\mathbf G)=\frac{4\pi e^2\widetilde\rho_{\mathbf G}}
{\Omega|\mathbf G|^2}\quad(\mathbf G\ne0),
\qquad
D_{ij}^I=D_{ij}^{I,\mathrm{ION}}
+\int V_{\mathrm{eff}}(\mathbf r)Q_{ij}^{I,\mathrm{DEP}}(\mathbf r-\mathbf R_I)\,d^3r.
$$

这里 `rho~_G` 是 HAPPY 的电子数归一化密度系数，`Omega` 是胞体积，`e^2` 是
eV/angstrom 单位制中的静电换算因子。`rho_core` 是 POTCAR 中仅用于 NLCC 的
部分芯密度，不会重复加入 Hartree 密度。`DION` 是 POTCAR 固定 onsite 项，
`QDEP` 才给出势依赖的 MIMIC_US 校正；`VH(G=0)` 是势规范并设为零。

### HALF 当前实际实现的范围

CPU 与 CUDA 路径现已实现上式完整的 Gamma 点 MIMIC_US 项；`--uspp-dij` 用于启用势依赖
校正。GPU 先对周期有效势进行三次 B 样条预滤波，再围绕每个原子在球面网格取样，
投影到实球谐函数，并用 VASP 的双球贝塞尔补偿函数完成径向积分，最后与 AE−PS
多极矩收缩得到每个原子的 $D_{ij}^I$。显式多 k 点能带和对称性约化总能量已经
实现；矩阵自由算符、力和波函数导出仍是计划项。单个非 Gamma k 点可用
`--kpoint KX KY KZ`。权威状态见
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md)。

对已实现的 DION 问题，`S` 正定，可将 `S = L L^H` 作 Cholesky 分解，化为
`L^-1 H L^-H y = epsilon y`，再以 `c = L^-H y` 恢复广义本征矢。CPU 的 MKL 和
GPU 的 cuSOLVER 都执行这一稠密广义厄米求解。

## 代码逻辑与工作流

```text
CHGCAR + POTCAR
  -> 解析晶体结构、平滑密度和多元素 PAW 数据集
  -> 按 ENCUT 选择 Gamma 平面波基
  -> 以 FFT 构造 Veff：Hartree、局域离子势、XC 与 NLCC
  -> 生成倒空间 PAW 投影子，构造 DION/QPAW 矩阵
  -> 可选：GPU 构造 QDEP，得到逐原子的 D = DION + DeltaD
  -> 组装稠密复数 H、S 矩阵并进行厄米化
  -> 广义厄米对角化：CPU 使用 MKL，GPU 使用 cuSOLVER
  -> 输出本征值及 JSON 验证报告
```

CUDA 路径将 FFT 势构造、投影子计算、H/S 组装和本征值求解全部保留在设备端，
避免传输完整的复数 H/S 稠密矩阵。平面波数增加后，三次标度的稠密广义本征
求解会成为主要耗时。主要
源码模块为：

| 模块 | 职责 |
| --- | --- |
| `src/half_chgcar.F90`、`src/half_potcar.F90` | 输入解析 |
| `src/half_basis.F90`、`src/half_fft.F90` | 平面波基与 FFT |
| `src/half_potential.F90`、`src/half_cuda_potential.cuf` | 有效势 |
| `src/half_paw.F90`、`src/half_cuda_assembly.cuf` | PAW 项与 H/S 组装 |
| `src/half_dense_solver.F90`、`src/half_cuda_solver.cuf` | CPU/GPU 广义本征求解 |
| `app/half_cli.F90` | 统一 CLI |

## 构建

需要 Linux x86-64、CMake 3.24+ 与 NVIDIA HPC SDK（`nvfortran`）。CPU 完整
Gamma 求解需要设置 `MKLROOT`；CUDA 路径使用 cuFFT 与 cuSOLVER。

```bash
# CPU
tools/half-cmake cpu all

# CUDA；GPU 架构由工具探测，也可显式设置
tools/half-cmake cuda13 all
HALF_GPU_CC=89 tools/half-cmake cuda13 all
```

也可以使用 CMake preset：

```bash
cmake --preset cpu-release
cmake --build --preset cpu-release

cmake --preset cuda13-release
cmake --build --preset cuda13-release
```

## CLI

```bash
# 检查平滑密度输入
half inspect CHGCAR.smooth --encut 400

# Gamma 固定密度重建；写出 JSON 报告
half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --xc pbe --backend cuda \
  --uspp-dij --reference-eigenval EIGENVAL --output gamma_validation.json

# 任意分数倒空间 k 点
half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --backend cuda --uspp-dij \
  --kpoint 0.125 0.25 0.375 --output kpoint_validation.json

# 显式使用 CPU 后端
half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --xc pbe --backend cpu

# 从显式倒空间 KPOINTS 重建多 k 点能带
half bands CHGCAR.smooth POTCAR KPOINTS \
  --encut 400 --bands 12 --backend cuda --uspp-dij --output bands.json

# 或生成与 ASE 一致的立方晶系默认高对称路径
half bands CHGCAR.smooth POTCAR \
  --encut 400 --bands 12 --npoints 60 --backend cuda --output-prefix bands

# 使用 spglib 不可约 k 网格计算固定密度 Harris 总能量
half energy CHGCAR.smooth POTCAR \
  --encut 400 --kspacing 0.5 --bands 12 --backend cuda \
  --vaspwave-h5 vaspwave.h5 --forces --force-step 0.001 --output energy.json

# 也可在显式网格上读取匹配的 VASP EIGENVAL
half energy CHGCAR.smooth POTCAR \
  --kpoints-file KPOINTS --reference-eigenval EIGENVAL --bands 60 --output energy.json
```

`half-validate-gamma`、`half-bands` 与 `half-energy` 都是对应子命令的兼容别名。

## HfO2 大矩阵速度

条件：HfO2 平滑 CHGCAR、匹配的 PBE `Hf_pv+O` POTCAR、PBE、520 eV、60 bands、
Gamma 点，对应 3407×3407 的 complex128 广义本征问题。CPU 测试将进程固定在
一个逻辑核，所有 OpenMP/BLAS 线程控制均设为 1。

| 实现 | 资源 | wall time |
| --- | --- | ---: |
| HALF CUDA（`--uspp-dij`） | RTX PRO 6000 | 1.696 s 中位数 |
| HALF CPU Fortran（`--uspp-dij`） | 单逻辑核 | 43.790 s |
| HAPPY Python（`--uspp-dij`） | 单逻辑核 | 64.23 s |

数值等价的 CPU Fortran 比 HAPPY 快 **1.47×**；CUDA 比 HAPPY 快 **37.87×**，
比 CPU Fortran 快 **25.82×**。五次 CUDA 新进程运行的
中位数为 `1.696 s`，其中有效势 `0.258 s`、包含 QDEP 的 H/S 组装 `0.349 s`、
cuSOLVER `1.061 s`。更早的 4090 数据是 DION-only 历史记录，不再作为
MIMIC_US 加速比引用。完整原始时间、环境和约束见
[`docs/validation/hfo2_single_core_benchmark.json`](docs/validation/hfo2_single_core_benchmark.json)。

## 验证状态

- Si Gamma/MIMIC_US：与 HAPPY 的前 8 条本征值最大偏差约 `3e-12 eV`。
- HfO2 Gamma/MIMIC_US：60 条本征值最大偏差约 `7.9e-12 eV`。
- band path/k 网格、总能量/力与 `vaspwave.h5` 输出仍在移植中。

详细验证记录见 [`docs/VALIDATION.md`](docs/VALIDATION.md)，功能状态见
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md)。
