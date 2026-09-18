# DeePAW–HALF

**HALF**（**H**arris **A**ssociative **L**inearized **A**ugmented **P**lane
**W**ave **F**ortran）是 DeePAW–HAPPY 的 CPU/CUDA Fortran 固定密度电子结构
重建实现。它不是新的机器学习哈密顿量：输入是 DeepAW 或 VASP 产生的平滑
`CHGCAR` 密度以及匹配的 PAW `POTCAR`，输出为固定密度下的本征值与验证报告。

```text
结构 -> DeepAW -> 平滑 CHGCAR -> HALF -> H(k), S(k), 能带/本征值
```

HAPPY 是移植过程中的数值 oracle。当前完整验证通过的是 Si 的 Gamma/DION
固定密度路径；HfO2 的大矩阵性能已经记录，但其与 HAPPY 的逐本征值一致性仍是
独立验证门槛。

## 数值模型：从固定密度到本征值

HALF 是**固定密度**重建：它读取 `CHGCAR` 中的平滑价电子密度 `rho~(r)`，不做
自洽场（SCF）迭代。在给定 ENCUT 下，波函数用平面波 `|k+G>` 展开，并求解：

```text
H(k) c_n = epsilon_n S(k) c_n
```

在原子单位制下，Gamma 路径实际组装的稠密矩阵为：

```text
H_GG' = |k+G|^2/2 * delta_GG' + Veff(G-G')
        + sum_(I,i,j) beta_i^I(G) D_ij^I beta_j^I*(G')

S_GG' = delta_GG' + sum_(I,i,j) beta_i^I(G) Q_ij^I beta_j^I*(G')
```

第一项是动能；`Veff(G-G')` 是固定密度有效局域势的傅里叶系数。最后两项是
PAW/USPP 增强：`beta` 是倒空间 PAW 投影子，`Q` 修正重叠矩阵，`D` 是 onsite
哈密顿量矩阵。

有效势由固定平滑密度推出：

```text
Veff(r) = Vion_local(r) + VH[rho~](r) + Vxc[rho~ + rho_core](r)
VH(G)   = 4*pi*rho~(G)/|G|^2,  G != 0
D_ij^I = DION_ij^I + integral Veff(r) QDEP_ij^I(r) dr
```

`rho_core` 是 POTCAR 中用于 NLCC 的部分芯密度；`DION` 是 POTCAR 固定 onsite
项；`QDEP` 给出势依赖的 MIMIC_US 校正。`VH(G=0)` 由选定的静电规范处理。由于
`S` 正定，可将 `S = L L^H` 作 Cholesky 分解，把问题化为
`L^-1 H L^-H y = epsilon y`，再以 `c = L^-H y` 恢复广义本征矢。CPU 的 MKL 和
GPU 的 cuSOLVER 都执行这一稠密广义厄米求解。

## 代码逻辑与工作流

```text
CHGCAR + POTCAR
  -> 解析晶体结构、平滑密度和多元素 PAW 数据集
  -> 按 ENCUT 选择 Gamma 平面波基
  -> 以 FFT 构造 Veff：Hartree、局域离子势、XC 与 NLCC
  -> 生成倒空间 PAW 投影子，构造 D/Q 矩阵
  -> 组装稠密复数 H、S 矩阵并进行厄米化
  -> 广义厄米对角化：CPU 使用 MKL，GPU 使用 cuSOLVER
  -> 输出本征值及 JSON 验证报告
```

CUDA 路径将 FFT 势构造、投影子计算、H/S 组装和本征值求解全部保留在设备端，
避免传输完整的复数 H/S 稠密矩阵。平面波数增加后，三次标度的稠密广义本征
求解会成为主要耗时。主要源码模块为：

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
  --reference-eigenval EIGENVAL --output gamma_validation.json

# 显式使用 CPU 后端
half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --xc pbe --backend cpu
```

`half-validate-gamma` 是 `half gamma` 的兼容别名。`half-bands` 与
`half-energy` 已保留兼容命令名，但在任意 k 点与总能量数值闭合前会明确报告
该功能尚未完成。

## HfO2 大矩阵速度

条件：HfO2 平滑 CHGCAR、匹配的 PBE `Hf_pv+O` POTCAR、PBE、520 eV、60 bands、
Gamma 点，对应 3407×3407 的 complex128 广义本征问题。CPU 测试将进程固定在
一个逻辑核，所有 OpenMP/BLAS 线程控制均设为 1。

| 实现 | 资源 | wall time |
| --- | --- | ---: |
| HALF CUDA | RTX PRO 6000 | 1.60 s |
| HALF CUDA | RTX 4090 | 1.795 s |
| HALF CPU Fortran | 单逻辑核 | 43.94 s |
| HAPPY Python（`--uspp-dij`） | 单逻辑核 | 64.23 s |

HALF CPU 单核比 HAPPY Python 快 **1.46×**；RTX 4090 与 RTX PRO 6000 相比
单核 HALF CPU 分别快 **24.47×** 与 **27.41×**。同一矩阵规模下，PRO 6000
约比 4090 快 **12%**。完整原始时间、环境和约束见
[`docs/validation/hfo2_single_core_benchmark.json`](docs/validation/hfo2_single_core_benchmark.json)。

## 验证状态

- Si Gamma/DION：与 HAPPY 本征值最大偏差约 `1.6e-11 eV`。
- HfO2：已用于大矩阵性能测试，但尚未完成与 HAPPY 的逐带数值验收。
- 任意 k 点、势依赖 MIMIC_US `D`、总能量/力与 `vaspwave.h5` 输出仍在移植中。

详细验证记录见 [`docs/VALIDATION.md`](docs/VALIDATION.md)，功能状态见
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md)。
