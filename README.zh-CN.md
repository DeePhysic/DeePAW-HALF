# DeePAW–HALF

**HALF**（**H**arris **A**ssociative **L**inearized **A**ugmented **P**lane
**W**ave **F**ortran）是 DeePAW 的固定密度电子结构重建方法，包含 Python、
CPU Fortran 和 GPU/CUDA Fortran 三种实现。它不是新的机器学习哈密顿量：输入是
DeepAW 或 VASP 产生的平滑 `CHGCAR` 密度以及匹配的 PAW `POTCAR`，输出为固定
密度下的本征值与验证报告。

```text
结构 -> DeepAW -> 平滑 CHGCAR -> HALF -> H(k), S(k), 能带/本征值
```

三种实现之间的相互一致性构成内部数值交叉验证。当前 GPU/CUDA Fortran
Gamma/MIMIC_US 固定密度路径已在 Si 与 HfO2 上通过逐本征值验证；CPU Fortran
QDEP（默认采用 VASP 式直接 FFT 网格 SETDIJ）、任意 k 点、多 k 点能带、spglib
不可约 k 网格、占据数、Ewald 和 Harris 固定密度总能量也已闭合。Si 总能量及
各分量在 Python 与 Fortran 实现之间的差小于 `4e-12 eV`。旧的中心差分力只
保留为解析力的数值 oracle，不再作为正式力功能。原生解析力目前已经实现
Ewald、倒空间局域势以及
广义 PAW `D-epsilon Q` Hellmann--Feynman 导数，以及 NLCC/`FORCOR` partial-core
导数。现在还会从 POTCAR 的 AE/PS partial waves 与补偿多极矩重建输出态的
augmentation density，通过 `dD_ij/dR` 计入其显式位移力，并把带 augmentation
的输出密度用于力泛函。DeepAW 平滑输入密度保持冻结；Harris 响应包含
输出—输入密度差产生的 Hartree+XC kernel 场，并按 POTCAR `PSPRHO` 的移动
导数收缩。球形原子 PAW double counting
也会由 POTCAR 的 AE/PS partial waves、原子占据、core density、`DEXC` 与补偿
电荷自动重建，不读取 CHGCAR augmentation 尾部，也不依赖 VASP 输出数值。
独立 CLI 的 `--forces` 全程不移动原子；相对 `LMAXMIX=-1` 的完全自洽 VASP，
Si/HfO2 的能量差为 `-1.360/+9.837 meV/atom`，力分量 MAE 为
`0.414/10.988 meV/Angstrom`。原生
`vaspwave.h5` 输出以及 HDF5 电荷/结构/内嵌 POTCAR 输入已经支持。

PAW 无矩阵算符、全带约束最小化、残差预条件、S 度量正交化、重启式
Rayleigh-Ritz、复杂度以及 VASP 直接内存接入的完整推导见
[Harris 无矩阵全带加速原理与推导](docs/HARRIS_ACC_THEORY.zh-CN.md)。
L 通道、Harris Hartree+XC 响应及 Si/HfO2 力验证见
[Harris 力的 L 通道与密度响应优化](docs/validation/HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.zh-CN.md)。
项目对外采用的完全自洽 VASP 对比见
[HALF 能量和力与完全自洽 VASP 的对比](docs/validation/HALF_VS_VASP_SCF_ENERGY_FORCE.zh-CN.md)。
VASP 接入区分两种模式：one-step approximate SCF 只做一次电子更新并得到近自洽
能量；ACC-SCF 继续迭代到 `EDIFF`，获得完全自洽的密度、能量和力。

## 数值模型：从固定密度到本征值

DeePAW-HALF 采用**固定密度**重建：它读取
`CHGCAR` 中的平滑价电子密度 `rho~(r)`，不做自洽场（SCF）迭代。在给定 ENCUT
下，波函数用平面波 `|k+G>` 展开，并求解：

$$
H(\mathbf{k})\,\mathbf{c}_n
=\varepsilon_n S(\mathbf{k})\,\mathbf{c}_n .
$$

对于完整的 DeePAW-HALF `MIMIC_US` 算符，稠密矩阵为：

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

这里 `rho~_G` 是电子数归一化密度系数，`Omega` 是胞体积，`e^2` 是
eV/angstrom 单位制中的静电换算因子。`rho_core` 是 POTCAR 中仅用于 NLCC 的
部分芯密度，不会重复加入 Hartree 密度。`DION` 是 POTCAR 固定 onsite 项，
`QDEP` 才给出势依赖的 MIMIC_US 校正；默认实现把补偿核直接放在 FFT 网格上，
并从 POTCAR partial waves 在对数径向网格上统一重建 `QPAW(i,j,L)`；`VH(G=0)`
是势规范并设为零。

计算力时，HALF 从各 k 点波函数重建 PAW onsite 占据矩阵与 POTCAR
augmentation density：

$$
P_{ij}^{I}=\sum_{n\mathbf k}w_{\mathbf k}f_{n\mathbf k}
\langle\widetilde\psi_{n\mathbf k}|\beta_i^I\rangle
\langle\beta_j^I|\widetilde\psi_{n\mathbf k}\rangle,
\qquad
\rho_{\mathrm{aug}}(\mathbf r)=\sum_{Iij}P_{ij}^{I}
Q_{ij}^{I}(\mathbf r-\mathbf R_I).
$$

$Q_{ij}^{I}$ 的多极矩来自 POTCAR 的 AE−PS partial waves 与双球贝塞尔补偿
函数。令
$\rho_{\mathrm{out}}=\widetilde\rho_{\mathrm{wave}}+\rho_{\mathrm{aug}}$，
局域势力和显式 augmentation 力分别为

$$
\mathbf F_I^{\mathrm{loc}}=-\int\rho_{\mathrm{out}}(\mathbf r)
\frac{\partial V_{\mathrm{loc}}^I}{\partial\mathbf R_I}\,d^3r,
\qquad
\mathbf F_I^{\mathrm{aug}}=-\sum_{n\mathbf k}w_{\mathbf k}f_{n\mathbf k}
\mathbf c_{n\mathbf k}^{\dagger}
\frac{\partial D^I}{\partial\mathbf R_I}\mathbf c_{n\mathbf k},
$$

其中

$$
\frac{\partial D_{ij}^{I}}{\partial\mathbf R_I}
=\int V_{\mathrm{eff}}(\mathbf r)
\frac{\partial Q_{ij}^{I}(\mathbf r-\mathbf R_I)}
{\partial\mathbf R_I}\,d^3r.
$$

这是 PAW 耦合乘积两侧的不同导数，二者都必须保留；CLI 将后一项单列为
`paw_aug`。

### HALF 当前实际实现的范围

CPU 与 CUDA 路径现已实现上式完整的 Gamma 点 MIMIC_US 项；`--uspp-dij` 用于启用势依赖
校正。GPU 先对周期有效势进行三次 B 样条预滤波，再围绕每个原子在球面网格取样，
投影到实球谐函数，并用 VASP 的双球贝塞尔补偿函数完成径向积分，最后与 AE−PS
多极矩收缩得到每个原子的 $D_{ij}^I$。显式多 k 点能带和对称性约化总能量已经
实现；H/S 组装、矩阵应用和波函数导出也已实现。解析 `--forces` 会重建 POTCAR
augmentation occupancy 与 augmentation density；有限差分只保留为验证 oracle。
单个非 Gamma k 点可用
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
  -> 参考路径：组装稠密复数 H、S，使用 MKL/cuSOLVER 对角化
  -> ACC 路径：分块 H*Psi/S*Psi，全带 Rayleigh-Ritz，只求解小型子空间
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
| `src/half_potcar_interp.F90` | 与 POTCAR 表格约定一致的局域势、`PSPCOR`、`PSPRHO` 插值 |
| `src/half_potential.F90`、`src/half_cuda_potential.cuf` | 有效势 |
| `src/half_paw.F90`、`src/half_cuda_assembly.cuf` | PAW 项与 H/S 组装 |
| `src/half_dense_solver.F90`、`src/half_cuda_solver.cuf` | CPU/GPU 广义本征求解 |
| `src/half_cuda_iterative_solver.cuf` | CUDA 无矩阵全带迭代与提前停止 |
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

# CPU MPI（显式启用；普通串行构建没有 MPI 启动开销）
cmake -S . -B build/cpu-mpi -DHALF_ENABLE_CUDA=OFF \
  -DHALF_ENABLE_MPI=ON -DHALF_MKL_ROOT="$MKLROOT"
cmake --build build/cpu-mpi -j

# Intel oneAPI 2023：直接使用 MPI wrapper；该版本推荐 mpiifort
module load intel/oneapi2023
cmake -S . -B build/cpu-oneapi-mpi -DCMAKE_Fortran_COMPILER=mpiifort \
  -DHALF_ENABLE_CUDA=OFF -DHALF_ENABLE_MPI=ON \
  -DHALF_MKL_ROOT="$MKLROOT"
cmake --build build/cpu-oneapi-mpi -j

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

# 无矩阵 Harris 初始波函数：分块 H*Psi/S*Psi、全带 Ritz 优化和残差提前停止
half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 64 --backend cuda --solver acc --uspp-dij \
  --acc-tol 1e-3 --acc-max-iter 40 --acc-block-size 16 \
  --output gamma_acc.json

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
  --vaspwave-h5 vaspwave.h5 --output-prefix energy

# 独立 HALF 解析 PAW/Harris 力；不移动原子
half energy CHGCAR.smooth POTCAR \
  --encut 400 --kspacing 0.5 --bands 12 --backend cuda \
  --forces --output-prefix forces

# 也可在显式网格上读取匹配的 VASP EIGENVAL
half energy CHGCAR.smooth POTCAR \
  --kpoints-file KPOINTS --reference-eigenval EIGENVAL --bands 60 --output energy.json

# 4 个单线程 CPU rank 按 k 点并行
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun --bind-to core -np 4 build/cpu-mpi/half bands \
  CHGCAR.smooth POTCAR KPOINTS --backend cpu --bands 60 --output bands.json
```

`--solver acc` 是大基组的主加速路径。它不构造 `NPL x NPL` 稠密 H/S：局域势
FFT、动能、PAW 投影收缩、残差、预条件、S 正交化及块旋转均留在 GPU；cuSOLVER
只求解 `NBANDS` 或 `2*NBANDS` 的 Rayleigh-Ritz 子空间问题。JSON 会记录迭代次数
和最终最大残差。`--solver evx` 可作为稠密部分谱参考。VASP 中设置
`HALF_MODE=ACC` 即选用相同路径并直接使用 VASP 的 `NBANDS`；默认残差阈值为
`1e-4 eV`。

在 20 原子 CsPbBr3 验证体系（`NPL=20640--20780`、12 个 k 点、105 bands）中，
VASP `HALF_MODE=ACC` 与稠密 `VASP_LIKE` 都使用 5 个 SCF LOOP，最终能量相同；
总耗时从 1861.0 秒降至 275.5 秒，即加速 6.75 倍。记录见
[RTX PRO 6000 验证数据](docs/validation/cspbbr3_vasp_acc_pro6000.json)。

MPI 作用于 CPU 的 `bands` 和 `energy` 路径。rank `r` 负责
`r+1, r+1+nranks, ...` 这些 k 点，集合通信恢复原顺序的本征值与平面波计数，
只有 rank 0 写 JSON。可选的有限差分验证 oracle 对每个位移结构复用同一分发；
生产 `--forces` 绝不移动原子。通常应设置
`OMP_NUM_THREADS=1`、`MKL_NUM_THREADS=1`，避免每个 rank 再生成一组线程。
多 rank CUDA 与多 rank `vaspwave.h5` 会被明确拒绝。

HfO2 的 36-k 点 CPU 实测中，1/2/4 ranks 分别为
116.59/58.44/29.87 秒，即 1.995×/3.903× 加速；2/4 ranks 的 JSON 本征值与
串行逐位相同。记录见
[`docs/validation/hfo2_cpu_mpi.json`](docs/validation/hfo2_cpu_mpi.json)。

`half-validate-gamma`、`half-bands` 与 `half-energy` 都是对应子命令的兼容别名。

`bands --output-prefix NAME` 同时生成 `NAME.json`、`NAME.csv`、`NAME.npz` 与
无外部绘图库依赖的 `NAME.png`，其中包含 VBM、CBM、采样带隙、Gamma 直接带隙、
路径距离、基组规模和重叠矩阵诊断。`energy --output-prefix NAME` 生成 JSON 与
兼容旧格式的 NPZ，后者包含 k 点、权重、本征值、占据数和力。

## 库 API 与 VASP 接入

同一实现可安装为 `libhalf.so`。稳定 C ABI 使用不透明句柄和类似 libxc 的运行时
整数选择器；`half_api.f90` 为 Fortran 调用方提供不依赖编译器 `.mod` 的
`bind(C)` 接口。API 覆盖基组查询、稠密 H/S 组装、成组 H/S 应用和 k 点直接
求解；ABI v1 的 CPU context 支持三层计算，CUDA context 支持全设备驻留的直接
求解。创建 context 时即可选择 CPU/CUDA 与 EVD/EVJ。安装后的 CMake 工程链接
`HALF::half`，也可使用 `half.pc`。

VASP adapter 还会把晶格、分数坐标、物种编号和致密电荷网格尺寸直接传入 HALF
内存。HALF 深拷贝并保存这份请求元数据；当前求解器暂不使用它，从而为下一步直接
调用 DeepAW API 留出稳定接口。

完整生命周期、数据布局和 VASP adapter 方案见
[中文 API 指南](docs/API.zh-CN.md)与
[VASP 6.6.0 接入实测](docs/VASP_INTEGRATION.zh-CN.md)。独立调用示例见
[C 示例](examples/api/half_c_example.c)，VASP 侧桥接模块见
[`examples/vasp/half_vasp_init.F`](examples/vasp/half_vasp_init.F)。
从 DeepAW 密度直接重建能带、计算能量和力、HALF 介入 VASP 加速 SCF 的完整文章，
以及与已报道机器学习势的对比、HfO₂ 逐点能带验证、MP-85 平均迭代加速和大基组
ACC 加速结果见
[验证与性能报告](docs/validation/HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md)。

在私有分支 `deepaw-half-acc` 中，VASP 6.6.0 完整源码位于
`vendor/vasp-6.6.0`，可在 Pro 6000 节点用一条命令完成 HALF、MKL FFTW wrapper
和 `vasp_std` 的配置、编译及 HALF 测试：

```bash
tools/half-cmake vasp all
```

脚本自动探测 `/opt/nvidia/hpc_sdk/Linux_x86_64` 下最新 NVHPC、其自带 CUDA
版本、GPU compute capability、HPC-X MPI，以及 `/opt/intel/oneapi` 的 MKL。
最终程序位于：

```text
build/vasp-cuda<CUDA>-cc<CC>-release/vasp-6.6.0/bin/vasp_std
```

需要让 VASP 本身也使用其 OpenACC GPU port 时可运行：

```bash
tools/half-cmake vasp all -- -DHALF_VASP_ENABLE_OPENACC=ON
```

默认保持已在 Pro 6000 验证的配置：VASP 主程序使用 NVHPC/HPC-X CPU 路径，
HALF 的势、H/S 组装和广义本征求解使用 CUDA GPU。

## HfO2 大矩阵速度

条件：HfO2 平滑 CHGCAR、匹配的 PBE `Hf_pv+O` POTCAR、PBE、520 eV、60 bands、
Gamma 点，对应 3407×3407 的 complex128 广义本征问题。CPU 测试将进程固定在
一个逻辑核，所有 OpenMP/BLAS 线程控制均设为 1。

| 实现 | 资源 | wall time |
| --- | --- | ---: |
| DeePAW-HALF GPU/CUDA Fortran（`--uspp-dij`） | RTX PRO 6000 | 1.696 s 中位数 |
| DeePAW-HALF CPU Fortran（`--uspp-dij`） | 单逻辑核 | 43.790 s |
| DeePAW-HALF Python（`--uspp-dij`） | 单逻辑核 | 64.23 s |

CPU Fortran 实现比 Python 实现快 **1.47×**；GPU/CUDA Fortran 比 Python
实现快 **37.87×**，比 CPU Fortran 快 **25.82×**。五次 CUDA 新进程运行的
中位数为 `1.696 s`，其中有效势 `0.258 s`、包含 QDEP 的 H/S 组装 `0.349 s`、
cuSOLVER `1.061 s`。更早的 4090 数据是 DION-only 历史记录，不再作为
MIMIC_US 加速比引用。完整原始时间、环境和约束见
[`docs/validation/hfo2_single_core_benchmark.json`](docs/validation/hfo2_single_core_benchmark.json)。

## 验证状态

- Si Gamma/MIMIC_US：Python 与 GPU/CUDA Fortran 实现的前 8 条本征值最大
  偏差约 `3e-12 eV`。
- HfO2 Gamma/MIMIC_US：60 条本征值最大偏差约 `7.9e-12 eV`。
- band path/k 网格、总能量/力与 `vaspwave.h5` 已有逐项验证记录。

详细验证记录见 [`docs/VALIDATION.md`](docs/VALIDATION.md)，功能状态见
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md)。
