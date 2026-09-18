# DeePAW–HALF

**HALF** 是 DeePAW–HAPPY 固定密度电子结构重构器的 CPU/CUDA Fortran 实现。项目名
沿用 **H**arris **A**ssociative **L**inearized Augmented Plane Wave **F**ortran。

科学计算契约与 HAPPY 保持一致：

```text
周期结构 -> DeepAW -> 平滑 CHGCAR -> HALF -> H(k), S(k), 能带、波函数、能量
```

HALF 不是新的神经网络 Hamiltonian。它从预测的平滑密度和严格匹配的 VASP PAW
数据重构局域势与线性化 PAW/USPP-like 算符。移植期间，HAPPY 是逐项数值基准。

## 当前里程碑

0.3 版本已经闭合 Si Gamma 固定密度 PBE 路径：

- 读取 CHGCAR 的结构头和平滑密度，并在恰好 `NGX*NGY*NGZ` 个数值后停止；
- 晶格、倒格矢和平面波基数据结构；
- 与 HAPPY 相同截断约定的 Gamma 点平面波筛选；
- 平滑电子数检查；
- 多数据集文本 POTCAR 读取；
- 完整 FFT 网格上的 Hartree、离子势、NLCC、LDA/PBE；
- PAW 倒空间投影子和 DION/QPAW 重叠矩阵；
- MKL CPU 与 cuSOLVER GPU 广义本征求解；
- 完全独立的 CPU Fortran 与 CUDA Fortran 后端；
- CUDA 12.4 与 CUDA 13.0 构建配置；
- 带 checksum 校验的 CPU/GPU 性能基准；
- `half-inspect` 命令和面向数值对齐的测试。

Si Gamma/DION 路径与 HAPPY 的本征值误差约为 `1.6e-11 eV`。任意 k 点、
势依赖 MIMIC_US D、能量/力和 `vaspwave.h5` 的进度见
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md)，因此当前版本还不是 HAPPY
所有工作流的完整替代品。

## 构建与远端验证

推荐使用仓库内统一的 CMake 工具；它会配置独立构建目录、编译并运行 CTest：

```bash
tools/half-cmake cpu all
tools/half-cmake cuda12 all
tools/half-cmake cuda13 all
HALF_GPU_CC=89 tools/half-cmake cuda12 all
```

设置 `MKLROOT`（或传入 `-DHALF_MKL_ROOT=...`）可启用完整网格 FFT 和稠密
Gamma 求解器。CUDA 后端要求 NVIDIA HPC SDK 的 `nvfortran`。

对应的原始 CMake 命令为：

```bash
# 纯 CPU Fortran，不编译 .cuf，也不链接 CUDA runtime
cmake --preset cpu-release
cmake --build --preset cpu-release

# CUDA 12.4：NVHPC 24.5 或更新版本
cmake --preset cuda12-release
cmake --build --preset cuda12-release

# CUDA 13.0：需要 NVHPC 25.9 或更新版本
cmake --preset cuda13-release
cmake --build --preset cuda13-release
```

本项目配置的 RTX 3080 环境为 `hz.icqms.group:8122`：

```bash
scripts/remote_build.sh hz.icqms.group 8122
scripts/remote_benchmark.sh hz.icqms.group 8122 5000
scripts/remote_cuda13_build.sh hz.icqms.group 8122 5000
```

可选的第四个参数用于覆盖 GPU 计算能力；默认通过 `nvidia-smi` 自动检测（例如
RTX 3080 为 `86`，RTX 4090 为 `89`）。不同架构使用独立构建目录，避免共享用户
目录时互相覆盖。

脚本仅同步 HALF 源码和不含许可数据的 Si 平滑 CHGCAR，不会复制 POTCAR。
CUDA 12 使用集群已有的 NVHPC 24.5 构建。CUDA 13 脚本完全不使用容器：它调用
用户已有的 mamba，在共享用户安装下创建持久的 `half-cuda13` 环境，验证 NVIDIA
官方 RHEL/Rocky NVHPC 25.9 RPM，并将 NVHPC 安装到
`/share/home/limusen/app`，不在 `/tmp` 中安装 mamba 环境。

## 验收原则

每个生产模块必须先在完全相同输入上通过 HAPPY 对照，再做性能优化。Si 是第一
道验收门，多元素 HfO2 是声明多元素支持前的必要门。误差阈值见
[`docs/VALIDATION.md`](docs/VALIDATION.md)。
