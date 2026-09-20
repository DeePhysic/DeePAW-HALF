# HfO2 VASP-like NBANDS 验证

日期：2026-09-20

`deepaw-half-vasp` 分支新增 `HALF_MODE=VASP_LIKE`（本次验证使用其兼容简写
`VASP`）。VASP 将 `WDES%NB_TOT` 作为
`nbands` 通过已有 C ABI 传给 HALF；HALF 随后使用 CUDA 广义 Hermitian 区间
求解器 `cusolverDnZhegvdx`，只计算最低的指定数量本征对。默认值
`HALF_MODE=TRADITIONAL` 执行传统全谱求解（兼容简写为 `DENSE`）。

测试在 NVIDIA RTX PRO 6000 Blackwell 上完成，编译环境为 NVHPC 26.5 及其
CUDA 13.2。两次运行使用完全相同的 HfO2 输入：

- 36 个 k 点、52 条能带、每个 k 点约 3700 个平面波
- `ENCUT=521`、`KSPACING=0.35`、`EDIFF=1E-4`、`ALGO=All`
- `ISPIN=1`、`LMAXMIX=-1`、`ICHARG=1`
- 使用 DeepAW CHGCAR；`LHALF_INIT=.TRUE.`、`LHALF_API=.FALSE.`
- 单 MPI rank、CPU 单线程；使用同一张 GPU 和同一个可执行文件

| HALF 模式 | VASP 总耗时 | SCF LOOP 数 | 最终 EIGENVAL |
|---|---:|---:|---:|
| `VASP` | 105.990 s | 6 | 参考值 |
| `DENSE` | 135.934 s | 6 | 最大绝对差 0.0 eV |

这组端到端测试中，VASP-like 模式加速 **1.283 倍**，总耗时降低 **22.0%**；
SCF 轨迹和 LOOP 数不变。36 个 k 点的输出均明确记录
`requested NBANDS=52 (VASP-like partial spectrum)`。

运行文件保存在 Pro 6000 的
`/data/limusen/half-vasp-like-validation/{vasp,dense}`，不纳入仓库。
