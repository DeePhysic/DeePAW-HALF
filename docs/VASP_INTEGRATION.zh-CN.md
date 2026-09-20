# VASP 6.6.0 直接调用 DeePAW-HALF 初始化波函数

本文描述进程内集成方案。该私有分支按仓库所有者授权保存了用户提供的
`vasp-edge-release.6.6.0` 源码；不得从此私有仓库向无 VASP 许可的第三方再分发。
下述插入点已在该源码上核对并测试。

私有分支 `vasp-6.6-half-integration` 已包含完整源码和一键构建入口：

```bash
git switch vasp-6.6-half-integration
tools/half-cmake vasp all
```

CMake 会复制私有源码到构建目录，生成匹配当前 NVHPC/CUDA/cc 的
`makefile.include`，从 oneMKL 源码构建 FFTW3 wrapper，先构建 `libhalf.so`，再
构建并链接 `vasp_std`。源目录不会产生 `.o`、`.mod` 或可执行文件。

## VASP 启动时实际做了什么

传统重启分成两个阶段：

1. 波函数分配前，`main.F` 调用 `INWAV_HEAD`，读取并校验 `NKPTS`、`NBANDS`、
   `ENCUT` 和晶格等头信息。
2. `ALLOCW` 之后，`main.F` 调用 `INWAV_FAST`。它逐自旋、k 点和能带读取全局串行
   平面波系数，必要时按新晶格/截断能重映射，再用 `DIS_PW_BAND` 写入分布式
   `W%CPTWFP`，同时填写 `W%CELTOT` 和 `W%FERTOT`。

HALF 不应导入或伪造 PAW 投影系数。后续 VASP 会统一调用 `PROALL` 和 `ORTHCH`，
用自己的投影子重建 `W%CPROJ`，并按 PAW 广义重叠矩阵正交化输入轨道。

电荷密度初始化和波函数初始化是两条独立路径。`ICHARG=2/3` 会进入 VASP 的重叠
原子密度（SAD）路径。因此只返回 HALF 波函数并不能绕过 SAD。DeepAW 预测密度应先
写成 `CHGCAR`，用 `ICHARG=1` 读取（固定密度则用 `ICHARG=11`），HALF 再基于同一
密度构造哈密顿量并求初始轨道。

## 直接调用接口

求解前，VASP 先通过内存传递已经解析好的宿主请求：

```c
int half_set_request_geometry(
    half_handle handle, const half_request_geometry_v1 *geometry,
    char *error, int error_capacity);
```

adapter 把 `LATT_CUR%A` 打包为逐晶格矢量的行主序，把 `T_INFO%POSION` 打包为
逐原子的分数坐标，直接传递从 1 开始的 `T_INFO%ITYP`，并以致密电荷网格
`GRIDC%NGPTAR` 作为 `grid[3]`。HALF 深拷贝所有字段；目前只保存、不影响文件输入
求解，为下一步直接调用 DeepAW API 准备。`half_get_request_geometry()` 提供可审计
的查询/复制路径。

波函数传递使用 ABI v1 的兼容扩展：

```c
int half_solve_kpoint_mapped(
    half_handle handle, const double kpoint[3], int nbands,
    int64_t host_npw, const int32_t *host_gvec,
    double *eigenvalues_eV, double _Complex *eigenvectors, int64_t ld,
    double *overlap_min, double *overlap_max,
    char *error, int error_capacity);
```

`host_gvec` 按宿主需要的顺序存放 `(Gx,Gy,Gz)`。HALF 会验证双方基组完全相同，
再把每个本征矢直接排列成宿主顺序。基组大小不同、缺失 G 向量或宿主向量重复都会
明确报错，绝不会静默截断。

VASP 侧用 `wavedes1%PL_INDEX`、`%PL_COL` 及本地 `%IGX/%IGY/%IGZ` 构造串行
G 向量顺序；每个全局能带得到一列系数后，调用 VASP 现有的 `DIS_PW_BAND` 分发，
并把本征值写入 `W%CELTOT`。占据数初始化以及 `PROALL/ORTHCH` 仍由 VASP 完成。

## 建议的 VASP 控制流

应新增独立逻辑标签（例如 `LHALF_INIT`），不要伪装成 `ISTART=1`；后者还会打开
与 WAVECAR 重启有关的其他语义。

在 `ALLOCW` 之后加入分支：

```text
if LHALF_INIT:
    拒绝尚未支持的计算模式
    从 DeepAW CHGCAR 和 POTCAR 创建 HALF context
    从内存传递晶格、分数坐标、物种编号和致密网格
    对每个不可约 k 点：
        构造 VASP 串行 G 向量表
        调用 half_solve_kpoint_mapped
        对每个全局能带调用 DIS_PW_BAND
        填写 W%CELTOT
    销毁 HALF context
    标记波函数已经初始化（只跳过 WFINIT）
else if ISTART > 0:
    保留原 WAVECAR 路径
else:
    保留原 WFINIT 路径
```

DeepAW 的平滑 CHGCAR 通常只有实空间密度网格，不包含 VASP 写出的 PAW onsite
occupancy 尾部。VASP 6.6.0 的标准 `READCH` 会把缺少该尾部视为整份 CHGCAR
失败，并回退到 SAD。HALF 模式需要一个受控例外：网格与结构校验成功后保留已经
读入的 DeepAW 平滑密度，同时保留 VASP 在 `READCH` 前通过
`DEPATO/SET_RHO_PAW` 初始化的原子 PAW onsite 矩阵。不能把 onsite 矩阵清零，也
不能在普通 `ICHARG=1` 路径上放宽此检查。

推荐输入为：

```text
LHALF_INIT = .TRUE.
ISTART     = 0
ICHARG     = 1
INIWAV     = 1
```

新标签只能阻止 `WFINIT`，不能阻止 `CHGCAR` 读取，也不能跳过后续的 `PROALL` 和
`ORTHCH`。

如需让 HALF 直接从 DeePAW-eSCN 获取重构密度，在启动环境中设置：

```bash
export HALF_ESCN_URL=http://127.0.0.1:8265
```

adapter 会把现有的晶格、分数坐标、从 1 开始的物种编号和致密网格传给
`half_create_from_escn()`。原子序数由匹配的 POTCAR 自动确定，因此无需修改
VASP-to-HALF descriptor。HALF 只请求密度，将 API 的 C-order 网格转换为内部
顺序，并显式归一到 POTCAR 价电子总数。未设置该变量时继续使用原 CHGCAR 构造器。
该开关目前只改变 HALF 的重构密度；在另行实现 VASP host charge-grid 内存注入前，
VASP 自己的 `ICHARG=1` 密度仍会读取 CHGCAR。

可复制进合法 VASP 源码树的独立 adapter 位于
[`examples/vasp/half_vasp_init.F`](../examples/vasp/half_vasp_init.F)。它只包含
HALF 项目自己的桥接代码，不包含或再分发 VASP 源码。

## MPI 与 GPU 约束

- 首个正确性版本限定为全复数、共线、`ISPIN=1`、`KPAR=1`。通过
  `DIS_PW_BAND` 分发后，`NCORE>1` 可以支持，不应直接写内部局部分块。
- 每个 k 点组只允许一个根进程执行 HALF 稠密求解；所有进程仍共同进入 VASP 的
  分发例程。继续前广播 HALF 状态和本征值。
- 扩展 `KPAR>1` 时，每个 k 点组根进程各建一个 HALF context，只求本组 k 点；
  不能让每个 MPI rank 重复启动同一个 GPU 稠密求解。
- 在明确实现自旋密度、旋量和共轭配对约定前，应拒绝 `ISPIN=2`、
  `LNONCOLLINEAR`、SOC 和 Gamma 压缩存储。
- CUDA VASP 应使用相同主版本的 NVHPC/CUDA 工具链链接 `libhalf`，HALF 必须遵循
  VASP 已为该 rank 选择的设备。当前 ABI 经主机内存返回最终系数，功能正确但尚非
  零拷贝。

## Pro 6000 实测

2026-09-19 在 RTX PRO 6000 Blackwell、NVHPC 26.5、其自带 CUDA 13.2、VASP
6.6.0 full-complex 构建上完成了 HfO2 实测：520 eV、54×54×54 密度网格、3407
个平面波、56 条能带、Gamma 点。

- 1 rank：HALF 初始化和一轮 VASP Davidson 总 wall time 4.90 s，首轮 rms
  `2.81e-4`。
- 2 ranks、`NCORE=2`：4.66 s；能量和首轮 rms 与单 rank 一致，证明
  `DIS_PW_BAND` 分发路径可用。
- 原生 SAD + 随机 `WFINIT` 基线：2.63 s，但首轮 rms 为 158。HALF 的首轮残差
  小约 `5.6e5` 倍；其额外耗时是
  一次性启动成本，不能用这个单 SCF 步测试代表完整收敛总耗时。
- HALF 两次运行均读得 88.0000237 个电子，日志中没有
  `charge density of overlapping atoms calculated`，并完成能量和力计算。

机器可读记录见
[`docs/validation/vasp_hfo2_pro6000.json`](validation/vasp_hfo2_pro6000.json)。
编译时使用 `/opt/intel/oneapi/mkl/2026.0` 的 MKL FFTW3 wrapper；CUDA 编译工具包
实际为 NVHPC 26.5 自带的 13.2，不把仅存在文档的 `/usr/local/cuda-13.4` 误报为
可用工具链。

进一步的完整 SCF 对照采用 `KSPACING=0.35`、36 个不可约 k 点、
`LMAXPAW=-1`、`ALGO=All`、`EDIFF=1E-4`，并让 `PREC=High` 自动选择 500 eV。
HALF 用 4 个电子迭代收敛，SAD 需要 16 个；包含 HALF 36 个 k 点初始化的端到端
wall time 分别为 193.44 s 和 266.96 s。详细参数、逐步收敛轨迹、校验和及限制见
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md`](validation/HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md)。
不设置任何 `LMAX*` 标签时，VASP 报告自动模式 `LMAXPAW=-100`、默认
`LMAXMIX=2`；HALF/SAD 分别用 5/16 步收敛。该消融及跨轮 timing 限制也已收入
同一报告。

## 后续必须通过的验证

1. 消除任意能带相位后，与离线 HALF 结果逐列比较映射后的系数。
2. 在 VASP `ORTHCH` 后检查 `C^H S C = I`。
3. 对同一个 DeepAW `CHGCAR` 比较随机 `WFINIT`、WAVECAR 重启和 HALF 初始化的
   首步残差及最终收敛能量。
4. 已确认 `LHALF_INIT + ICHARG=1` 的输出中没有进入重叠原子平滑密度路径。
5. `NCORE=1` 与 `NCORE=2` 已通过；仍需单独实现并测试 `KPAR>1`。

当前文件构造器仍会重新解析 `CHGCAR/POTCAR`。后续可增加内存构造器，直接接收
VASP 已解析的晶格、密度网格和 PAW 数据，以避免重复 I/O；这属于性能优化，不是
正确性前提。
