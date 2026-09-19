# HfO₂：VASP 6.6.0 中 HALF 与 SAD 初始化的收敛对比

## 摘要

本报告比较 VASP 6.6.0 在同一 12 原子 HfO₂ 体系上的两条启动路径：

1. DeepAW 平滑密度（`ICHARG=1`）加 DeePAW-HALF 初始波函数；
2. VASP 原生叠加原子密度 SAD（`ICHARG=2`）加默认初始波函数。

最终测试采用 `PREC=High`、`KSPACING=0.35`、`KGAMMA=.TRUE.`、
`LMAXPAW=-1`、`ALGO=All`、`EDIFF=1E-4` 和 `NELMDL=0`。目录中不存在
`KPOINTS` 文件，因此 VASP 确实从 `KSPACING` 生成网格。VASP 得到 4×4×4
Gamma-centered 网格，经对称性约化后为 36 个不可约 k 点。

HALF 路径用 4 个电子迭代达到收敛，SAD 路径需要 16 个；电子迭代减少 75%。
包含 HALF 对全部 36 个 k 点生成初始波函数在内，端到端 wall time 从
266.96 s 降到 193.44 s，即 1.380× 加速、wall time 减少 27.54%。两条路径
最终 `TOTEN` 相差 `3.813×10⁻⁶ eV/cell`。

## 重要参数更正

本报告只采用下面的最终参数；此前探索性单 Gamma 测试不计入结论。

- 正确标签是 `LMAXPAW=-1`，不是 `LMAXMIX=-1`。前者是 HAPPY HfO₂
  MIMIC_US 参考所用的 VASP onsite-density 设置。
- 显式 `KPOINTS` 会覆盖 `KSPACING`。最终测试已删除 `KPOINTS`，由
  `KSPACING=0.35` 和 `KGAMMA=.TRUE.` 生成网格。
- 最终测试没有写 `ENCUT`。对本 POTCAR，最大 `ENMAX` 为 400 eV；VASP 6.6.0
  在 `PREC=High` 下把它提高 25%，OUTCAR 中实际值为 **500 eV**。这不是原
  HAPPY 参考的 520 eV；如果要求严格复现 520 eV 参考，必须显式写
  `ENCUT=520`。

## 软件与硬件

| 项目 | 配置 |
|---|---|
| 主机 | `hpc.icqms.group:2222`，BH000，Rocky Linux |
| GPU | NVIDIA RTX PRO 6000 Blackwell Workstation Edition，ECC 关闭 |
| 驱动 / 运行时 | 595.84 / CUDA 13.2 |
| CPU | 2× Intel Xeon Platinum 8488C，96 物理核、192 逻辑 CPU |
| 编译器 | NVHPC `nvfortran 26.5` |
| VASP | 6.6.0，full-complex `vasp_std`，单 MPI rank，`KPAR=1` |
| HALF | 分支 `vasp-6.6-half-integration`，提交 `7c39293` |
| HALF 后端 | CUDA，全 GPU MIMIC_US/QDEP 初始求解 |
| VASP 后续求解 | VASP 原生 PAW，`ALGO=All`（OUTCAR：`IALGO=58`） |

VASP 主程序构建时未启用 OpenACC；HALF 的势、MIMIC_US/QDEP、H/S 组装和
广义本征求解在 GPU 上执行。报告中的 wall time 是完整 VASP 进程的单次实测，
未做多次采样统计。

## 体系与输入完整性

| 项目 | 数值 |
|---|---:|
| 原子数 | 12（4 Hf + 8 O） |
| 价电子数 | 88 |
| 能带数 | 60 |
| VASP FFT 网格 | 40×40×40 |
| 致密电荷网格 | 54×54×54 |
| 自动 ENCUT | 500 eV |
| 全 k 网格 | 4×4×4 Gamma-centered |
| 不可约 k 点 | 36 |
| 每个 k 点平面波数 | 3184–3243 |

关键输入校验和：

```text
POSCAR  9824f8a27d3908071c10427a1e376d90cb083ecd3243c84a87213bda7a21276d
POTCAR  cc2e3ef4fe34f6ce782810bea297345f8bedd1ea5a55dc81b775427355a3c41e
CHGCAR  4d66ae56633352553c8d59ae17797e22e8b325e0ee266b7a5ff109b664fe888b
```

POTCAR 顺序为 `PAW_PBE Hf_pv 06Sep2000` 和 `PAW_PBE O 08Apr2002`。

## 最终 INCAR

HALF 路径：

```text
SYSTEM = HfO2 HALF KSPACING 0.35 ALGO All
LHALF_INIT = .TRUE.
ISTART = 0
ICHARG = 1
ISPIN = 1
ISYM = 2
PREC = High
KSPACING = 0.35
KGAMMA = .TRUE.
NBANDS = 60
NELM = 80
NELMDL = 0
EDIFF = 1E-4
LMAXPAW = -1
ALGO = All
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
NWRITE = 2
```

SAD 路径只改变：

```text
SYSTEM = HfO2 native SAD KSPACING 0.35 ALGO All
LHALF_INIT = .FALSE.
ICHARG = 2
```

## 运行时确认

HALF 日志确认 VASP 已把结构和致密网格直接传入 HALF context：

```text
HALF retained VASP host system: 12 ions, 2 types, grid 54 54 54
```

OUTCAR 对两条路径均给出：

```text
Grid dimensions derived from KSPACING:
k-points NKPTS = 36
ENCUT = 500.0 eV
LMAXPAW = -1
IALGO = 58
```

HALF 为全部 36 个不可约 k 点返回初始波函数。不同 k 点包含 3184–3243 个
平面波，重叠矩阵最小本征值约为 0.5274，最大本征值约为 2.126–2.137。

## 收敛结果

| 指标 | HALF | SAD | 比较 |
|---|---:|---:|---:|
| 电子迭代数 | **4** | 16 | 减少 75% |
| Wall time | **193.44 s** | 266.96 s | 1.380× 加速 |
| Wall time 减少 | — | — | 27.54% |
| 最终 TOTEN | -121.10188267 eV | -121.10187886 eV | 差 3.813×10⁻⁶ eV |
| 最后一步 `dE` | -3.768×10⁻⁵ eV | -4.693×10⁻⁵ eV | 均满足 `EDIFF=1E-4` |

### HALF 电子迭代

```text
SDA: 1  E=-121.097377345 eV  rms=2.26E-02
CGA: 2  E=-121.101599216 eV  dE=-4.2219E-03  rms=7.21E-04
CGA: 3  E=-121.101844991 eV  dE=-2.4577E-04  rms=9.99E-05
CGA: 4  E=-121.101882671 eV  dE=-3.7680E-05  rms=2.27E-05
```

### SAD 电子迭代

SAD 首步残差为 `1.49×10²`。前 5 步是默认波函数的 DAV 松弛；随后
`ALGO=All` 的 SDA/CGA 从第 6 步继续，直到第 16 步：

```text
DAV:  1  E= 535.796475637 eV  rms=1.49E+02
DAV:  5  E=-143.425843006 eV  rms=9.74E-02
SDA:  6  E= -83.254783716 eV  rms=1.71E+02
CGA: 15  E=-121.101831932 eV  dE=-1.9284E-04  rms=1.22E-04
CGA: 16  E=-121.101878858 eV  dE=-4.6926E-05  rms=2.87E-05
```

## 默认 LMAX 消融测试

为了区分 `LMAXPAW=-1` 对结果的影响，另做一组完全不设置 `LMAXPAW` 和
`LMAXMIX` 的测试。其余参数与主测试相同：`PREC=High`、`KSPACING=0.35`、
`KGAMMA=.TRUE.`、`ALGO=All`、`EDIFF=1E-4`、`NELMDL=0`，仍为 36 个不可约
k 点和自动 500 eV 截断能。

VASP 6.6.0 在 OUTCAR 中报告的实际默认值是：

```text
LMAXPAW = -100   # VASP 自动模式的内部哨兵值
LMAXMIX = 2
```

默认 LMAX 下的同轮 HALF/SAD 对比为：

| 指标 | HALF | SAD | 比较 |
|---|---:|---:|---:|
| 电子迭代数 | **5** | 16 | 减少 68.75% |
| Wall time | **145.76 s** | 262.70 s | 1.802× 加速 |
| Wall time 减少 | — | — | 44.51% |
| 最终 TOTEN | -121.040359303 eV | -121.040376497 eV | 差 1.7194×10⁻⁵ eV |
| 最后一步 `dE` | -9.3122×10⁻⁵ eV | -4.7930×10⁻⁵ eV | 均满足 `EDIFF=1E-4` |

不同 onsite 设置的迭代数汇总：

| VASP onsite 设置 | HALF 电子迭代 | SAD 电子迭代 |
|---|---:|---:|
| `LMAXPAW=-1` | **4** | 16 |
| 不设置 `LMAX*`（默认自动模式） | **5** | 16 |

因此 `LMAXPAW=-1` 让 HALF 再少一个电子迭代；SAD 在两种设置下均为 16 步。
默认设置下 HALF 仍将迭代数降低 68.75%。

`LMAXPAW=-1` 与默认自动模式之间的最终绝对能量相差约 `0.0615 eV/cell`，这是
PAW onsite-density 处理不同导致的能量规范变化。只能在相同 LMAX 设置内比较
HALF 与 SAD，不能把两种设置的绝对 `TOTEN` 混为同一回归基准。

默认 LMAX 运行目录：

```text
HALF: /data/limusen/half-vasp-kspacing035-algoall-half-default-lmax
SAD:  /data/limusen/half-vasp-kspacing035-algoall-sad-default-lmax
```

默认 LMAX 的 HALF 虽然比 `LMAXPAW=-1` 多一个电子步，但该次 wall time 反而更短。
这说明两轮之间存在节点负载、缓存或频率状态波动；不能用跨轮单样本 wall time
推断默认 LMAX 更快。LMAX 消融应优先比较电子迭代数，同一轮 HALF/SAD wall time
只作为端到端观测值。

## 解释

HALF 的优势主要表现为更接近自洽子空间的初始波函数。对于同样的 36 个 k 点和
60 条能带，SAD 路径需要先消除默认初始波函数的大残差，然后才能进入稳定的
`ALGO=All` 收敛区间；HALF 从第一步就直接处于这一收敛区间。

最终能量在微电子伏特每晶胞量级一致，说明两条路径收敛到同一个 VASP 原生 PAW
解。HALF 只负责初始波函数；波函数交给 VASP 后，投影系数、PAW overlap、onsite
项以及后续 SCF 都由 VASP 原生 PAW 路径处理。

## 公平性与限制

- 这是两条完整启动工作流的比较：DeepAW 密度 + HALF 波函数，对比 SAD 密度 +
  默认波函数。它衡量“绕过 SAD”的实际收益，不是只替换波函数、保持初始密度
  完全相同的消融实验。
- HALF wall time包含 36 次 GPU 稠密广义本征求解及系数传回；没有从总时间中扣除
  初始化成本。
- 每条路径只运行一次；wall time 适合判断本次端到端差异，不是严格的统计性能
  分布。尤其不能直接比较两轮不同 LMAX 测试的 wall time。电子迭代数和相同
  LMAX 下的最终能量是更稳定的比较指标。
- 本测试自动使用 500 eV。它不能直接替代 HAPPY 文档中显式 520 eV 的绝对能量
  回归结果。
- 当前限制仍为共线 `ISPIN=1`、full-complex、`KPAR=1`；尚未覆盖 SOC、非共线和
  `KPAR>1`。

## 可复现位置

```text
HALF: /data/limusen/half-vasp-kspacing035-algoall-half
SAD:  /data/limusen/half-vasp-kspacing035-algoall-sad
VASP: /data/limusen/half-vasp-oneclick/build/vasp-cuda13.2-cc120-release/vasp-6.6.0/bin/vasp_std
```

机器可读摘要见
[`hfo2_vasp_half_vs_sad_kspacing035.json`](hfo2_vasp_half_vs_sad_kspacing035.json)。
