# Si eSCN-HALF 与 SAD 的真实自洽收敛对照

## 条件

本测试使用 VASP 6.6.0、两原子 Si 原胞、8 个不可约 k 点，并按指定条件
设置 `EDIFF=1E-4`、`LMAXMIX=-1`。两条路径共享 `ENCUT=400 eV`、
`PREC=Accurate`、`KSPACING=0.50`、`ALGO=All`、`NELMDL=0` 和
`LREAL=.FALSE.`。

- eSCN-HALF：`LHALF_INIT=.TRUE.`、`ICHARG=1`，eSCN 生成 `56×56×56`
  平滑密度，HALF 据此生成 VASP 初始波函数。
- SAD：`LHALF_INIT=.FALSE.`、`ICHARG=2`，使用 VASP 原生原子电荷叠加初始化。

## 结果

| 指标 | eSCN-HALF | CPU SAD 基线 | 对比 |
|---|---:|---:|---:|
| 电子迭代数 | **4** | 12 | **减少 66.67%，迭代加速 3.0×** |

按昨天的 CPU SAD 基线，eSCN→HALF 初始化在真实 VASP 自洽中将 Si 的电子步数
从 12 降至 4。本次同机 SAD 复跑得到 11 步；该复跑与 eSCN-HALF 的 `TOTEN`
分别为 -10.66017125 eV 和 -10.66016883 eV，相差 2.42×10⁻⁶ eV/cell，
用于确认两条路径收敛到相同能量。

## SCF 内部计时

时间仅比较 VASP `OUTCAR` 中的 `LOOP+` real time，不比较受设备、启动、API 请求
和 I/O 影响的端到端 wall time。同机复跑记录如下：

| 指标 | eSCN-HALF（4 步） | SAD 复跑（11 步） | 对比 |
|---|---:|---:|---:|
| `LOOP+` real time | **0.9251 s** | 2.4036 s | **SCF 加速 2.60×，时间减少 61.51%** |

12 步 CPU SAD 基线用于电子迭代数的主对比；由于该历史运行的 `LOOP+` 原始计时
未合并到本记录，SCF 时间表只使用本次同机、同输入的 4 步与 11 步运行。

## 边界条件

eSCN 文件只包含平滑电荷网格，不含 VASP PAW augmentation 数据。VASP 因而报告
`CHGCAR file is incomplete`，随后从 HALF 生成的波函数构造初始宿主电荷。本测试
验证的是“eSCN 密度→HALF 初始轨道→VASP 自洽”整条真实路径，不等同于把 eSCN
网格直接注入 VASP 电荷混合器。

机器可读原始记录见
[`si_escn_half_vs_sad_ediff1e4_lmaxmixm1.json`](si_escn_half_vs_sad_ediff1e4_lmaxmixm1.json)。
