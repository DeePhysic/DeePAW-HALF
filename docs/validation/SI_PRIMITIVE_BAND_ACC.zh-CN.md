# Si 原胞 CLI 能带：EVD 与矩阵自由 ACC

## 配置

- 节点：RTX PRO 6000，单 GPU
- 输入：2 原子 Si 原胞、DeepAW 平滑电荷密度、PAW POTCAR
- `ENCUT=520 eV`，12 条带，60 个路径点
- 路径：`GXWKGLUWLK,UX`
- 横坐标在不连续的 `K,U` 边界断开；能量零点为 EVD 的 VBM

![Si 原胞能带及逐带 ACC 误差](assets/si_primitive_bandstructure_evd_acc.png)

## 结果

| 求解模式 | 墙钟时间 | 最低 4 条占据带最大误差 | 最低 8 条带最大误差 | 12 条带最大误差 |
|---|---:|---:|---:|---:|
| 稠密 EVD | 7.37 s | 参考 | 参考 | 参考 |
| ACC，最多 40 次 | 8.90 s | `2.49e-13 eV` | `3.79e-12 eV` | `1.820 eV` |
| ACC，最多 200 次 | 31.69 s | `2.96e-13 eV` | `2.96e-13 eV` | `0.464 eV` |

默认 ACC 与 EVD 的间接带隙分别为 `0.5635976766544735 eV` 和
`0.5635976766545410 eV`。因此 HALF--VASP 使用的 ACC 路径已经正确得到
占据子空间、低能空带和带隙；大误差只出现在当前初始化用途不需要的最高几条
空带。默认 40 次 ACC 比 EVD 慢 21%，原因是这个 Si 原胞每个 k 点只有约
1100 个平面波，稠密问题太小，迭代、FFT 和正交化开销不能摊薄。ACC 的目标是
避免大体系的 `NPL x NPL` 稠密矩阵，并为 VASP 提供足够准确的初始波函数；它
不是这个小算例的计时优势路径。

`half bands` 现已接通 `--acc-tol`、`--acc-max-iter` 和
`--acc-block-size`，可按目标空带范围调节精度。原始结果为
[`EVD`](assets/si_primitive_bands_evd.json)、
[`ACC-40`](assets/si_primitive_bands_acc.json) 和
[`ACC-200`](assets/si_primitive_bands_acc_strict.json)。

## 与 VASP 6.6.0 的逐点对比

另以 VASP `ALGO=All` 自洽得到的 Si CHGCAR 作为共同固定密度输入。VASP
使用 `ICHARG=11`、`LMAXMIX=-1`，HALF 使用稠密 EVD；两边均为 520 eV、
12 条带和完全相同的 60 个显式 k 点，并分别按自己的 VBM 对齐。

![HALF-EVD 与 VASP ALGO=All 能带](assets/si_half_evd_vs_vasp_bands.png)

VASP 与 HALF 的间接带隙分别为 `0.604772 eV` 和 `0.604765817 eV`，差
`-0.00618 meV`。占据带 MAE/最大误差为 `0.00367/0.01672 meV`；最低
8 条带为 `0.00439/0.05765 meV`。第 11--12 条最高空带在 VASP 的能量
收敛停止条件下没有达到同等级的本征残差，因而不用于低能能带精度结论。
机器可读统计见
[`si_half_evd_vs_vasp_bands.json`](assets/si_half_evd_vs_vasp_bands.json)。
