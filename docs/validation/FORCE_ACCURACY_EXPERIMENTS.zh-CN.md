# HALF 解析力精度实验

本报告只记录 `force-accuracy-experiments` 分支上的实验，不改变稳定分支的默认
结果。目标是分别回答：

1. HALF 解析力是否等于同一个 HALF Harris 能量的负梯度；
2. HALF 与 VASP `LMAXMIX=-1` MIMIC_US 力的剩余差异来自数值求积还是 onsite
   PAW 模型。

## 1. 可切换的 QDEP 球面求积

原实现把 QDEP 球面求积固定为 12 个 Gauss-Legendre 极角点和 24 个均匀方位角点。
实验分支增加两个 CMake cache 参数：

```bash
-DHALF_QDEP_NTHETA=24
-DHALF_QDEP_NPHI=48
```

CPU Fortran 与 CUDA Fortran 使用同一参数，因此能量中的 `Dij`、解析力中的
`dDij/dR` 和 projector contraction 不会因后端使用不同角网格而失配。稳定默认值
仍是 `12×24`，便于做严格 A/B 对比。

## 2. 测试设置

- GPU：RTX Pro 6000；测试时存在另一 GPU 进程，因此 wall time 只作竞争环境记录，
  不用于加速结论。
- NVHPC 26.5，CUDA 13.2，complex128/float64。
- PBE，`ISPIN=1`，`SIGMA=0.02 eV`，`KSPACING=0.35`。
- Si：`ENCUT=520 eV`，24 bands，28 个 full-mesh k 点。
- HfO2：`ENCUT=521 eV`，80 bands，使用 VASP `IBZKPT`。
- 中心差分位移：atom 1、x 方向、`h=0.001 Angstrom`。
- CHGCAR 只提供固定平滑密度；augmentation 与原子 PAW 项全部来自 POTCAR。

## 3. 结果

### 3.1 与同一 HALF 能量的导数一致性

| 体系 | QDEP 网格 | 解析力 (eV/Angstrom) | 中心差分 (eV/Angstrom) | 绝对误差 |
|---|---:|---:|---:|---:|
| Si atom 1/x | 12×24 | 0.006100663 | 0.004351093 | 0.001749570 |
| Si atom 1/x | 24×48 | -0.000148095 | -0.000191646 | 0.000043550 |
| HfO2 atom 1/x | 12×24 | -0.104255268 | -0.107329684 | 0.003074416 |
| HfO2 atom 1/x | 24×48 | -0.103968772 | -0.103712466 | 0.000256306 |

Si 的导数误差改善约 40.2 倍，HfO2 改善约 12.0 倍。Si 全结构最大力从
`0.006975` 降到 `0.000727 eV/Angstrom`。这证明 12×24 球面积分不足是此前
augmentation/QDEP 解析力离散不一致的主要来源。

### 3.2 与 VASP 初始 MIMIC_US 状态比较

| 体系 | QDEP 网格 | 力分量 MAE (eV/Angstrom) | 最大分量差 |
|---|---:|---:|---:|
| Si | 12×24 | 0.005320861 | 0.007808088 |
| Si | 24×48 | 0.000692031 | 0.001559912 |
| HfO2 | 12×24 | 0.068045473 | 0.319542268 |
| HfO2 | 24×48 | 0.068019383 | 0.319255772 |

Si 的 VASP 力 MAE 改善约 7.69 倍。HfO2 的 HALF 内部导数已经明显收敛，但
VASP MAE 几乎不变，说明其主误差不是球面求积，而是 HALF 与 VASP 的 onsite
PAW/QDEP/one-centre Hamiltonian 尚未完全同构。

### 3.3 能量变化

| 体系 | 12×24 HALF free energy (eV) | 24×48 (eV) | 变化 (eV) |
|---|---:|---:|---:|
| Si | -43.314024964 | -43.314030732 | -0.000005769 |
| HfO2 | -121.886526811 | -121.885912868 | 0.000613942 |

力比总能量对 QDEP 角向求积更敏感，特别是高角动量、多元素体系。

## 4. 下一步实验

1. 测试 `18×36`、`24×48`、`32×64`，以解析力–中心差分误差和额外 QDEP 时间
   选择最低充分阶数，而不是直接把最大网格设为默认。
2. 对所有对称不等价原子的三个方向使用 `h=0.002, 0.001, 0.0005 Angstrom`
   检查中心差分收敛。
3. 将 QDEP 球面求积时间从总 eigensolver 时间中单独计时；无 GPU contention
   后再做性能比较。
4. 在 HfO2 上逐项比较 onsite density matrix、`Dij`、AE/PS one-centre energy
   和 VASP 对应量。提高角网格不能替代这一步。

原始输出位于：

- `/data/limusen/deepaw_half_force_compare/Si/qdep_24x48`
- `/data/limusen/deepaw_half_force_compare/HfO2/qdep_24x48`

