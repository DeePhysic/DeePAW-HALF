# HfO2 的 VASP/HALF `NWAVE` 边界不一致

## 摘要

MP-85 批量测试的 `mp-352 / HfO2` 在 `ENCUT=520 eV`、第 34 个 k 点触发：

```text
HALF: host and HALF plane-wave bases have different sizes
DeePAW-HALF k-point solve failed on the root rank
```

VASP 和 HALF 使用相同的截断公式 `E(G+k) = HSQDTM |G+k|^2 < ENCUT`，
`HSQDTM`、`2*pi` 和严格小于号也一致。差异来自用于计算 `NWAVE` 的晶格精度。

## 复现与证据

- 原始目录：`/data/limusen/deepawchgs/054_mp-352_HfO2`
- `KSPACING=0.35`，`KGAMMA=.TRUE.`，`EDIFF=1E-4`
- `ISPIN=1`，`LMAXMIX=-1`，`ENCUT=520 eV`
- 失败点：k 点 34，`(0.0, 0.5, 0.5)`
- VASP：3672 个平面波
- HALF：3676 个平面波

VASP 从 POSCAR 使用约 16 位精度晶格：

```text
 5.0279225900000002  0.0000000000000000 -0.8189099300000000
 0.0000000000000000  5.1506447000000000  0.0000000000000000
-0.0283774200000000  0.0000000000000000  5.2761341799999997
```

DeepAW CHGCAR 头部只保留六位小数：

```text
 5.027923  0.000000 -0.818910
 0.000000  5.150645  0.000000
-0.028377  0.000000  5.276134
```

用两套晶格和同一公式重新枚举可精确复现 3672/3676。六位小数晶格额外包含
`(-1,9,-1)`、`(-1,-10,-1)`、`(1,9,0)`、`(1,-10,0)` 四个 G 向量；其计算动能
均为 `519.9999852228437 eV`，刚好落入 520 eV 截断球。

## 代码原因

VASP 在 `wave.F` 的 `GEN_LAYOUT/GEN_INDEX` 路径中用 `LATT_CUR%B` 构造基组；
HALF 在 `half_basis.F90::build_plane_wave_basis` 中使用 `crystal%reciprocal`。

adapter 已通过 `half_set_request_geometry()` 传入 VASP 的高精度晶格，但当前
`half_library.F90::context_make_basis()` 仍然从文件构造器解析的 `self%crystal`
构造基组：

```fortran
call build_plane_wave_basis(self%crystal,self%charge%shape,self%encut,kpoint,basis)
```

因此 request geometry 虽已保存在内存中，却尚未用于文件模式的实际求解；HALF
继续采用精度较低的 CHGCAR 晶格。

## 临时规避验证

重跑目录为：

```text
/data/limusen/deepawchgs_retry_encut521/054_mp-352_HfO2
```

仅把 `ENCUT` 改为 521 eV 后，两套晶格在第 34 个 k 点均产生 3680 个平面波。
HALF 初始化全部 36 个 k 点、接受 DeepAW CHGCAR，并在 6 个 SCF `LOOP` 后达到
`EDIFF=1E-4`；作业正常结束，总耗时 132 秒。

## 正式修复要求

521 eV 只是诊断和临时规避。正式修复应：

1. 在 mapped solve 中以宿主的高精度晶格和 G 向量集合为权威输入。
2. 存在 request geometry 时，用 `request_lattice` 和 `request_grid` 构造基组。
3. PAW 原子位置、投影子及局域势装配也统一使用 request geometry。
4. 错误信息同时打印 `host_npw`、`half_npw` 和首个不同的 G 向量。
5. 增加 520 eV 回归测试：无需更改 ENCUT 即得到与 VASP 相同的 3672 个 G 向量。

