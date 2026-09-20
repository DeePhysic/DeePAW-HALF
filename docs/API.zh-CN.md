# libhalf 集成 API

`libhalf` 是供 VASP 和其他电子结构软件调用的稳定、编译器无关接口。其使用方式
类似 libxc：调用方以整数 ID 选择功能，创建不透明 context，查询能力与基组，再
通过 C ABI 执行计算。公开 ABI 用 `HALF_ABI_VERSION` 单独管理版本。

## 构建与链接

```bash
cmake -S . -B build -DHALF_ENABLE_CUDA=OFF -DHALF_BUILD_SHARED_LIBRARY=ON
cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

下游 CMake 工程直接使用导出目标：

```cmake
find_package(HALF 0.6 CONFIG REQUIRED)
target_link_libraries(mycode PRIVATE HALF::half)
```

安装内容还包括 `half.pc`、`half.h` 与可直接随调用方编译的 Fortran 接口源文件
`half_api.f90`。完整独立示例见 `examples/api/half_c_example.c` 和
`examples/cmake/CMakeLists.txt`。

CUDA 调用方最终链接时必须使用 `nvfortran -cuda`（CMake 中设
`LINKER_LANGUAGE Fortran` 并采用安装时的 CUDA 版本/CC 参数），以便 NVHPC
注册 CUDA Fortran device code；仓库示例会自动设置。纯 CPU `libhalf` 可由普通
C 编译器链接。

## 生命周期和三层接口

1. 用 `half_get_abi_version()`、`half_get_capabilities()` 检查 ABI 和能力。
2. 用 `half_create_from_files_backend()` 创建文件 context，或用
   `half_create_from_escn()` 从 DeePAW-eSCN 获取密度；运行时可选 LDA/PBE、
   CPU/CUDA/AUTO、EVD/EVJ/VASP、ENCUT 及 MIMIC_US QDEP。VASP selector 是
   CUDA 区间求解路径，求解数量取调用方传入的 `NBANDS`。简化函数
   `half_create_from_files()` 等价于 AUTO + EVD，并保持 ABI 稳定。
3. 宿主程序可用 `half_set_request_geometry()` 附加内存中已经解析好的晶格、
   分数坐标、从 1 开始的物种编号和三维实空间网格。
4. 对每个分数倒空间 k 点，用 `half_get_basis_size()` 和可选的
   `half_get_basis()` 得到 HALF 的精确 G 向量顺序。
5. 按集成深度选择：`half_assemble_hs()` 返回稠密 H/S，`half_apply_hs()` 对一组
   态施加算符，`half_solve_kpoint()` 直接返回本征对。
6. 用 `half_destroy()` 释放不透明句柄。

矩阵与态块均为 complex128、BLAS 列主序；能量单位为 eV，长度为 Angstrom，k 点
为分数倒空间坐标。把本征矢指针设为空即可只求本征值。每个函数都返回状态码，
并写入调用方提供的错误缓冲区。

`half_set_request_geometry()` 接收 `half_request_geometry_v1` 描述符，并在返回前
深拷贝所有数组，调用方随后即可释放原缓冲区。`lattice` 采用“晶格矢量、笛卡尔
分量”的行主序，`positions_fractional` 采用逐原子布局。调用
`half_get_request_geometry()` 可查询或复制保存的值。文件路径会保留这份元数据用于
审计；eSCN 构造器则直接使用同一个结构。

## DeePAW-eSCN 远程密度

0.6 版提供两层远程推理接口：

- `half_escn_health()` 和 `half_escn_predict()` 直接封装 HTTP 服务。调用方传入
  原子序数、Å 单位笛卡尔坐标、行主序晶胞和 `[nx,ny,nz]`；返回服务器原始 C-order
  float32 密度，并可选返回 `nu`、`alpha`、`beta`、`risk`。
- `half_create_from_escn()` 接收现有 `half_request_geometry_v1`，从匹配 POTCAR
  自动得到原子序数和价电子数，调用 `/v1/predict`，把 C-order 网格转换成 HALF
  内部兼容 CHGCAR 的顺序，并创建可直接求解的普通 context。

HALF 数值路径应选择 `HALF_DENSITY_NORMALIZE_VALENCE`，它把预测网格均值显式
归一到 POTCAR ZVAL 总和。只有模型输出已经满足下游约定时才使用
`HALF_DENSITY_RAW`。两种模式都不会截断负密度。

`HALF_CAP_ESCN_API` 表示该能力可用。CMake 优先使用 libcurl；Unix 环境没有
libcurl 开发包时，会自动采用无外部依赖的 HTTP transport。内置 transport 支持
SSH 隧道所需的 `http://`，HTTPS 需要 libcurl。可用
`-DHALF_ENABLE_ESCN_API=OFF` 显式关闭。

按服务文档建立隧道后可直接验证：

```bash
ssh -N -L 8265:127.0.0.1:8265 cmu-pro6000
./build/cpu-release/half-escn-example http://127.0.0.1:8265
```

`half_escn_predict()` 会检查 2,000,000 网格点和 2 MiB 请求限制。高层构造器只
请求密度，因为 HALF Hamiltonian 当前不消费不确定性；需要 NIG 数组的程序应调用
低层函数并设置 `HALF_ESCN_INCLUDE_UNCERTAINTY`。

CUDA 构建中 AUTO 选择 CUDA；`half_solve_kpoint()` 的势构造、MIMIC_US/QDEP、
H/S 组装和广义本征求解全部在 GPU 上执行，只把请求的本征对传回。ABI v1 的
主机 H/S 导出与 H/S 应用目前是 CPU 接口；CUDA context 对这两个调用返回
`HALF_ERROR_UNAVAILABLE`。后续可在不改变已有符号的情况下增加 device-pointer
扩展。CPU 构建中 AUTO 选择 oneMKL；EVJ 和 VASP 部分谱模式仅支持 CUDA，EVD
同时支持两者。

## 接入 VASP 的方式

最小 VASP 补丁在固定平滑密度和匹配 POTCAR 就绪后创建一次 context，并把
`LATT_CUR%A`、`T_INFO%POSION`、`T_INFO%ITYP`、`GRIDC%NGPTAR` 从内存直接传给
HALF。对每个 VASP k 点，先查询 `npw` 与 G 向量映射，只做一次系数顺序转换，然后
可选择：

- 使用 CPU context，在 VASP 的迭代本征求解器中调用 `half_apply_hs()`；或
- 调用 `half_solve_kpoint()`，直接使用 HALF 的 CPU 或全 GPU 稠密 EVD。
- 自带平面波排序的宿主程序可调用 `half_solve_kpoint_mapped()` 并传入
  `(Gx,Gy,Gz)`；HALF 会验证基组完全一致，再按宿主顺序返回系数。VASP 接入见
  [VASP_INTEGRATION.zh-CN.md](VASP_INTEGRATION.zh-CN.md)。

同一 context 可跨所有 k 点复用。密度、结构、POTCAR、ENCUT、XC 或后端变化时应
销毁并重建。每个 MPI 进程使用独立 context，k 点分发仍由调用软件负责；同一个
context 不应被两个线程并发执行可变操作。

文件构造器读取 CHGCAR/vaspwave.h5 与 POTCAR（也可用 vaspwave.h5 内嵌的
POTCAR）。远程构造器使用现有 VASP geometry descriptor 与 POTCAR，不需要修改
`vasp2half` 结构。VASP adapter 通过 INCAR 显式选择远程路径：

```text
LHALF_INIT = .TRUE.
LHALF_API = .TRUE.
HALF_MODE = VASP_LIKE
HALF_ESCN_URL = http://127.0.0.1:8265
```

`HALF_MODE=TRADITIONAL` 是默认值。`HALF_MODE=VASP_LIKE` 选择 CUDA 部分谱求解器，
并通过现有 mapped k-point ABI 传入 VASP 的 `NBANDS`，使 HALF 只返回所需的最低
本征对。`DENSE` 和 `VASP` 是兼容简写。

`LHALF_API=.FALSE.` 为默认值，此时即使环境中存在 `HALF_ESCN_URL`，HALF 也
继续使用 CHGCAR。启用 `LHALF_API` 后，INCAR 中的 URL 优先；只有 INCAR URL
为空时才把同名环境变量作为兼容回退。两处都没有地址会直接报错，CHGCAR 缺失
不会隐式触发联网。当前 VASP 集成中，这只切换 HALF 内部用于重构初始波函数的
密度；VASP 自己的 `ICHARG=1` host density 仍从 CHGCAR 读取。后者若也要取消
文件输入，需要另做 VASP charge-grid 内存注入，不能隐含在此构造器中。
