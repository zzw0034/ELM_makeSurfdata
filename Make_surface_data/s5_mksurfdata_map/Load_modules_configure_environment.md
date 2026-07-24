# Load the relevant modules and configure the environment

本文是 **Pathfinder (ORNL, `pflogin*.ornl.gov`)** 上编译、运行 mksurfdata_map 的环境配置说明。只需用下面**一节「一键复制」**即可；其余为与教程对照、排错与注意事项。

> **已在 Pathfinder 验证（2026-07-12）**：下面的一键复制可直接编出可运行的二进制。验证过的完整脚本见
> `<E3SM>/components/elm/tools/mksurfdata_map/build_pathfinder.sh`（`bash build_pathfinder.sh` 一键编译）。

---

## 一键复制（推荐）

打开终端，**整段复制粘贴执行**即可。

```bash
# 1) 编译器与 MPI
module purge
module load gcc/12.4.0
module load openmpi/5.0.5

# 2) I/O 库
module load hdf5/1.14.5-mpi
module load netcdf-c/4.9.2-mpi-h5f
module load netcdf-fortran/4.6.1-mpi-h5f

# 3) NETCDF_HOME / HDF5_HOME
#    Pathfinder 上这套 module 不设 NETCDF_DIR、HDF5_DIR、OLCF_* 等变量，
#    所以直接用 nf-config / nc-config 取安装路径（已验证解析到正确的 spack 目录）。
export NETCDF_HOME=$(nf-config --prefix)
export HDF5_HOME=$(nc-config --prefix)

# 4) Makefile 用
export LIB_NETCDF=$NETCDF_HOME/lib
export INC_NETCDF=$NETCDF_HOME/include
export USER_FC=mpifort
export USER_CC=mpicc
export USER_LDFLAGS="$(nc-config --libs) $(nf-config --flibs)"
# gfortran 12 编这套较老的 E3SM 源码需要三个兼容 flag（逐个报错试出，见下方排错表）：
#   -fallow-invalid-boz       : nanMod.F90 的八进制 BOZ 字面量
#   -fallow-argument-mismatch : mkncdio.F90 里 nf_def_var/nf_put_var_int 一次传数组一次传标量（rank 不匹配）
#   -ffree-line-length-none   : mkfileMod.F90 属性字符串超过 132 列被截断
export USER_FFLAGS="-fallow-invalid-boz -fallow-argument-mismatch -ffree-line-length-none"
```

编译：

```bash
cd <E3SM>/components/elm/tools/mksurfdata_map/src
gmake clean
gmake -j 8
```

（登录节点上 `gmake -j 8` 即可，构建很快。生成的可执行文件在上一级目录：`../mksurfdata_map`。）

---
---

## 排错与分步理解（可选）

标注说明：**✅ 本次实测** = 2026-07-12 在 Pathfinder 编译时亲自撞到、导致编译停下的错误；**通用** = 一般性排查建议或被 flag 提前挡掉的问题，本次并未作为错误遇到。

- **✅ 本次实测 — mkncdio.F90: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)**：`nf_def_var`（及 `nf_put_var_int`）在同一文件里一次用数组、一次用标量做同一个参数，gfortran 10+ 默认当**错误**。加 `-fallow-argument-mismatch` 降为警告即可（一键复制已含）。这是第 1 次编译（只带 `-fallow-invalid-boz`）真正停下的地方。
- **✅ 本次实测 — mkfileMod.F90: Line truncated [-Werror=line-truncation] / Unterminated character constant**：自由格式源码里有超过 132 列的长属性字符串，被截断导致字符常量未闭合。加 `-ffree-line-length-none` 取消行长限制（一键复制已含）。这是第 2 次编译停下的地方。
- **通用 — nanMod.F90: BOZ literal constant ... is neither a data-stmt-constant nor an actual argument to INT, REAL, DBLE, or CMPLX**：gfortran 12 对八进制/十六进制字面量更严格。本次因一键复制里从头就带了 `-fallow-invalid-boz`，它只以 **warning** 出现、未阻塞编译；若去掉该 flag 会报成错误。若要改源码：把 `inf`/`nan`/`bigint` 的 BOZ 改为标准写法，例如 `TRANSFER(INT(O'0777600000000000000000', KIND=8), 1.0_r8)`（按实际类型与位数调整）。
- **通用 — 确认 module 是否生效**：执行 `which nc-config nf-config`、`nc-config --prefix`、`nf-config --prefix`。两 prefix 不同（netcdf-c 与 netcdf-fortran 分装），用 nf-config 作 NETCDF_HOME 即可（一键复制已体现）。本次未出问题。
- **通用 — 编译报错 `cannot find -l...` 或 `undefined reference`**：多为 netcdf-fortran 未 load 或 nf-config 不在 PATH；先确认第 2 步三个 I/O 库 module 都已 load。本次链接一次即通，未遇到。
- **小结**：本次实测只需在原有 `-fallow-invalid-boz` 基础上再补 `-fallow-argument-mismatch`、`-ffree-line-length-none` 两个 flag 即可编过；三者是这套源码在 gfortran 12 下的完整集合，缺一个就会在对应文件停下。若换更老的编译器或更新过源码，可能不全需要。

---

## 注意事项

- 使用 **hdf5/1.14.5-mpi**，不要用 serial 版 `hdf5/1.14.5`，否则与 netcdf-*-mpi-h5f 库不匹配。
- **编译**在登录节点做即可（小活儿）。上述环境变量只在**当前 shell**有效；换 shell 或提交作业脚本时，要在脚本里重新 `module load` 那几行（编译期的 USER_* flag 则不必再设，二进制已经编好）。
- **Pathfinder 运行时注意（cgroup / `/dev/shm`）**：Pathfinder 的 Slurm 用 cgroup 把 `/dev/shm` 用量计入作业内存配额，MPI 默认走共享内存传输时可能被 SIGBUS 杀掉。若以 MPI 并行方式运行 mksurfdata_map（或任何 MPI 程序），在作业脚本里加
  `export UCX_TLS=rc,ud,dc,self,tcp`
  强制通信走 InfiniBand + TCP、禁用共享内存传输。
- 配置好后，按 `your_workflow_status_and_step5.md` 的 5.2–5.6 继续（生成 namelist、设 mksrf_fdynuse、运行 mksurfdata_map）。
