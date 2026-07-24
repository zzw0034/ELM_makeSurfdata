# Pathfinder Slurm 速查：`-p` / `-q` / `--mem` 怎么配

> 面向 ORNL Pathfinder、账户 `hpcl-cli185`（用户 zw5）。
> 数据用 `sinfo` / `sacctmgr` 实测于 2026-07-09，如集群扩容请重新核对。
> 相关：`esmf_pathfinder_setup.md`（ESMF/mapping 作业完整设置）。

---

## 0. 一句话

`-p`（分区）决定用哪批节点，`-q`（QOS）必须和分区配套，`--mem`（每节点内存）
是节点门槛——**写得越小，能排上的节点越多**。Pathfinder 强制显式写
`-q -n -c --mem`，少一个提交就被拒。

---

## 1. 你能用的四种 `-p` + `-q` 组合

账户授权了两个 QOS：`hpcl-cli185`（专属）和 `normal`（公共）。可用组合：

| `-p` 分区 | 配套 `-q` | 节点 | 墙钟上限 | 排队 | 用途 |
|---|---|---|---|---|---|
| `hpcl-cli185` | `hpcl-cli185` | 20 × 128核 × 478GB | **无限制** | 几乎不排队 | ✅ 首选，mapping/日常 |
| `parallel` | `normal` | 206（混合，见 §3） | **24h** | 看节点档位 | 多节点大作业 |
| `serial` | `normal` | 同上 206 | 24h | 看节点档位 | 单节点任务 |
| `gpu` | `normal` | 38 × A2 GPU | 24h | — | GPU，单节点，加 `--gres=gpu:N`（N≤6） |

**配对铁律**：专属分区配专属 QOS，公共分区配 `normal`。
`-p hpcl-cli185 -q normal` 会被拒（实测）。

官方文档只列公共的 serial/parallel/gpu + normal；`hpcl-*` 是 condo 专属分区，
文档不列，但你的账户可用。

---

## 2. `--mem` 基本规则

- 单位：`--mem=470g` 或 `--mem=481280`（MB）。是**每节点**值，不是总量。
- 含义：Slurm 只把作业放到「可分配内存 ≥ 该值」的节点上。
- 必填：Pathfinder 不接受 `--mem=0`（Baseline 老写法）。
- **总内存 = 节点数 × --mem**。同样的总内存，可以少节点多要，也可多节点少要。

---

## 3. parallel 分区的内存档位（关键）

parallel 的 206 个节点是**混合规格**，`--mem` 有四个临界点，跨过一个就少一批节点：

| `--mem` ≤ | 可用节点数 | 累计包含的节点规格 |
|---|---|---|
| `220g`  (231906 MB) | **206（全部）** | 226GB×68 + 478GB×72 + 755GB×38 + 2.3TB×28 |
| `470g`  (489954 MB) | 138 | 478GB×72 + 755GB×38 + 2.3TB×28 |
| `750g`  (773052 MB) | 66  | 755GB×38 + 2.3TB×28 |
| `2200g` (2321280 MB) | 28  | 2.3TB×28（大内存节点） |

**别顶格写**：门槛是「≥」，写的值不能超过节点标称可分配内存。
例：755GB 档精确上限 773052 MB，写 `755g`(=773120) 就超 68 MB，这批节点全失效。
所以每档留几 GB 余量：想要 478GB 档写 `470g`，想要 755GB 档写 `750g`。

---

## 4. 想用更多节点 → 把 `--mem` 调小

作业要的是总内存。少节点每个多要 = 只能挑大节点、排队久；
多节点每个少要 = 几乎全池可选、排队快。前提是程序能把内存摊到更多节点
（MPI 程序如 E3SM、ESMF_RegridWeightGen 都可以）。

例（E3SM spinup，总内存需求 ~8TB）：

```bash
# A) 只能等 28 个大内存节点
#SBATCH -N 8
#SBATCH --mem=1024gb        # 8 × 1TB = 8TB，仅 2.3TB 档满足

# B) 等价总量，138 个节点可选，排队快得多
#SBATCH -N 18
#SBATCH --mem=470g          # 18 × 470GB ≈ 8.5TB
```

⚠️ 改节点数不是只改 SBATCH 头：E3SM 要重新 `case.setup` 调整 PE layout，
MPI 作业要相应改 `-n` / layout。

---

## 5. 专属分区没有档位问题

`-p hpcl-cli185` 的 20 节点规格统一（478GB），`--mem=470g` 就是唯一答案，
不用纠结档位。§3/§4 的档位思维只在去 `parallel` 借公共节点时才用到。

---

## 6. 可直接抄的模板

**mapping / 日常（专属，首选）：**

```bash
#SBATCH -p hpcl-cli185
#SBATCH -A hpcl-cli185
#SBATCH -q hpcl-cli185
#SBATCH --nodes=10
#SBATCH -n 120
#SBATCH -c 1
#SBATCH --mem=470g
#SBATCH -t 4:00:00
```

**多节点大作业（公共，需大内存或 >20 节点）：**

```bash
#SBATCH -p parallel
#SBATCH -A hpcl-cli185
#SBATCH -q normal
#SBATCH --nodes=18
#SBATCH -n <总任务数>
#SBATCH -c 1
#SBATCH --mem=470g          # 见 §3 选档位
#SBATCH -t 24:00:00         # normal 硬上限；更长靠 restart 分段
```

---

## 7. 自查命令

```bash
sacctmgr show associations user=zw5 format=account,qos      # 我能用哪些 QOS
sacctmgr show qos format=name,MaxWall,MaxTRESPU             # 各 QOS 墙钟/资源上限
shownodes parallel                                          # 分区节点实时状态
sinfo -N -p parallel -o "%m %c" | sort | uniq -c            # 各内存档位节点数
sbatch --test-only job.sh                                   # 提交前校验（不真跑）
```
