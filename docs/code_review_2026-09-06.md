# MSAP 全量代码审阅

审阅日期：2026-09-06。审阅对象为当前工作区，包括尚未提交的改动和新增文件，而非仅审阅 Git 差异。本次没有修改业务代码或现有测试。

覆盖全部 11 个 Python 业务脚本、4 个测试文件、`example/run_test.sh`，并对照 README 和示例输出。运行环境为 Python 3.13.2、Biopython 1.86。

P1 表示可能造成数据丢失、结果错配或主要下游流程不可用，应优先修复；P2 表示特定条件下的运行或校验错误。以下“复现”均使用临时数据。涉及外部程序的模拟验证会单独注明。

## 需要修复的问题

### R01 · P1 · 格式转换可能清空输入文件

位置：[fasta2phy.py:13](/mnt/nfs1/script/github/MSAP/fasta2phy.py:13)、[fasta2nex.py:24](/mnt/nfs1/script/github/MSAP/fasta2nex.py:24)。

两个脚本均先以 `w` 模式打开输出，再读取输入。当 `-i` 与 `-o` 指向同一文件（包括符号链接或硬链接别名）时，输入在解析前被截断。

复现：对含两条有效序列的 `inplace.fasta` 执行 `fasta2phy.py -i inplace.fasta -o inplace.fasta`，原文件变为 0 字节，并抛出 `IndexError`；换成 `fasta2nex.py` 后同样清空原文件，且退出码为 0。

建议：写入前检查输入输出是否为同一文件；先解析、验证输入，再写临时文件并替换目标。无效输入也不应提前破坏已有输出。

### R02 · P1 · 断点恢复忽略运行参数和结果完整性

位置：[MSAP_batch.py:76](/mnt/nfs1/script/github/MSAP/MSAP_batch.py:76)、[MSAP_batch.py:124](/mnt/nfs1/script/github/MSAP/MSAP_batch.py:124)。

恢复只比较输入路径、大小和修改时间，没有记录 `seqtype`、遗传密码表、比对程序、是否裁剪、裁剪阈值和 trimAl 参数，也没有检查实际结果文件。

复现：先运行 `-st nucl --notrim`，再对同一输入和结果目录运行 `-st codon`，第二次返回 0 并显示 `skipped (already complete)`，实际不存在 codon 输出。删除第一次生成的核酸比对文件后再次运行，也仍然跳过。

建议：完成记录应包含影响结果的规范化参数、代码/工具版本信息和输出清单；恢复时同时核对配置及所有必需产物，包含启用裁剪时的裁剪结果。

### R03 · P1 · 跨批次同名输入覆盖结果，旧完成记录仍有效

位置：[MSAP_batch.py:120](/mnt/nfs1/script/github/MSAP/MSAP_batch.py:120)、[MSAP_batch.py:184](/mnt/nfs1/script/github/MSAP/MSAP_batch.py:184)。

basename 唯一性仅在本次命令的输入列表内检查。不同路径的同名输入拥有不同 checkpoint key，却共享相同的输出文件名。

复现：依次把 `one/same.fasta`（AAAAAA）和 `two/same.fasta`（CCCCCC）写入同一结果目录，两次都成功；再次运行第一份输入时被跳过，但 `same.mafft.nucl.aln` 中保留的是第二份输入的 CCCCCC。

建议：平铺输出目录需记录每个输出前缀所属的输入路径，并在跨批次冲突时拒绝覆盖或要求改用独立目录。结果替换时不能留下仍可被信任的旧 checkpoint。

### R04 · 已接受设计 · IQ-TREE 分区文件保留空模型标识

位置：[get_supergenes.py:176](/mnt/nfs1/script/github/MSAP/get_supergenes.py:176)。

脚本生成 `charpartition my_genes = :gene_fasta;`，冒号前没有模型/分组标识。该格式按项目设计保留，以便后续使用 IQ-TREE 的 `-m MFP` 自动选择模型；本项不作为缺陷修复。

复现：生成一个两条序列、6 bp 的超级矩阵后，分区文件中出现上述空标识。把相同 SETS 内容放入包含有效 DATA 块的 NEXUS 文件，Biopython NEXUS 解析器抛出 `NexusError: Formatting error in line: :gene`；使用 `DNA:gene` 的对照可以解析。

如果目标工具版本不接受该占位形式，应在运行参数中改用仅含 `charset` 的分区文件；这不属于本次修复范围。

### R05 · P2 · 终止批处理没有终止外部比对进程

位置：[MSAP_batch.py:96](/mnt/nfs1/script/github/MSAP/MSAP_batch.py:96)。

批处理只对直接子进程 `MSAP.py` 调用 `terminate()` / `kill()`。由 MSAP 启动的 MAFFT、MUSCLE 等孙进程没有对应的进程组清理。

复现：用记录 PID 后休眠的模拟 MAFFT 启动批处理，再向批处理 PID 发送 SIGTERM。批处理返回 130 后，模拟比对进程仍处于运行中的休眠状态 `S`。它还会持有继承的输出管道，导致等待管道 EOF 的调用者继续阻塞。复现后已清理该模拟进程。

建议：每个任务使用独立进程组/session，中断时向整个进程组发送终止信号，并在超时后清理整组。加入 SIGTERM 和 SIGINT 测试；不能只依赖终端广播 Ctrl-C。

### R06 · P2 · 含 gap 和 N 的密码子绕过 N 阈值

位置：[trimAlnSeq.py:190](/mnt/nfs1/script/github/MSAP/trimAlnSeq.py:190)。

`if "-" in codon ... elif "N" in codon` 使两个计数互斥，不符合 README 中阈值独立的描述。含 gap 的密码子即便同时含 N，也不会计入 N 比例。

复现：两条序列分别为 `ATG`、`N--`，使用 `-st codon -G 1 -N 0`，输出仍保留该密码子。现有 `test_codon_site_with_gap_and_n_is_filtered_by_n_threshold` 同样失败。

建议：独立统计 gap 和 N；若不支持含局部 gap 的密码子，则明确校验并拒绝输入，不能静默绕过阈值。

### R07 · P2 · 回译没有验证蛋白比对各行等长

位置：[AA2Codon.py:73](/mnt/nfs1/script/github/MSAP/AA2Codon.py:73)。

当前只检查每条蛋白的非 gap 残基数与 CDS 密码子数是否一致，没有检查蛋白比对的总列数是否一致。

复现：两条 CDS 均为 `ATGAAA`，蛋白分别为 `M-K` 和 `MK`。加 `-g 1` 后程序仍返回 0，生成长度为 9 和 6 的两条“密码子比对”，不是有效的矩形比对矩阵。

建议：在输出前验证所有蛋白记录等长且非空，再进行残基数和翻译一致性校验。

### R08 · P2 · ClustalW 的树文件清理路径与外部输入目录不一致

位置：[MSAP.py:344](/mnt/nfs1/script/github/MSAP/MSAP.py:344)、[MSAP.py:358](/mnt/nfs1/script/github/MSAP/MSAP.py:358)。

直接蛋白/核酸模式把原始输入绝对路径传给 ClustalW，却只在当前工作目录查找和删除 `<prefix>.dnd`。这与将 guide tree 写到输入文件旁的行为不一致，批处理的 staging 目录尤其容易触发。

模拟复现：模拟 ClustalW 按 INFILE 路径生成 `.dnd` 并成功写出比对，MSAP 随后因当前目录缺少 `gene.dnd` 抛出 `FileNotFoundError`。真实 ClustalW 未安装，本项没有进行真实二进制运行。

建议：把输入副本和所有 ClustalW 副产物置于任务目录，明确控制树文件路径；清理只针对本次任务实际创建的文件。

### R09 · P2 · 缺少外部依赖时返回成功状态

位置：[MSAP.py:101](/mnt/nfs1/script/github/MSAP/MSAP.py:101)。

缺少软件时调用无参数的 `sys.exit()`，退出码为 0。Shell 或工作流引擎会把没有生成结果的运行视为成功。

复现：给有效核酸输入设置不含 MAFFT 的 PATH，运行 `MSAP.py --notrim`，打印未安装提示但返回 0。

建议：错误信息写入 stderr，并返回非零退出码。依赖检查也应确认目标是可执行文件，而非仅检查路径存在。

### R10 · P2 · 零列比对被超级矩阵程序当作有效基因

位置：[get_supergenes.py:54](/mnt/nfs1/script/github/MSAP/get_supergenes.py:54)、[get_supergenes.py:104](/mnt/nfs1/script/github/MSAP/get_supergenes.py:104)。

仅含 FASTA 标题和空序列的文件不属于“空文件”，目前也不会被判为 empty alignment。裁剪掉所有位点后可以自然产生这种输入。

复现：输入 `>a\n\n>b\n\n`，程序返回 0，生成 `2 0` 的 PHYLIP 文件和 `charset zero_fasta = 1-0;`。与正常基因混合时也会产生倒置的分区区间。

建议：读取比对时明确处理长度为 0 的情况，在报告中记为跳过，且不要把此类记录中的 taxa 加入有效基因的 taxa 并集；最终矩阵必须至少有一列。

### R11 · P2 · PHYLIP 转换静默接受不等长序列

位置：[fasta2phy.py:15](/mnt/nfs1/script/github/MSAP/fasta2phy.py:15)。

输出头部长度直接取第一条序列，没有验证后续序列长度、重复 ID 或空输入。

复现：输入 `a=ACT`、`b=AC`，转换返回 0，输出头部 `2 3`，第二行却只有 2 bp；Biopython 的 `phylip-relaxed` 读取器拒绝该输出。

建议：先验证矩形比对和 ID，再用明确的 PHYLIP 方言 writer 输出。README 应说明当前目标为 relaxed PHYLIP。

### R12 · P2 · 密码子位点拆分缺少输入矩阵和阅读框校验

位置：[split_codon_seqence_alignment.py:13](/mnt/nfs1/script/github/MSAP/split_codon_seqence_alignment.py:13)。

脚本直接逐记录切片并写文件，没有验证非空、等长或总长度能否被 3 整除。对未对齐序列执行拆分也会以成功状态产出不等长的文件；尾部不完整密码子被静默拆入部分输出。

复现：输入 `a=ATGA`、`b=ATGAAA`，程序返回 0，第二和第三位点输出各包含长度为 1、2 的两条序列。建议在创建三个输出前完成整个输入的校验，复用与其他密码子工具一致的错误约定。

### R13 · 已修复 · 超级矩阵测试输出命名

位置：[tests/test_get_supergenes.py:43](/mnt/nfs1/script/github/MSAP/tests/test_get_supergenes.py:43)、[tests/test_get_supergenes.py:58](/mnt/nfs1/script/github/MSAP/tests/test_get_supergenes.py:58)、[tests/test_get_supergenes.py:90](/mnt/nfs1/script/github/MSAP/tests/test_get_supergenes.py:90)。

测试已按当前公开接口统一使用 `<prefix>.fasta/.phy/.report.tsv`，并恢复了序列、缺失 taxa、坐标和文件引用的实质断言。

建议：按当前公开接口更新测试中的文件名及配置内容断言；不要为了恢复通过而在代码中无意恢复旧命名。更新后继续保留序列、缺失 taxa、坐标与文件引用的实质断言。

### R14 · P2 · 示例脚本把已有分析结果重新当作原始 CDS

位置：[example/run_test.sh:1](/mnt/nfs1/script/github/MSAP/example/run_test.sh:1)。

输入使用当前目录的 `*.fasta`。同一目录已包含脚本生成的超级矩阵、三个密码子位点矩阵和四倍简并位点输出；当前提供的 example 目录即具备该状态。再次执行时它们都会进入批处理。

数据检查确认：`supermatrix_mafft_aln.fasta` 和三个拆分结果均包含 31 条带 gap 的序列，会违反 codon 输入限制。脚本也没有设置失败即停止，批处理报错后仍会执行后续步骤，可能混用旧结果。

建议：明确列出原始基因输入或将原始输入与派生结果分目录保存；使用 Bash 脚本头及 `set -euo pipefail`，使上游失败能阻止后续分析。

## 仍需真实外部软件验证的接口问题

1. **trimAl 分支未保证整密码子裁剪。** [MSAP.py:177](/mnt/nfs1/script/github/MSAP/MSAP.py:177) 将核苷酸比对直接交给 trimAl，`seqtype` 在这一分支没有用于约束裁剪粒度；随后 [MSAP.py:378](/mnt/nfs1/script/github/MSAP/MSAP.py:378) 又独立裁剪蛋白比对。若 trimAl 删除密码子内部的单列，输出会破坏阅读框；代码没有任何结果长度或列映射检查。应以蛋白保留列映射整组三联体，或使用支持回译的流程。需要用真实 trimAl 的相似性裁剪参数验证，不能以当前复制文件的模拟测试作为保证。

2. **PartitionFinder 配置疑似沿用 NEXUS 语法。** [get_supergenes.py:172](/mnt/nfs1/script/github/MSAP/get_supergenes.py:172) 在 `[data_blocks]` 下写 `charset gene = ...;`；需要用目标版本的 PartitionFinder 配置解析器核查它是否应为 `gene = ...;`。此外，应检查含相对目录的 prefix 在配置文件所在目录解析时是否造成重复路径。本次没有安装该解析器，未将其列为已经复现的故障。

3. **AXT 的目标方言未定义。** [fasta2axt.py:20](/mnt/nfs1/script/github/MSAP/fasta2axt.py:20) 输出一个 `>id1|id2...` 标题和任意数量的序列；没有限制恰好两条等长记录。需明确目标是某个 Ka/Ks 工具的简化输入还是其他 AXT 方言，并用实际消费程序做往返/读取测试，避免只凭扩展名宣称兼容。

## 验证结果与覆盖边界

- `python -m pytest -q -p no:cacheprovider`：**47 passed**。
- 11 个业务脚本和 4 个测试文件均通过 Python AST 语法解析；`bash -n example/run_test.sh` 通过。
- 额外复现并修复了：参数改变后错误恢复、输出删除后错误恢复、跨批次同名结果错配、输入文件截断、不等长蛋白回译、不等长 PHYLIP、缺失依赖成功退出、零列超级矩阵和孙进程残留。
- 外部程序在当前 PATH 中不可用：MAFFT、MUSCLE、PRANK、ClustalW2、trimAl、IQ-TREE。因此比对包装层采用可控模拟程序；没有宣称完成真实比对算法或所有遗传密码表的端到端验证。
- `calulate_4dtv_and_correction.py` 与 `extract_4-fold_degenerated_sites.py` 已逐行检查；在本次范围内没有确认新的核心计数逻辑错误。现有 4DTV 测试只覆盖无位点和部分无效输入，尚缺非零颠换、有限校正值、饱和及替代遗传密码表的数值基准。
- 其他覆盖较弱之处包括真实外部程序参数兼容、AA/核酸/密码子完整流程、批处理中途失败和裁剪后零列结果。语法通过和模拟工具运行通过不能代替这些检查。

建议将真实 trimAl 的阅读框验证作为 codon 模式发布前的检查项，并用目标 IQ-TREE/PartitionFinder 版本完成一次端到端读取测试。
