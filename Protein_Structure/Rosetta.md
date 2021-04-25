### 加载Rosetta的环境变量

```shell
export ROSETTA=/Share/home/zhangqf/lipan/usr/rosetta/rosetta_bin_linux_2019.12.60667_bundle
export RNA_TOOLS=$ROSETTA/tools/rna_tools
export PATH=$RNA_TOOLS/bin/:$ROSETTA/main/source/bin:$PATH
export PYTHONPATH=$PYTHONPATH:$RNA_TOOLS/bin
export ROSETTA3=$ROSETTA/main
export ROSETTA3_SRC=$ROSETTA/main/source
export ROSETTA3_MAIN=$ROSETTA/main
export ROSETTA3_DB=$ROSETTA/main/database
```

### Rosetta的命令行

https://new.rosettacommons.org/docs/latest/full-options-list

```
-hb_cen_soft 							# Use softer version of cen_hb term (Default: false)
-relax:default_repeats 5 	# Default number of repeats done by FastRelax. Has no effect if a custom script is used!
-default_max_cycles 200 	# max cycles for MinimizerOptions (Default: 2000)
-out:level 100 						# Specified hierarchical mute levels for individual channels in following format: -levels all:300 core.pose:500. Numeric values could be substituted with mute level names like: debug, info, error etc. Please note that all: is synonymous to -level:
```

### 打分函数

https://www.rosettacommons.org/demos/latest/tutorials/scoring/scoring

```shell
# 由于PDB文件中国存在一个P，所以使用ignore_unrecognized_res忽略这个原子
score_jd2 -in:file:s 3tdm.pdb -ignore_unrecognized_res
# -out:pdb表示输出实际打分的PDB文件
```

分数存储在一个叫做`score.sc`的文件中，PDB存在`3tdm_0001.pdb`文件中。由于Rosetta随机地重构丢失的侧链，所以每次运行都可能会产生不同的分数和PDB文件。注意，score_jd2的输出每次都是append到`score.sc`文件中。

```shell
cd /Share/home/zhangqf/lipan/usr/rosetta/rosetta_bin_linux_2019.12.60667_bundle/demos/tutorials/scoring
score_jd2 @flag_docking # @表示读入文件内容作为参数
```

文件flag_docking中的内容是：

```
-in:file:l input_files/pdblist

-score:weights docking

-out:file:scorefile output_files/score_docking.sc
```

ref2015是一套常用的全原子打分函数参数，所在目录为：

```
/Users/lee/anaconda3/lib/python3.7/site-packages/pyrosetta/database/scoring/score_functions/hbonds/ref2015_params
```

在其他特殊条件下，可以使用其他特殊的打分函数。ref2015中的分数主要分为：
* 物理力：静电相互作用力、范德华相互作用
* 基于统计得到的项：来自Ramachandran空间的Torsion angles概率

talaris2013和talaris2014的权重文件：
https://www.rosettacommons.org/demos/latest/tutorials/scoring/scoring

