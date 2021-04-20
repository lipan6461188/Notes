## 测试hhblits速度



测试的代码为：

```shell
hhblits=(/Share2/home/zhangqf/lipan/usr/hhsuite-3.3.0-SSE2-Linux/bin/hhblits /Share2/home/zhangqf/lipan/usr/hhsuite-3.3.0-AVX2-Linux/bin/hhblits /Share2/home/zhangqf/lipan/usr/hh-suite-master/build/bin/hhblits)
cpu=(1 5 10)

for id1 in $(seq 0 2); do
    for id2 in $(seq 0 2); do
        echo "=======> ${hhblits[$id1]} CPU=${cpu[$id2]} <=========="
        cmd="time ${hhblits[$id1]} \
            -i /BioII/zhangqf7/lipan/CPI/AlphaFold1/prepare_input/4hhjA01/input.fa \
            -d /150T/zhangqf/GenomeAnnotation/UniClust30/UniRef30_2020_06 \
            -o hhblits.hhr \
            -ohhm hhblits.hmm \
            -opsi hhblits.psi \
            -oa3m hhblits.a3m \
            -n 3 \
            -e 0.001 \
            -v 0 \
            -cpu ${cpu[$id2]}"
        eval $cmd
    done
done
```

三种hhsuite所在目录为：

```
/Share2/home/zhangqf/lipan/usr/hhsuite-3.3.0-SSE2-Linux/bin/hhblits
/Share2/home/zhangqf/lipan/usr/hhsuite-3.3.0-AVX2-Linux/bin/hhblits
/Share2/home/zhangqf/lipan/usr/hh-suite-master/build/bin/hhblits
```

使用MPI测试：

```shell
/Share2/home/zhangqf/usr/anaconda/bin/mpirun \
	-np 4 /Share2/home/zhangqf/lipan/usr/hh-suite-master/build/bin/hhblits \
		-i /BioII/zhangqf7/lipan/CPI/AlphaFold1/prepare_input/4hhjA01/input.fa \
    -d /150T/zhangqf/GenomeAnnotation/UniClust30/UniRef30_2020_06 \
    -o hhblits.hhr \
    -ohhm hhblits.hmm \
    -opsi hhblits.psi \
    -oa3m hhblits.a3m \
    -n 3 \
    -e 0.001 \
    -v 0 \
    -cpu [1,5,10]
```

在ZIO01节点上测试的结果如下：

| 程序                                 | -cpu 1    | -cpu 5    | -cpu 10   |
| ------------------------------------ | --------- | --------- | --------- |
| hhsuite-3.3.0-SSE2                   | 4m35.096s | 2m27.720s | 2m4.295s  |
| hhsuite-3.3.0-AVX2                   | 3m18.960s | 1m56.787s | 2m19.234s |
| hhsuite-3.3.0-Compile                | 4m43.762s | 2m15.691s | 2m13.547s |
| hhsuite-3.3.0-Compile + MPI (-np 4)  |           | 2m35.929s | 2m31.772s |
| hhsuite-3.3.0-Compile + MPI (-np 10) |           |           |           |

