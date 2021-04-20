<center><h1>AlphaFold1数据准备</h1></center>

<center>李盼 2021-04-11</center><br>

### 1. 训练数据

> Our models are trained on structures extracted from the PDB. We extract non-redundant domains by utilizing the <mark>CATH 35% sequence similarity cluster representatives</mark>. This generated 31,247 domains, which were split into train and test sets (29,427 and 1,820 proteins, respectively), keeping all domains from the same homologous superfamily (H-level in the CATH classification) in the same partition. The CATH superfamilies of FM domains from CASP11 and CASP12 were also excluded from the training set. From the test set, we took—at random—a single domain per homologous superfamily to create the 377 domain subset used for the results presented here. We note that accuracies for this set are higher than for the CASP13 test domains.

CATH的数据可以从[这里](http://www.cathdb.info/wiki/doku/?id=data:index)下载，文件`cath-b-s35-newest.gz`包含了训练的数据：

```shell
wget ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/cath-b-s35-newest.gz -O cath-b-s35-newest.gz
gzip -d cath-b-s35-newest.gz
wc -l cath-b-s35-newest
#    32694 cath-b-s35-newest
# 因此我下载的文件包含了32694个domain
```

注意作者下载的版本是： `CATH 16 March 2018`。`cath-b-s35-newest`文件的内容是：

```
1oaiA00 v4_2_0 1.10.8.10.1 561-619:A
3frhA01 v4_2_0 1.10.8.10.2 -1-56:A
4g3oA00 v4_2_0 1.10.8.10.3 446-498:A
4heoA00 v4_2_0 1.10.8.10.4 2-56:A
4wp2F00 v4_2_0 1.10.8.10.5 600-655:F
4pc3C01 putative 1.10.8.10.6 2-54:C
```

该文件分为四列，每一列分别表示：

* 前四位是PDB的ID，第五位是Chain ID，最后两位是domain number，如果只有一个domain，那就是00

* 表示定义该domain的版本

* CATH分类的代码，分别表示：

  * **C**：Class
  * **A**：Architecture
  * **T**：Topology
  * **H**：Homologous Superfamily
  * **S**：Sequence Family (S35)

  可以参考[这里](http://www.cathdb.info/wiki/doku/?id=faq)查看更多的信息。

* 表示氨基酸序列的范围，可以对应到PDB文件中

AlphaFold的[Github](https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13)目录下有`train_domains.txt`和`test_domains.txt`两个文件，可以直接使用其中的ID。

### 2. 解析MSA特征

> For each training sequence, we searched for and aligned to the training sequence similar protein sequences in the <mark>Uniclust30 dataset</mark> with <mark>HHblits</mark> and used the returned MSA to generate profile features with the position-specific substitution probabilities for each residue as well as covariation features—the parameters of a regularized pseudolikelihood-trained Potts model similar to CCMpred. CCMpred uses the Frobenius norm of the parameters, but we feed both this <mark>norm (1 feature) </mark>and the <mark>raw parameters (484 features)</mark> into the network for each residue pair *ij*. In addition, we provide the network with features that explicitly represent gaps and deletions in the MSA. To make the network better able to make predictions for shallow MSAs, and as a form of data augmentation, we take a sample of <mark>half the sequences</mark> from the the HHblits MSA before computing the MSA-based features. Our training set contains <mark>10 such samples for each domain</mark>. We extract additional profile features using PSI-BLAST.

#### 1. HHblits搜索同源序列

Uniclust30 dataset可以从[这里](http://gwdu111.gwdg.de/~compbiol/uniclust/2020_06/)下载，HHblits可以从[这里](https://github.com/soedinglab/hh-suite)下载，使用下面的命令来搜索：

```bash
hhblits \
    -i input.fasta \
    -d /150T/zhangqf/GenomeAnnotation/UniClust30/UniRef30_2020_06 \
    -o hhblits.hhr \
    -ohhm hhblits.hmm \
    -scores hhblits.score \
    -opsi hhblits.psi \
    -oa3m hhblits.a3m \
    -n 3 \
    -e 0.001 \
    -cpu 20
```

`.hmm`是显著比对的MSA的HHM信息；`.score`是一个矩阵，表示每一个匹配的各项分数；`.hhr`是输出的匹配的结果

#### 2. PSI-Blast搜索同源序列

首先需要准备blast nr索引文件：

```shell
lftp ftp.ncbi.nlm.nih.gov:/blast/db
mget nr.* # 下载解压以后就可以了，不需要合并
```

然后运行psiblast搜索同源序列，这一步应该非常慢：

```shell
# hhblits在生成psi文件时，把序列的ID进行了压缩，所以有些ID是一样的，会导致psiblast读取时报错，所以要先把ID做独一话处理，并转成fas格式文件
awk '{print ">"$1"_"NR"\n"$2}' hhblits.psi > hhblits.fas

# 运行psiblast
psiblast -in_msa hhblits.fas \
	-out_ascii_pssm psiblast.mtx \
	-out_pssm psiblast.chk \
	-out psiblast.out \
	-db /150T/zhangqf/GenomeAnnotation/INDEX/blast/nr/nr \
	-evalue 0.001 \
	-num_iterations 3 \
	-num_threads 20
```

`psiblast.chk`是psiblast输出的checkpoint文件，可以作为[ChkParse](https://github.com/psipred/psipred/blob/master/src/chkparse.c)的输入。这里的`psiblast.mtx`文件就是输入：

```bash
chkparse psiblast.chk > chkparse.out
```

#### 3. Potts model parameters

CCMpred可以从[这里](https://github.com/soedinglab/CCMpred)下载。

```shell
# 转成ccmpred可以接受的输入格式，【注意这个reformat.pl可能不太准】
# reformat.pl a3m psi output.a3m output.psi
# 获得Frobenius norm of 484 channels of pseudolikelihood
awk '{print $2}' hhblits.psi > ccmpred.in
ccmpred -r ccmpred.raw -n 500 ccmpred.in ccmpred.out -t 20
```

参数`-r`表示原始的prediction matrix，`-n`表示最大迭代数，`output.mat`是一个484*484的矩阵。

#### 4. profile_with_prior

> A profile computed using HHBlits which takes into account priors and Blosum matrix. See equation 5 in https://doi.org/10.1093/nar/25.17.3389.

添加`-nocontxt`参数可以使用替换矩阵，而不是上下文特异的伪计数。

`-nocontxt   use substitution-matrix instead of context-specific pseudocounts `

```shell
hhblits \
    -i input.fasta \
    -d /150T/zhangqf/GenomeAnnotation/UniClust30/UniRef30_2020_06 \
    -o hhblits_prior.hhr \
    -ohhm hhblits_prior.hmm \
    -scores hhblits_prior.score \
    -opsi hhblits_prior.psi \
    -oa3m hhblits_prior.a3m \
    -n 3 \
    -e 0.001 \
    -nocontxt \
    -cpu 20
```

### 3. 输入数据的维度

> * Number of HHblits alignments (scalar).
>
> * Sequence-length features: 1-hot amino acid type (21 features); profiles: PSI-BLAST (21 features), HHblits profile (22 features), non-gapped profile (21 features), HHblits bias, HMM profile (30 features), Potts model bias (22 features); deletion probability (1 fea- ture); residue index (integer index of residue number, consecutive except for multi-segment domains, encoded as 5 least-significant bits and a scalar).
>
> * Sequence-length-squared features: Potts model parameters (484 features, fitted with 500 iterations of gradient descent using Nesterov momentum 0.99, without sequence reweighting); Frobenius norm (1 feature); gap matrix (1 feature).

为了显示这些数据的输入维度，我们下载[alphafold](https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13)代码，放在一个叫做`alphafold_casp13`的目录下，然后下载输入文件：http://bit.ly/alphafold-casp13-data (43.5 GB)，取其中一个例子`T0953s1.tfrec`：

```python
import tensorflow as tf
from alphafold_casp13 import contacts_dataset

features = ['num_alignments','aatype','profile','hhblits_profile','non_gapped_profile',
            'hmm_profile','pseudo_bias','deletion_probability','residue_index',
            'pseudolikelihood','pseudo_frob','gap_matrix']

dataset = contacts_dataset.create_tf_dataset('T0953s1.tfrec',tuple(features))
dataset = dataset.batch(1)
iterator = tf.compat.v1.data.make_one_shot_iterator(dataset)
_input_batch = iterator.get_next()
for key in _input_batch:
  print(f"{key}: {_input_batch[key].shape}")
  print(_input_batch[key])
```

##### 一维数据

> Number of HHblits alignments (scalar)

`num_alignments`: (1, 64, 1)
    整数：这个例子都是120

##### 二维数据

> 1-hot amino acid type (21 features)

`aatype`: (1, 64, 21)
    0和1，one-hot编码

> PSI-BLAST profiles (21 features)

`profile`: (1, 64, 21)
    浮点数，(-15.1, 13.5)

> HHblits profiles (22 features)

`hhblits_profile`: (1, 64, 22)
    0和大于0的浮点数, (0.0, 1.0)

> non-gapped profile (21 features)

`non_gapped_profile`: (1, 64, 21)
    0和大于0的浮点数, (0.0, 1.0)

> HHblits bias

这个不存在，应该是文章写错了~

> HMM profile (30 features)

`hmm_profile`: (1, 64, 30)
    0和大于0的浮点数, (0.0, 1.0)

> Potts model bias (22 features)

`pseudo_bias`: (1, 64, 22)
    浮点数，均值为0，但是带有null

> deletion probability (1 feature)

`deletion_probability`: (1, 64, 1)
    大于0的浮点数, (0, 0.998)

> residue index (integer index of residue number, consecutive except for multi-segment domains, encoded as 5 least-significant bits and a scalar).

`residue_index`：(1, 64, 1)

元素是0-63

##### 三维数据

> Potts model parameters (484 features, fitted with 500 iterations of gradient descent using Nesterov momentum 0.99, without sequence reweighting)

`pseudolikelihood`: (1, 64, 64, 484)
    浮点数，均值为0，但是带有null

> Frobenius norm (1 feature)

`pseudo_frob`: (1, 64, 64, 1)
    浮点数, (-3.28, 2.25)

> gap matrix (1 feature)

`gap_matrix`: (1, 64, 64, 1)
    大于0的浮点数, (0, 0.999)

### 4. 准备输入数据

下面以`T0953s1`为例，来准备它的输入数据。把它命名为`input.fasta`，使用上面的`hhblits`和`psiblast`命令寻找多序列比对。

```
>T0953s1 (72 residues)
GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS
```

#### aatype

对于氨基酸序列，就是一个one-hot编码，注意`X`表示任意氨基酸。

```python
def sequence_to_onehot(sequence):
  """Maps the given sequence into a one-hot encoded matrix."""
  mapping = {aa: i for i, aa in enumerate('ARNDCQEGHILKMFPSTWYVX')}
  num_entries = max(mapping.values()) + 1
  one_hot_arr = np.zeros((len(sequence), num_entries), dtype=np.int32)
  
  for aa_index, aa_type in enumerate(sequence):
    aa_id = mapping[aa_type]
    one_hot_arr[aa_index, aa_id] = 1
  
  return one_hot_arr

seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
sequence_to_onehot(seq).shape # (72, 21)
```

反向把one-hot转成氨基酸序列：

```python
aaall = 'ARNDCQEGHILKMFPSTWYVX'
aaa = tf.argmax(_input_batch['aatype'][0],1).numpy()
print( "".join([ aaall[i] for i in aaa ]) )
# GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAI
```

#### profile

psi-blast的profile。有两种产生方式，一种是在`psiblast`命令中加`-out_ascii_pssm`参数直接生成mtx文件；另外一种是加`-out_pssm`参数生成checkpoint文件，然后使用[ChkParse](https://github.com/psipred/psipred/blob/master/src/chkparse.c)处理checkpoint文件。

```shell
chkparse psiblast.chk > chkparse.out
```

```python
aaall = 'ARNDCQEGHILKMFPSTWYVX'
ncbicodes = "*A*CDEFGHIKLMNPQRSTVWXY*****"

profile = []
for line in open("chkparse.out"):
    if line.startswith('-32768'):
        data = line.strip().split()
        mtx = [ int(data[ncbicodes.find(alp)]) for alp in aaall ]
        profile.append( mtx )

profile = np.array(profile, dtype='float32')
profile = profile / profile.std() # 这里做一个归一化，但是不清楚AlphaFold是不是这么做的
profile = np.expand_dims(profile,0)
print(profile.shape) # (1, 72, 21)
```

#### hhblits_profile

比对谱中22种残基的分布：20种残基+`X`+`-`

```python
qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
aaall = 'ARNDCQEGHILKMFPSTWYVX-'
hhblits_profile = np.zeros( [len(qeury_seq), len(aaall)] )

msa = [line.strip().split()[1] for line in open('hhblits.psi')]
for seq in msa:
    for col, alpha in enumerate(seq):
        hhblits_profile[col,aaall.find(alpha)] += 1

hhblits_profile /= len(msa)
hhblits_profile = np.expand_dims(hhblits_profile,0)
print(hhblits_profile.shape) # (1, 72, 22)
```

#### deletion_probability

表示缺失状态的序列比例。这里读入`.a3m`文件，`a3m`文件中的每条序列都和第一条进行比较，小写字母表示一个gap。

```python
hhblits_a3m_sequences = []
for line in open('hhblits.a3m'):
    if line[0] != '>':
        hhblits_a3m_sequences.append( line.strip() )

deletion_matrix = []
for msa_sequence in hhblits_a3m_sequences:
  deletion_vec = []
  deletion_count = 0
  for j in msa_sequence:
    if j.islower():
      deletion_count += 1
    else:
      deletion_vec.append(deletion_count)
      deletion_count = 0
  deletion_matrix.append(deletion_vec)

deletion_matrix = np.array(deletion_matrix)
deletion_matrix[deletion_matrix != 0] = 1.0
deletion_probability = deletion_matrix.sum(axis=0) / len(deletion_matrix)
deletion_probability = deletion_probability[np.newaxis, :, np.newaxis]
print(deletion_probability.shape) # (1, 72, 1)
```

#### gap_matrix

gap共发生的概率分布，是一个$L \times L$的矩阵$M$，$L$表示序列的长度，每一个元素$M[i,j]$表示比对中第$i$位和第$j$位共同是gap的序列数量。

```python
MSA = []
for line in open('hhblits.psi'):
    MSA.append( line.strip().split()[1] )

gap_count = [ [ s=="-" for s in seq ] for seq in MSA ]
gap_count = np.array(gap_count, dtype='float32')
gap_matrix = np.matmul(gap_count.T, gap_count)
gap_matrix /= len(MSA)
gap_matrix = gap_matrix[ np.newaxis, :, :, np.newaxis ]
print( gap_matrix.shape ) # (1, 72, 72, 1)
```

#### hmm_profile

```python
def extract_hmm_profile(hhm_file, sequence, asterisks_replace=0.0):
  """Extracts information from the hmm file and replaces asterisks."""
  profile_part = hhm_file.split('#')[-1]
  profile_part = profile_part.split('\n')
  whole_profile = [i.split() for i in profile_part]
  # This part strips away the header and the footer.
  whole_profile = whole_profile[5:-2]
  gap_profile = np.zeros((len(sequence), 10))
  aa_profile = np.zeros((len(sequence), 20))
  count_aa = 0
  count_gap = 0
  for line_values in whole_profile:
    if len(line_values) == 23:
      # The first and the last values in line_values are metadata, skip them.
      for j, t in enumerate(line_values[2:-1]):
        aa_profile[count_aa, j] = (
            2**(-float(t) / 1000.) if t != '*' else asterisks_replace)
      count_aa += 1
    elif len(line_values) == 10:
      for j, t in enumerate(line_values):
        gap_profile[count_gap, j] = (
            2**(-float(t) / 1000.) if t != '*' else asterisks_replace)
      count_gap += 1
    elif not line_values:
      pass
    else:
      raise ValueError('Wrong length of line %s hhm file. Expected 0, 10 or 23'
                       'got %d'%(line_values, len(line_values)))
  hmm_profile = np.hstack([aa_profile, gap_profile])
  assert len(hmm_profile) == len(sequence)
  return hmm_profile

qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
hhm_file = open('hhblits.hmm').read()
hmm_profile = extract_hmm_profile( hhm_file,  qeury_seq)
hmm_profile = hmm_profile[np.newaxis, :]
print(hmm_profile.shape) # (1, 72, 30)
```

#### non_gapped_profile

不考虑gap。

```python
qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
aaall = 'ARNDCQEGHILKMFPSTWYVX'
non_gapped_profile = np.zeros( [len(qeury_seq), len(aaall)] )

msa = [line.strip().split()[1] for line in open('hhblits.psi')]
for seq in msa:
    for col, alpha in enumerate(seq):
        if alpha != '-':
            non_gapped_profile[col,aaall.find(alpha)] += 1

non_gapped_profile /= len(msa)
non_gapped_profile = np.expand_dims(non_gapped_profile,0)
print(non_gapped_profile.shape) # (1, 72, 21)
```

#### reweighted_profile

对每条序列分配不同的权重，罕见的序列权重更高。首先计算序列权重向量：

```python
def sequence_weights(sequence_matrix):
  """Compute sequence reweighting to weight rarer sequences higher."""
  num_rows, num_res = sequence_matrix.shape
  cutoff = 0.62 * num_res
  weights = np.ones(num_rows, dtype=np.float32)
  for i in range(num_rows):
    for j in range(i + 1, num_rows):
      similarity = (sequence_matrix[i] == sequence_matrix[j]).sum()
      if similarity > cutoff:
        weights[i] += 1
        weights[j] += 1
  return 1.0 / weights

msa = np.array([list(line.strip().split()[1]) for line in open('hhblits.psi')])
seq_weight = sequence_weights(msa)
print(seq_weight.shape) # (230,)
```

然后重新计算`reweighted_profile`：

```python
qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
aaall = 'ARNDCQEGHILKMFPSTWYVX-'
reweighted_profile = np.zeros( [len(qeury_seq), len(aaall)] )

msa = [line.strip().split()[1] for line in open('hhblits.psi')]
for i,seq in enumerate(msa):
    for col, alpha in enumerate(seq):
        reweighted_profile[col,aaall.find(alpha)] += seq_weight[i]

reweighted_profile /= len(msa)
reweighted_profile = np.expand_dims(reweighted_profile,0)
print(reweighted_profile.shape) # (1, 72, 22)
```

#### profile_with_prior

`HHBlits`计算的prior，但是会考虑priors and Blosum matrix。

```python
qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'
aaall = 'ARNDCQEGHILKMFPSTWYVX-'
profile_with_prior = np.zeros( [len(qeury_seq), len(aaall)] )

msa = [line.strip().split()[1] for line in open('hhblits_prior.psi')]
for seq in msa:
    for col, alpha in enumerate(seq):
        profile_with_prior[col,aaall.find(alpha)] += 1

profile_with_prior /= len(msa)
profile_with_prior = np.expand_dims(profile_with_prior,0)
print(profile_with_prior.shape) # (1, 72, 22)
```

#### pseudo_bias/pseudo_frob/pseudolikelihood

计算pseudo-likelihood时，作者应该没有使用ccmpred，因为CCMPred计算的频率和作者的维度不一样。

```python
qeury_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'

pseudo_bias = []
pseudolikelihood = np.zeros([len(qeury_seq), len(qeury_seq), 21, 21])

IN = open('ccmpred.raw')
pseudolikelihood_pair = []
is_bias = True
left, right = 0, 0
for line in IN:
    if line[0] != '#':
        if is_bias:
            pseudo_bias.append([float(d) for d in line.strip().split()])
        else:
            pseudolikelihood_pair.append([float(d) for d in line.strip().split()])
    else:
        is_bias = False
        if pseudolikelihood_pair:
            pseudolikelihood[left, right] = pseudolikelihood_pair
            pseudolikelihood_pair = []
        left, right = line.strip().split()[1:]
        left, right = int(left), int(right)

pseudo_bias = np.array(pseudo_bias)
pseudo_bias = pseudo_bias[np.newaxis, :, :]
print(pseudo_bias.shape) # (1, 72, 20)

pseudolikelihood = np.reshape(pseudolikelihood, [1, len(qeury_seq), len(qeury_seq), 21*21])
print(pseudolikelihood.shape) # (1, 72, 72, 441)

pseudo_frob = np.array([ [ float(d) for d in line.strip().split() ] for line in open('ccmpred.out') ])
pseudo_frob = pseudo_frob[np.newaxis, :, :, np.newaxis]
print(pseudo_frob.shape) # (1, 72, 72, 1)
```

