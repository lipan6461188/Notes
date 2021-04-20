<center><h1>准备AlphaFold1原始TF代码的输入</h1></center>

### 1. 测试tf.data.Dataset

> 由于AlphaFold1默认把一个`.tfrec`文件读取成tf.data.Dataset对象，首先来研究一下该文件该如何读取

导入TensorFlow，如果是2.x版本，就导入版本1：

```python
# 版本1.x
import tensorflow as tf
# 版本2.x
import tensorflow.compat.v1 as tf
```

定义其他变量：

```python
NUM_RES = 'num residues placeholder'
FEATURES = {
    'aatype': (tf.float32, [NUM_RES, 21]),
    'alpha_mask': (tf.int64, [NUM_RES, 1]),
    'alpha_positions': (tf.float32, [NUM_RES, 3]),
    'beta_mask': (tf.int64, [NUM_RES, 1]),
    'beta_positions': (tf.float32, [NUM_RES, 3]),
    'between_segment_residues': (tf.int64, [NUM_RES, 1]),
    'chain_name': (tf.string, [1]),
    'deletion_probability': (tf.float32, [NUM_RES, 1]),
    'domain_name': (tf.string, [1]),
    'gap_matrix': (tf.float32, [NUM_RES, NUM_RES, 1]),
    'hhblits_profile': (tf.float32, [NUM_RES, 22]),
    'hmm_profile': (tf.float32, [NUM_RES, 30]),
    'key': (tf.string, [1]),
    'mutual_information': (tf.float32, [NUM_RES, NUM_RES, 1]),
    'non_gapped_profile': (tf.float32, [NUM_RES, 21]),
    'num_alignments': (tf.int64, [NUM_RES, 1]),
    'num_effective_alignments': (tf.float32, [1]),
    'phi_angles': (tf.float32, [NUM_RES, 1]),
    'phi_mask': (tf.int64, [NUM_RES, 1]),
    'profile': (tf.float32, [NUM_RES, 21]),
    'profile_with_prior': (tf.float32, [NUM_RES, 22]),
    'profile_with_prior_without_gaps': (tf.float32, [NUM_RES, 21]),
    'pseudo_bias': (tf.float32, [NUM_RES, 22]),
    'pseudo_frob': (tf.float32, [NUM_RES, NUM_RES, 1]),
    'pseudolikelihood': (tf.float32, [NUM_RES, NUM_RES, 484]),
    'psi_angles': (tf.float32, [NUM_RES, 1]),
    'psi_mask': (tf.int64, [NUM_RES, 1]),
    'residue_index': (tf.int64, [NUM_RES, 1]),
    'resolution': (tf.float32, [1]),
    'reweighted_profile': (tf.float32, [NUM_RES, 22]),
    'sec_structure': (tf.int64, [NUM_RES, 8]),
    'sec_structure_mask': (tf.int64, [NUM_RES, 1]),
    'seq_length': (tf.int64, [NUM_RES, 1]),
    'sequence': (tf.string, [1]),
    'solv_surf': (tf.float32, [NUM_RES, 1]),
    'solv_surf_mask': (tf.int64, [NUM_RES, 1]),
    'superfamily': (tf.string, [1]),
}
i_features = [
      "profile",
      "hhblits_profile",
      "aatype",
      "pseudo_frob",
      "pseudolikelihood",
      "deletion_probability",
      "gap_matrix",
      "pseudo_bias",
      "profile_with_prior",
      "profile_with_prior_without_gaps",
      "reweighted_profile",
      "non_gapped_profile",
      "hmm_profile",
      "num_alignments",
      "seq_length"]
i_scalars = [ "num_effective_alignments", "resolution"]
i_target = [
      "sec_structure",
      "sec_structure_mask",
      "solv_surf",
      "solv_surf_mask",
      "beta_positions",
      "beta_mask",
      "domain_name",
      "chain_name",
      "resolution",
      "num_alignments",
      "superfamily",
      "profile",
      "hhblits_profile",
      "residue_index",
      "between_segment_residues"]
```

读取输入数据集：

```python
tf_dataset = tf.data.TFRecordDataset(filenames=['T0949.tfrec']) # 返回一个Dataset类型
```

#### 1.1 不解析直接查看数据

```python
tf_dataset = tf_dataset.batch(1) # 由于数据的维度不一样，所以一次只能取出一个example
iterator = tf.compat.v1.data.make_one_shot_iterator(tf_dataset)

_input_batch = iterator.get_next()
print(type(_input_batch)) 
# <class 'tensorflow.python.framework.ops.EagerTensor'>

# 解析这个特征
parsed_features = tf.io.parse_single_example(_input_batch[0], 
        { 'aatype': tf.io.FixedLenSequenceFeature(shape=(), 
                   dtype=FEATURES['aatype'][0], allow_missing=True) }
    )
print( parsed_features['aatype'].shape ) # (2688,)

# 解析的同时指定shape
parsed_features = tf.io.parse_single_example(_input_batch[0], 
        { 'aatype': tf.io.FixedLenSequenceFeature(shape=[128, 21], dtype=FEATURES['aatype'][0], allow_missing=True) }   # 这里注意，序列的长度是128，所以是128
)
print( parsed_features['aatype'].shape ) # (1, 128, 21)
```

#### 1.2 使用map函数进行解析

```python
features = i_features + i_scalars + i_target
required_features = ['aatype', 'sequence', 'seq_length']
features = list(set(required_features) | set(features))
features = {name: FEATURES[name] for name in features}

def shape(feature_name, num_residues):
    # 数据类型和shape
    unused_dtype, raw_sizes = FEATURES[feature_name]
    replacements = {NUM_RES: num_residues}
    sizes = [replacements.get(dimension, dimension) for dimension in raw_sizes]
    return sizes

def parse_tfexample(raw_data, features):
    feature_map = {
      k: tf.io.FixedLenSequenceFeature(shape=(), dtype=v[0], allow_missing=True) for k, v in features.items()
    }
    # 把数据按照格式解析
    parsed_features = tf.io.parse_single_example(raw_data, feature_map) 
   	# 每一个patch都有一个seq_length
    num_residues = tf.cast(parsed_features['seq_length'][0], dtype=tf.int32) 
    for k, v in parsed_features.items():
        new_shape = shape(feature_name=k, num_residues=num_residues)
        parsed_features[k] = tf.reshape(v, new_shape, name='reshape_%s' % k)
    return parsed_features

# 这里是分别一个一个地处理实例，每一个实例可能是那个蛋白的一个patch
dataset = tf_dataset.map(lambda raw: parse_tfexample(raw, features)) 
dataset = dataset.batch(1)
iterator = tf.compat.v1.data.make_one_shot_iterator(dataset)
_input_batch = iterator.get_next()

print( _input_batch['aatype'].shape ) # (1, 128, 21)

for f in features:
    print(f, type(_input_batch[f]), _input_batch[f].shape)

data_to_save = { f:_input_batch[f].numpy() for f in _input_batch }
```

### 2. 把读取的`.tfrec`写入新的`.tfrec`

首先定义几个解析函数：

```python
def _bytes_feature(value):
  """Returns a bytes_list from a string / byte."""
  if isinstance(value, type(tf.constant(0))):
    value = value.numpy() # BytesList won't unpack a string from an EagerTensor.
  return tf.train.Feature(bytes_list=tf.train.BytesList(value=value))

def _float_feature(value):
  """Returns a float_list from a float / double."""
  return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def _int64_feature(value):
  """Returns an int64_list from a bool / enum / int / uint."""
  return tf.train.Feature(int64_list=tf.train.Int64List(value=value))
```

然后一边读取，一边写入：

```python
writer = tf.python_io.TFRecordWriter('test.tfrec') # 打开一个文件
iterator = tf.compat.v1.data.make_one_shot_iterator(dataset)
for _input_batch in iterator: # 每次写入一个example
    data_to_save = { f:_input_batch[f].numpy() for f in _input_batch }
    feature_parse = {}
    for k in data_to_save:
        if FEATURES[k][0] == tf.string:
            #print(data_to_save[k])
            feature_parse[k] = _bytes_feature( data_to_save[k][0] )
        elif FEATURES[k][0] == tf.float32:
            feature_parse[k] = _float_feature( np.reshape(data_to_save[k],[-1]) )
        else:
            feature_parse[k] = _int64_feature( np.reshape(data_to_save[k],[-1]) )
    example = tf.train.Example(features=tf.train.Features( feature=feature_parse ))
    writer.write(example.SerializeToString())

writer.close()
```








