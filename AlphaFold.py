
import sys, os
import numpy as np
import collections
import Cluster
import General
import Colors

AATYPES = 'ARNDCQEGHILKMFPSTWYVX'

# Database index directories
UniClust30_db = '/150T/zhangqf/GenomeAnnotation/UniClust30/UniRef30_2020_06'
Blast_nr_db = '/150T/zhangqf/GenomeAnnotation/INDEX/blast/nr/nr'

plmDCAPATH = '/BioII/zhangqf7/lipan/CPI/plmDCA-master/plmDCA_asymmetric_v2'

threeToOne = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "ASX": "B",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLX": "Z",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}
oneToThree = { b:a for a,b in threeToOne.items() }

# [1, seq_len, 1]
def get_num_alignments(hhblits_psi_file, query_seq):
    """
    All values are align_len
    Return [1, seq_len, 1]
    """
    x = np.zeros([1, len(query_seq), 1])
    x[:] = len(open(hhblits_psi_file).readlines())
    return x

# [1, seq_len, 21]
def get_aatype(query_seq):
    """
    One-hot encode protein sequences
    Return [1, seq_len, 21]
    """
    mapping = {aa: i for i, aa in enumerate(AATYPES)}
    num_entries = max(mapping.values()) + 1
    one_hot_arr = np.zeros((len(query_seq), 21), dtype=np.int32)
    
    for aa_index, aa_type in enumerate(query_seq):
        aa_id = mapping[aa_type]
        one_hot_arr[aa_index, aa_id] = 1
    
    one_hot_arr = one_hot_arr[np.newaxis, :, :]
    return one_hot_arr

# [1, seq_len, 21]
def get_profile(chkparse_file):
    """
    PSI-BLAST profile
    Return [1, seq_len, 21]
    """
    ncbicodes = "*A*CDEFGHIKLMNPQRSTVWXY*****"
    ref_seq = "ACDEFGHIKLMNPQRSTVWXY"
    
    profile = []
    for line in open(chkparse_file):
        if line.startswith('-32768'):
            data = line.strip().split()
            mtx = [ int(data[ncbicodes.find(aa)]) for aa in ref_seq ]
            profile.append( mtx )
    
    profile = np.array(profile, dtype='float32')
    #profile = profile / profile.mean() # 这里做一个归一化，但是不清楚AlphaFold是不是这么做的
    profile = profile[np.newaxis, :, :]
    
    return profile

# [1, seq_len, 22]
def get_hhblits_profile(hhblits_psi_file, query_seq):
    """
    HHBlits profile
    Return [1, seq_len, 22]
    """
    _AATYPES = AATYPES+'-'
    hhblits_profile = np.zeros( [len(query_seq), len(_AATYPES)] )
    
    ref_seq_order = "ACDEFGHIKLMNPQRSTVWXY-"
    msa = [line.strip().split()[1] for line in open(hhblits_psi_file)]
    for seq in msa:
        for col, aa in enumerate(seq):
            hhblits_profile[col,ref_seq_order.find(aa)] += 1
    
    hhblits_profile /= len(msa)
    hhblits_profile = hhblits_profile[np.newaxis, :, :]
    
    return hhblits_profile

# [1, seq_len, 21]
def get_non_gapped_profile(hhblits_psi_file, query_seq):
    """
    HHBlits profile no gap
    Return [1, seq_len, 21]
    """
    ref_seq_order = "ACDEFGHIKLMNPQRSTVWXY"
    non_gapped_profile = np.zeros( [len(query_seq), len(ref_seq_order)] )
    
    msa = [line.strip().split()[1] for line in open(hhblits_psi_file)]
    for seq in msa:
        for col, alpha in enumerate(seq):
            if alpha != '-':
                non_gapped_profile[col,ref_seq_order.find(alpha)] += 1
    
    non_gapped_profile /= len(msa)
    non_gapped_profile = non_gapped_profile[np.newaxis, :, :]
    
    return non_gapped_profile

# [1, seq_len, 30]
def get_hmm_profile(hmm_file, query_seq, asterisks_replace=0.0):
    """
    Calculate the HMM profile
    Return [1, seq_len, 30]
    """
    profile_part = open(hmm_file).read().split('#')[-1]
    profile_part = profile_part.split('\n')
    whole_profile = [i.split() for i in profile_part]
    whole_profile = whole_profile[5:-2]
    gap_profile = np.zeros((len(query_seq), 10))
    aa_profile = np.zeros((len(query_seq), 20))
    count_aa = 0
    count_gap = 0
    for line_values in whole_profile:
        if len(line_values) == 23:
            # The first and the last values in line_values are metadata, skip them.
            for j, t in enumerate(line_values[2:-1]):
                aa_profile[count_aa, j] = (2**(-float(t) / 1000.) if t != '*' else asterisks_replace)
            count_aa += 1
        elif len(line_values) == 10:
            for j, t in enumerate(line_values):
                gap_profile[count_gap, j] = (2**(-float(t) / 1000.) if t != '*' else asterisks_replace)
            count_gap += 1
        elif not line_values:
            pass
        else:
            raise ValueError('Wrong length of line %s hhm file. Expected 0, 10 or 23 got %d'%(line_values, len(line_values)))
    
    hmm_profile = np.hstack([aa_profile, gap_profile])
    assert len(hmm_profile) == len(query_seq)
    
    hmm_profile = hmm_profile[np.newaxis, :]
    return hmm_profile

# [1, seq_len, seq_len, 441], [1, seq_len, 20]
def get_pseudolikelihood_ccmpred(ccmpred_raw_file, query_seq):
    """
    Get pseudolikelihood and pseudolikelihood bias
    Return [1, seq_len, seq_len, 441], [1, seq_len, 20]
    """
    pseudolikelihood = np.zeros([len(query_seq), len(query_seq), 21, 21])
    pseudo_bias = []
    
    pseudolikelihood_pair = []
    is_bias = True
    left, right = 0, 0
    compressed = False
    if ccmpred_raw_file.endswith('.gz'):
        import gzip
        IN = gzip.open(ccmpred_raw_file)
        compressed = True
    else:
        IN = open(ccmpred_raw_file)
    for line in IN:
        if compressed:
            line = line.decode()
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
    
    pseudolikelihood = np.reshape(pseudolikelihood, [1, len(query_seq), len(query_seq), 21*21])
    pseudo_bias = np.array(pseudo_bias)
    pseudo_bias = pseudo_bias[np.newaxis, :, :]
    
    return pseudolikelihood, pseudo_bias

# [1, seq_len, seq_len, 1]
def get_pseudolikelihood_frob_ccmpred(ccmpred_out_file):
    """
    Get pseudolikelihood frob
    Return [1, seq_len, seq_len, 1]
    """
    pseudo_frob = np.array([ [ float(d) for d in line.strip().split() ] for line in open(ccmpred_out_file) ])
    pseudo_frob = pseudo_frob[np.newaxis, :, :, np.newaxis]
    return pseudo_frob

# [1, seq_len, seq_len, 484], [1, seq_len, 22], [1, seq_len, seq_len, 1]
def get_pseudolikelihood(matlab_file):
    """
    The matlab_file is generated by plmDCA
    Get pseudolikelihood, pseudolikelihood bias, pseudolikelihood_frob
    Return [1, seq_len, seq_len, 484], [1, seq_len, 22], [1, seq_len, seq_len, 1]
    """
    from scipy.io import loadmat
    potts = loadmat(matlab_file)
    seq_len = potts['pseudo_bias'].shape[0]
    
    #index = np.array(list(range(1, 22)) + [0])
    #potts['pseudo_bias'] = potts['pseudo_bias'][:,index]
    #potts['pseudolikelihood'] = np.reshape(potts['pseudolikelihood'], [seq_len, seq_len, 22, 22])[:, :, :, index][:,:,index,:]
    #potts['pseudolikelihood'] = np.reshape(potts['pseudolikelihood'], [seq_len, seq_len, 484])
    
    pseudolikelihood = potts['pseudolikelihood'][np.newaxis, :, :, :]
    pseudo_bias = potts['pseudo_bias'][np.newaxis, :, :]
    pseudo_frob = potts['pseudo_frob'][np.newaxis, :, :, np.newaxis]
    
    return pseudolikelihood, pseudo_bias, pseudo_frob

# [1, seq_len, 1]
def get_deletion_probability(hhblits_a3m_file):
    """
    Calculate the insertion rate for each base
    Return [1, seq_len, 1]
    """
    hhblits_a3m_sequences = []
    for line in open(hhblits_a3m_file):
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
    
    return deletion_probability

# [1, seq_len, 1]
def get_residue_index(query_seq):
    """
    Get residual index array
    Return [1, seq_len, 1]
    """
    return np.reshape( range(len(query_seq)), [1, -1, 1] )

# [1, seq_len, seq_len, 1]
def get_gap_matrix(hhblits_psi_file):
    """
    Get gap matrix, M[i,j] count how many sequences with i,j are gaps
    Return [1, seq_len, seq_len, 1]
    """
    MSA = []
    for line in open(hhblits_psi_file):
        MSA.append( line.strip().split()[1] )
    
    gap_count = [ [ s=="-" for s in seq ] for seq in MSA ]
    gap_count = np.array(gap_count, dtype='float32')
    gap_matrix = np.matmul(gap_count.T, gap_count)
    gap_matrix /= len(MSA)
    gap_matrix = gap_matrix[ np.newaxis, :, :, np.newaxis ]
    return gap_matrix

# [1, seq_len, 22]
def get_reweighted_profile(hhblits_psi_file, query_seq):
    """
    Assign different weight for each sequence, higher weight for rarer sequence
    Return: [1, seq_len, 22]
    """
    ref_seq = "ACDEFGHIKLMNPQRSTVWXY-"
    msa = np.array([list(line.strip().split()[1]) for line in open(hhblits_psi_file)])
    num_rows, num_res = msa.shape
    cutoff = 0.62 * num_res
    weights = np.ones(num_rows, dtype=np.float32)
    for i in range(num_rows):
        for j in range(i + 1, num_rows):
            similarity = (msa[i] == msa[j]).sum()
            if similarity > cutoff:
                weights[i] += 1
                weights[j] += 1
    weights = 1.0 / weights
    
    #AATYPES_ = AATYPES + '-'
    reweighted_profile = np.zeros( [len(query_seq), len(ref_seq)] )
    
    #msa = [line.strip().split()[1] for line in open(hhblits_psi_file)]
    for i,seq in enumerate(msa):
        for j,aa in enumerate(seq):
            reweighted_profile[j,ref_seq.find(aa)] += weights[i]
    
    reweighted_profile /= len(msa)
    reweighted_profile = reweighted_profile[np.newaxis, :, :]
    
    return reweighted_profile

Input = collections.namedtuple(
    'Input', [
        'num_alignments', # [1, seq_len, 1]
        'aatype', # [1, seq_len, 21]
        'profile', # [1, seq_len, 21]
        'hhblits_profile', # [1, seq_len, 22]
        'non_gapped_profile', # [1, seq_len, 21]
        'hhblits_prior_profile', # [1, seq_len, 22]
        'non_gapped_prior_profile', # [1, seq_len, 21]
        'hmm_profile', # [1, seq_len, 30] 
        'pseudolikelihood', # [1, seq_len, seq_len, 484]
        'pseudo_bias', # [1, seq_len, 22]
        'pseudolikelihood_frob', # [1, seq_len, seq_len, 1]
        'deletion_probability', # [1, seq_len, 1]
        'residue_index', # [1, seq_len, 1]
        'gap_matrix', # [1, seq_len, seq_len, 1]
        'reweighted_profile' # [1, seq_len, 22]
    ])


####################################
### Functions for batch operations
####################################

def read_data_array_from_file(hhblits_psi_file, hhblits_prior_psi_file, hmm_file, chkparse_file, matlab_file, hhblits_a3m_file, query_seq):
    """
    First run prepare_alphafold1_input function
    
    Return collections.namedtuple Input object with 15 elements
    """
    num_alignments = get_num_alignments(hhblits_psi_file, query_seq) # [1, seq_len, 1]
    aatype = get_aatype(query_seq) # [1, seq_len, 21]
    profile = get_profile(chkparse_file) # [1, seq_len, 21]
    hhblits_profile = get_hhblits_profile(hhblits_psi_file, query_seq) # [1, seq_len, 22]
    non_gapped_profile = get_non_gapped_profile(hhblits_psi_file, query_seq) # [1, seq_len, 21]
    hhblits_prior_profile = get_hhblits_profile(hhblits_prior_psi_file, query_seq) # [1, seq_len, 22]
    non_gapped_prior_profile = get_non_gapped_profile(hhblits_prior_psi_file, query_seq) # [1, seq_len, 21]
    hmm_profile = get_hmm_profile(hmm_file, query_seq) # [1, seq_len, 30]
    #pseudolikelihood, pseudo_bias = get_pseudolikelihood(ccmpred_raw_file, query_seq) # [1, seq_len, seq_len, 441], [1, seq_len, 20]
    #pseudolikelihood_frob = get_pseudolikelihood_frob(ccmpred_out_file) # [1, seq_len, seq_len, 1]
    pseudolikelihood, pseudo_bias, pseudolikelihood_frob = get_pseudolikelihood(matlab_file) # [1, seq_len, seq_len, 484], [1, seq_len, 22], [1, seq_len, seq_len, 1]
    deletion_probability = get_deletion_probability(hhblits_a3m_file) # [1, seq_len, 1]
    residue_index = get_residue_index(query_seq) # [1, seq_len, 1]
    gap_matrix = get_gap_matrix(hhblits_psi_file) # [1, seq_len, seq_len, 1]
    reweighted_profile = get_reweighted_profile(hhblits_psi_file, query_seq) # [1, seq_len, 22]
    
    return Input(
      num_alignments=num_alignments,
      aatype=aatype,
      profile=profile,
      hhblits_profile=hhblits_profile,
      non_gapped_profile=non_gapped_profile,
      hhblits_prior_profile=hhblits_prior_profile,
      non_gapped_prior_profile=non_gapped_prior_profile,
      hmm_profile=hmm_profile,
      pseudolikelihood=pseudolikelihood,
      pseudo_bias=pseudo_bias,
      pseudolikelihood_frob=pseudolikelihood_frob,
      deletion_probability=deletion_probability,
      residue_index=residue_index,
      gap_matrix=gap_matrix,
      reweighted_profile=reweighted_profile)

get_alphafold1_input_array = read_data_array_from_file

# this function has been deprecated
def _reshape_input_array(instance):
    """
    Input a AlphaFold.Input object
    Reshape all object to [batch_size, channels, height, width]
    Return AlphaFold.Input object
    """
    assert isinstance(instance, Input)
    seq_len = instance[0].shape[1]
    
    # num_alignments: [1, seq_len, 1] -> [1, 1, seq_len, seq_len]
    num_alignments = np.repeat(instance[0][:, np.newaxis, :, :], seq_len, 3)
    # aatype: [1, seq_len, 21] -> [1, 21, seq_len, seq_len]
    aatype = np.repeat(np.transpose(instance[1], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # profile: [1, seq_len, 21] -> [1, 21, seq_len, seq_len]
    profile = np.repeat(np.transpose(instance[2], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # hhblits_profile: [1, seq_len, 22] -> [1, 22, seq_len, seq_len]
    hhblits_profile = np.repeat(np.transpose(instance[3], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # non_gapped_profile: [1, seq_len, 21] -> [1, 21, seq_len, seq_len]
    non_gapped_profile = np.repeat(np.transpose(instance[4], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # hhblits_prior_profile: [1, seq_len, 22] -> [1, 22, seq_len, seq_len]
    hhblits_prior_profile = np.repeat(np.transpose(instance[5], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # non_gapped_prior_profile: [1, seq_len, 21]-> [1, 21, seq_len, seq_len]
    non_gapped_prior_profile = np.repeat(np.transpose(instance[6], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # hmm_profile: [1, seq_len, 30]-> [1, 30, seq_len, seq_len]
    hmm_profile = np.repeat(np.transpose(instance[7], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # pseudolikelihood: [1, seq_len, seq_len, 484] -> [1, 484, seq_len, seq_len]
    pseudolikelihood = np.transpose(instance[8], [0,3,1,2])
    # pseudo_bias: [1, seq_len, 20] -> [1, 20, seq_len, seq_len]
    pseudo_bias = np.repeat(np.transpose(instance[9], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # pseudolikelihood_frob: [1, seq_len, seq_len, 1] -> [1, 1, seq_len, seq_len]
    pseudolikelihood_frob = np.transpose(instance[10], [0,3,1,2])
    # deletion_probability: [1, seq_len, 1] -> [1, 1, seq_len, seq_len]
    deletion_probability = np.repeat(np.transpose(instance[11], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # residue_index: [1, seq_len, 1] -> [1, 1, seq_len, seq_len]
    residue_index = np.repeat(np.transpose(instance[12], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    # gap_matrix: [1, seq_len, seq_len, 1] -> [1, 1, seq_len, seq_len]
    gap_matrix = np.transpose(instance[13], [0,3,1,2])
    # reweighted_profile: [1, seq_len, 22] -> [1, 22, seq_len, seq_len]
    reweighted_profile = np.repeat(np.transpose(instance[14], [0,2,1])[:,:,:,np.newaxis], seq_len, 3)
    
    return Input(
      num_alignments=num_alignments,
      aatype=aatype,
      profile=profile,
      hhblits_profile=hhblits_profile,
      non_gapped_profile=non_gapped_profile,
      hhblits_prior_profile=hhblits_prior_profile,
      non_gapped_prior_profile=non_gapped_prior_profile,
      hmm_profile=hmm_profile,
      pseudolikelihood=pseudolikelihood,
      pseudo_bias=pseudo_bias,
      pseudolikelihood_frob=pseudolikelihood_frob,
      deletion_probability=deletion_probability,
      residue_index=residue_index,
      gap_matrix=gap_matrix,
      reweighted_profile=reweighted_profile)

# prepare_alphafold1_input
def run_alphafold1_input_pipeline(query_fasta, outdir, cpu=20, verbose=True, noExecute=False, use_LSF=True, LSF_parameters={}):
    """
    Given a query fasta file, search homologous sequences with hhblits and psiblast
    
    Parameters:
        query_fasta         -- Input fasta file
        outdir              -- Output directory
        cpu                 -- How many cpu to use
        verbose             -- Output command information
        noExecute           -- True/False. If True, only print commands and not execute
        use_LSF             -- Use LSF or not
        LSF_parameters      -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmsearch', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    
    Output:
        These files are important:
            * hhblits.psi
            * hhblits_prior.psi
            * hhblits.hmm
            * chkparse.out
            * matlab.out
            * hhblits.a3m
    """
    # Excutable files
    hhblits_exe = General.require_exec("hhblits")
    psiblast_exe = General.require_exec("psiblast")
    matlab_exe = General.require_exec("matlab")
    chkparse_exe = General.require_exec("chkparse")
    
    if len(outdir) > 1:
        outdir = outdir.rstrip('/')
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        print( f"Create {outdir}..." )
        os.mkdir(outdir)
    hhblits_cmd = f"{hhblits_exe} -i {query_fasta} -d {UniClust30_db} -o {outdir}/hhblits.hhr -ohhm {outdir}/hhblits.hmm \
        -opsi {outdir}/hhblits.psi -oa3m {outdir}/hhblits.a3m -n 3 -e 0.001 -cpu {cpu}"
    hhblits_prior_cmd = f"{hhblits_exe} -i {query_fasta} -d {UniClust30_db} -o {outdir}/hhblits_prior.hhr -ohhm {outdir}/hhblits_prior.hmm \
        -opsi {outdir}/hhblits_prior.psi -oa3m {outdir}/hhblits_prior.a3m -n 3 -e 0.001 -nocontxt -cpu {cpu}"
    psi2fas_cmd = f"awk '{{print \">\"$1\"_\"NR\"\\n\"$2}}' {outdir}/hhblits.psi > {outdir}/hhblits.fas"
    # -in_msa {outdir}/hhblits.fas => -query {query_fasta}
    psiblast_cmd = f"{psiblast_exe} -query {query_fasta} -out_ascii_pssm {outdir}/psiblast.mtx -out_pssm {outdir}/psiblast.chk -out {outdir}/psiblast.out \
        -db {Blast_nr_db} -evalue 0.001 -num_iterations 3 -num_threads {cpu}"
    chkparse_cmd = f"{chkparse_exe} {outdir}/psiblast.chk > {outdir}/chkparse.out"
    #ccmpred_cmd = f"{ccmpred_exe} -r {outdir}/ccmpred.raw -n 500 {outdir}/ccmpred.in {outdir}/ccmpred.out -t {cpu}"
    matlab_cmd = f"{matlab_exe} -nodisplay -r \"cd '{plmDCAPATH}'; plmDCA('{outdir}/hhblits.fas', '{outdir}/matlab.mat'); exit;\""
    #awk_cmd = f"awk '{{print $2}}' {outdir}/hhblits.psi > {outdir}/ccmpred.in"
    
    if noExecute:
        verbose = True
        use_LSF = False

    if use_LSF:
        hhblits_job = Cluster.new_job(hhblits_cmd, queue=LSF_parameters.get('queue', 'Z-ZQF'), 
                cpu=LSF_parameters.get('cpu', 20), 
                job_name=LSF_parameters.get('job_name', 'AlphaFold1_hhblits'), 
                logFn=LSF_parameters.get('logFn', f'{outdir}/log'),
                errFn=LSF_parameters.get('errFn', f'{outdir}/error'))
        hhblits_prior_job = Cluster.new_job(hhblits_prior_cmd, queue=LSF_parameters.get('queue', 'Z-ZQF'), 
                cpu=LSF_parameters.get('cpu', 20), 
                job_name=LSF_parameters.get('job_name', 'AlphaFold1_hhblits_prior'), 
                logFn=LSF_parameters.get('logFn', f'{outdir}/log'),
                errFn=LSF_parameters.get('errFn', f'{outdir}/error'))
        psiblast_job = Cluster.new_job(psiblast_cmd, queue=LSF_parameters.get('queue', 'Z-ZQF'), 
                cpu=LSF_parameters.get('cpu', 20), 
                job_name=LSF_parameters.get('job_name', 'AlphaFold1_psiblast'), 
                logFn=LSF_parameters.get('logFn', f'{outdir}/log'),
                errFn=LSF_parameters.get('errFn', f'{outdir}/error'))
        #ccmpred_job = Cluster.new_job(ccmpred_cmd, queue=LSF_parameters.get('queue', 'Z-ZQF'), 
        #        cpu=LSF_parameters.get('cpu', 20), 
        #        job_name=LSF_parameters.get('job_name', 'AlphaFold1_ccmpred'), 
        #        logFn=LSF_parameters.get('logFn', f'{outdir}/log'),
        #        errFn=LSF_parameters.get('errFn', f'{outdir}/error'))
        
        if verbose: print(Colors.f( hhblits_job.get_submit_command(), 'yellow' ))
        hhblits_job.submit()
        if verbose: print(Colors.f( hhblits_prior_job.get_submit_command(), 'yellow' ))
        hhblits_prior_job.submit()
        
        hhblits_job.wait()
        #if verbose: print(Colors.f( awk_cmd, 'yellow' ))
        #os.system(awk_cmd)
        #if verbose: print(Colors.f( ccmpred_job.get_submit_command(), 'yellow' ))
        #ccmpred_job.submit()
        if verbose: print(Colors.f( psi2fas_cmd, 'yellow' ))
        os.system(psi2fas_cmd)
        if verbose: print(Colors.f( psiblast_job.get_submit_command(), 'yellow' ))
        psiblast_job.submit()
        if verbose: print(Colors.f( matlab_cmd, 'yellow' ))
        os.system(matlab_cmd)
        
        psiblast_job.wait()
        if verbose: print(Colors.f( chkparse_cmd, 'yellow' ))
        os.system(chkparse_cmd)
        
        hhblits_prior_job.wait()
        #ccmpred_job.wait()
    else:
        if verbose: print(Colors.f( hhblits_cmd, 'yellow' ))
        if not noExecute: os.system(hhblits_cmd)
        if verbose: print(Colors.f( hhblits_prior_cmd, 'yellow' ))
        if not noExecute: os.system(hhblits_prior_cmd)
        if verbose: print(Colors.f( psi2fas_cmd, 'yellow' ))
        if not noExecute: os.system(psi2fas_cmd)
        if verbose: print(Colors.f( psiblast_cmd, 'yellow' ))
        if not noExecute: os.system(psiblast_cmd)
        #if verbose: print(Colors.f( awk_cmd, 'yellow' ))
        #if not noExecute: os.system(awk_cmd)
        #if verbose: print(Colors.f( ccmpred_cmd, 'yellow' ))
        #os.system(ccmpred_cmd)
        if verbose: print(Colors.f( matlab_cmd, 'yellow' ))
        if not noExecute: os.system(matlab_cmd)
        if verbose: print(Colors.f( chkparse_cmd, 'yellow' ))
        if not noExecute: os.system(chkparse_cmd)

prepare_alphafold1_input = run_alphafold1_input_pipeline

def package_as_DM_AlphaFold_input(instance, chain_name, domain_name, query_seq):
    """
    Package AlphaFold.get_alphafold1_input_array(...) output object to a Deep Mind AlphaFold1 tensorflow input .tfrec file
    
    Parameters:
    instance             -- A object returned by AlphaFold.get_alphafold1_input_array
    chain_name           -- Chain name, 'T0949'
    domain_name          -- Domain name, 'T0949-l128_s55'
    query_seq            -- Query sequence
    
    Return:
        A dict with 33 elements
    
    Example:
        obj = AlphaFold.get_alphafold1_input_array(hhblits_psi_file, hhblits_prior_psi_file, hmm_file, 
            chkparse_file, matlab_file, hhblits_a3m_file, query_seq)
        feature_parse = package_DM_AlphaFold_input(obj)
        example = tf.train.Example(features=tf.train.Features( feature=feature_parse ))
        writer = tf.compat.v1.python_io.TFRecordWriter('output.tfrec')
        writer.write(example.SerializeToString())
        writer.close()
    """
    import tensorflow as tf
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
    
    seq_len = instance.aatype.shape[1]
    feature_parse = {
        'aatype':                   _float_feature( np.ndarray.flatten(instance.aatype) ),
        'gap_matrix':               _float_feature( np.ndarray.flatten(instance.gap_matrix) ), # obj['gap_matrix'],
        'hmm_profile':              _float_feature( np.ndarray.flatten(instance.hmm_profile) ), # obj['hmm_profile'],
        'non_gapped_profile':       _float_feature( np.ndarray.flatten(instance.non_gapped_profile) ), # obj['non_gapped_profile'],
        'pseudo_bias':              _float_feature( np.ndarray.flatten(instance.pseudo_bias) ), # obj['pseudo_bias'],
        'residue_index':            _int64_feature( np.ndarray.flatten(instance.residue_index) ), # obj['residue_index'],
        'profile_with_prior':       _float_feature( np.ndarray.flatten(instance.hhblits_prior_profile) ), # obj['hhblits_prior_profile'],
        'num_alignments':           _int64_feature( np.ndarray.flatten(instance.num_alignments.astype(np.int)) ), # obj['num_alignments'],
        'pseudolikelihood':         _float_feature( np.ndarray.flatten(instance.pseudolikelihood) ), # obj['pseudolikelihood'],
        'reweighted_profile':       _float_feature( np.ndarray.flatten(instance.reweighted_profile) ), # obj['reweighted_profile'],
        'deletion_probability':     _float_feature( np.ndarray.flatten(instance.deletion_probability) ), # obj['deletion_probability'],
        'hhblits_profile':          _float_feature( np.ndarray.flatten(instance.hhblits_profile) ), # obj['hhblits_profile'],
        'profile_with_prior_without_gaps': _float_feature( np.ndarray.flatten(instance.non_gapped_prior_profile) ), # obj['non_gapped_prior_profile'],
        'profile':                  _float_feature( np.ndarray.flatten(instance.profile)/100 ), # obj['profile'],
        'pseudo_frob':              _float_feature( np.ndarray.flatten(instance.pseudolikelihood_frob) ), # obj['pseudolikelihood_frob'],
        ### Added
        'beta_mask':                _int64_feature( np.ndarray.flatten( np.zeros( [seq_len, 1], dtype=np.int ) ) ),
        'beta_positions':           _float_feature( np.ndarray.flatten( np.zeros( [seq_len, 3] ) ) ),
        'between_segment_residues': _int64_feature( np.ndarray.flatten( np.zeros( [seq_len, 1], dtype=np.int ) ) ),
        'chain_name':               _bytes_feature( [str.encode(chain_name)] ),
        'domain_name':              _bytes_feature( [str.encode(domain_name)] ),
        'num_effective_alignments': _float_feature( instance.num_alignments[0][0] ),
        'resolution':               _float_feature( [0.0] ),
        'sec_structure':            _int64_feature( np.ndarray.flatten( np.tile([1,0,0,0,0,0,0,0], (seq_len,1)) ) ),
        'sec_structure_mask':       _int64_feature( np.ndarray.flatten( np.ones( [seq_len, 1], dtype=np.int ) ) ),
        'seq_length':               _int64_feature( [seq_len]*seq_len ),
        'sequence':                 _bytes_feature( [str.encode(query_seq)] ),
        'solv_surf':                _float_feature( np.ndarray.flatten( np.zeros( [seq_len, 1] ) ) ),
        'solv_surf_mask':           _int64_feature( np.ndarray.flatten( np.zeros( [seq_len, 1], dtype=np.int ) ) ),
        'superfamily':              _bytes_feature( [b''] ),
        
        'phi_angles':               _float_feature( np.ndarray.flatten( np.zeros( [seq_len, 1] ) ) ),
        'phi_mask':                 _int64_feature( np.ndarray.flatten( np.ones( [seq_len, 1], dtype=np.int ) ) ),
        'psi_angles':               _float_feature( np.ndarray.flatten( np.zeros( [seq_len, 1] ) ) ),
        'psi_mask':                 _int64_feature( np.ndarray.flatten( np.ones( [seq_len, 1], dtype=np.int ) ) )
    }
    return feature_parse

def get_cath_domain_distogram(cath_id, cath_seq, return_seq_pair=False):
    """
    Given a CATH ID and seq, return distogram correspoding to CATH sequence
    
    Return:
        if return_seq_pair == False:
            return a np.array with shape [ seq_len, seq_len ]. Missing values denoted with -1
        else:
            return np.array, [cath_seq, pdb_seq]
    
    Example:
        get_cath_domain_distogram("2j43A01", "ASHHLRXHFKTLPAGESLGSLGLWVWGDVDQPSKDWPNGAITXTKAKKDDYGYYLDVPLAAKHRQQVSYLINNKAGENLSKDQHISLLTPKXNEVWIDENY")
    """
    import atomium, Structure
    
    ## 分解
    code = cath_id[:4]
    chain = cath_id[4]
    domain = cath_id[5:7]
    
    ## 读取PDB文件
    pdb = atomium.fetch(code)
    chain = pdb.model.chain(chain)
    
    ## 获取PDB文件的序列和CATH数据库序列的关系
    present_sequence = "".join([ threeToOne.get(res.name, "X") for res in chain.residues() ])
    start = present_sequence.find(cath_seq)
    if start != -1:
        end = start + len(cath_seq)
        full_seq = present_sequence
        domain_seq = "-" * start + cath_seq + "-"*(len(present_sequence)-end)
        assert len(domain_seq) == len(full_seq)
    else:
        full_seq, domain_seq = Structure.multi_alignment( [present_sequence,cath_seq]  )
    
    domain_index = [ -1 ] * len(cath_seq)
    fi, di = 0, 0
    for f,d in zip(full_seq, domain_seq):
        if f == '-':
            di += 1
        elif d == '-':
            fi += 1
        else:
            domain_index[ di ] = fi
            di += 1
            fi += 1
    res_list = [ chain.residues()[di] if di!=-1 else None for di in domain_index ]
    
    ## 计算distogram
    distogram = np.zeros([len(res_list), len(res_list)])
    distogram[:] = -1
    for i in range(len(res_list)):
        if res_list[i] is None:
            continue
        if res_list[i].name == 'GLY':
            atom1 = res_list[i].atom(name="CA")
        else:
            atom1 = res_list[i].atom(name="CB")
        if atom1 is None:
            continue
        for j in range(i+1, len(res_list)):
            if res_list[j] is None:
                continue
            if res_list[j].name == 'GLY':
                atom2 = res_list[j].atom(name="CA")
            else:
                atom2 = res_list[j].atom(name="CB")
            if atom2 is None:
                continue
            distogram[i, j] = distogram[j, i] = atom1.distance_to(atom2)
    
    if return_seq_pair:
        pdb_seq = "".join([ threeToOne.get(res.name,'X') if res is not None else '?' for res in res_list ])
        return distogram, ( cath_seq, pdb_seq )
    else:
        return distogram


# input.fasta
# >T0953s1 (72 residues)
# GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS
# run_profile('input.fasta', 'aaaa', cpu=20)

# hhblits_psi_file = 'hhblits.psi'
# hhblits_prior_psi_file = 'hhblits_prior.psi'
# hmm_file = 'hhblits.hmm'
# chkparse_file = 'chkparse.out'
# ccmpred_out_file = 'ccmpred.out'
# ccmpred_raw_file = 'ccmpred.raw'
# hhblits_a3m_file = 'hhblits.a3m'
# query_seq = 'GALGSASIAIGDNDTGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSNGNLNQWGGGAIYCRDLNVS'

# alphafold1_input = get_alphafold1_input_array(hhblits_psi_file, hhblits_prior_psi_file, hmm_file, 
#     chkparse_file, ccmpred_out_file, ccmpred_raw_file, hhblits_a3m_file, query_seq)



