
import sys, os
import numpy as np
import _pickle
from pyrosetta import *
init(silent=True)

#############################
### Read and prepare Distance and Torsions
#############################

scriptdir = "/Users/lee/Project/代码/AlphaFold-post"

bg_distgram_file = scriptdir+'/T0949_bg.pickle'
distgram_file = scriptdir+'/T0949.pickle'
torsion_file = scriptdir+'/T0949.torsions'

distgram_bg             = _pickle.load( open(bg_distgram_file, 'rb') ) # dict
distgram                = _pickle.load( open(distgram_file, 'rb') ) # dict
seq                     = distgram['sequence']
torsion                 = _pickle.load( open(torsion_file, 'rb') ) # dict
torsion_distribution    = torsion = torsion['probs']

nres = len(seq)

dist = distgram['probs']
dist_bg = distgram_bg['probs']

# 每个torsion bin对应的角度
TBIN = 2*np.pi / 36
TBINS = np.array([ TBIN*i+TBIN/2-np.pi for i in range(36) ])

#############################
### Read and preprocess PDB file
#############################

# /Users/lee/Project/代码/AlphaFold-post/models/T0949_mode0_pcut0.05_noorient_Refined.pdb
pose = pose_from_pdb(sys.argv[1])

switch = SwitchResidueTypeSetMover("centroid") # fullatom
switch.apply(pose)

#############################
### 1. Calculate score2_smooth
#############################

score2_smooth_func = pyrosetta.create_score_function("score2_smooth")
V_score2_smooth = score2_smooth_func(pose)
V_score2_smooth /= nres
print( "V_score2_smooth:", V_score2_smooth )

#############################
### 2. Calculate V_Distance
#############################

#V_distance = - np.log((dist+0.00001) / (dist_bg+0.00001))
DSTEP = 0.28125

def get_pose_distogram(pose):
    nres = len(pose.sequence())
    #pose_dist = General.init_list_rect(nres, nres, 0.0)
    pose_dist = np.zeros([nres, nres], dtype=np.float32)
    for a in range(nres):
        atomA = 'CA' if pose.residue(a+1).name().startswith('GLY') else 'CB'
        xyzA = pose.residue(a+1).atom(atomA).xyz()
        for b in range(a+1, nres):
            atomB = 'CA' if pose.residue(b+1).name().startswith('GLY') else 'CB'
            xyzB = pose.residue(b+1).atom(atomB).xyz()
            pose_dist[a, b] = pose_dist[b, a] = xyzA.distance(xyzB)
    return pose_dist

pose_dist = get_pose_distogram(pose)
pose_bins = (pose_dist - 2.0) // DSTEP
pose_bins = pose_bins.astype(np.int32)

count = 0
V_distance = 0
for a in range(nres):
    for b in range(a+1, nres):
        if dist[a, b, :63].sum() < 0.2:
            continue
        count += 1
        bin_ = pose_bins[a, b]
        assert bin_ >= 0
        bin_ = np.clip(bin_, 0, 63)
        V_distance += - np.log((dist[a, b, bin_]+0.00001) / (dist_bg[a, b, bin_]+0.00001))

V_distance /= count
print("V_distance:", V_distance)


#############################
### 3. Calculate V_torsion
#############################

def get_pose_torsion(pose):
    nres = len(pose.sequence())
    pose_torsion = np.zeros([nres, 2], dtype=np.int32)
    for a in range(nres):
        pose_torsion[a, 0] = pose.phi(a+1)
        pose_torsion[a, 1] = pose.psi(a+1)
    return pose_torsion

pose_torsion = get_pose_torsion(pose)
pose_torsion_bins = (pose_torsion + 180) // 10 # phi, psi
torsion_ = torsion.reshape([-1, 36, 36]) # nres, phi, psi

torsion_prob = torsion_[ range(1, nres-1), pose_torsion_bins[1:-1,0], pose_torsion_bins[1:-1,1] ]
V_torsion = -np.sum(np.log( torsion_prob + 0.00001 ))

V_torsion /= nres - 1
print("V_torsion:", V_torsion)

#############################
### 4. Add together
#############################

V_total = V_score2_smooth + V_distance + V_torsion

print("V_total:", V_total)


#############################
### 5. Save
#############################

if os.path.exists(scriptdir + "/score.txt"):
    OUT = open(scriptdir + "/score.txt", 'a')
else:
    OUT = open(scriptdir + "/score.txt", 'w')
    print("Filename\tV_score2_smooth\tV_distance\tV_torsion\tV_total", file=OUT )

print( os.path.basename(sys.argv[1]), V_score2_smooth, V_distance, V_torsion, V_total, sep="\t", end="\n", file=OUT )
OUT.close()


