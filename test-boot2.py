file.close()
import pickle
file=open('rocs-bw-matrix.pickle', 'rb')
data=pickle.load(file)
file.close()
dbase=numpy.loadtxt('recent-all-top-agonist-dbase.list', dtype=str)
improt numpy
import numpy
dbase=numpy.loadtxt('recent-all-top-agonist-dbase.list', dtype=str)
dbase[:5]
states=numpy.loadtxt('agonist-states.txt', dtype=str)
states
'-'.join(x.split('-')[:2]
)
x='car-308-44209280-1'
'-'.join(x.split('-')[:2]
)
newstates=['-'.join(x.split('-')[:2]) for x in states]
newstates
states
dbase
x='apo-1156/stereo-15712586-4_000.mol2'
'-'.join(x.split('/')[0], x.split('stereo-')[0].split('_000')[0])
'-'.join([x.split('/')[0], x.split('stereo-')[0].split('_000')[0]])
'-'.join([x.split('/')[0], x.split('/stereo-')[1].split('_000')[0]])
newdbase=['-'.join([x.split('/')[0], x.split('/stereo-')[1].split('_000')[0]]) for x in dbase]
newdbase
numpy.savetxt('recent-all-top-agonist-reformat-dbase.list', newdbase)
numpy.savetxt('recent-all-top-agonist-reformat-dbase.list', newdbase, fmt='%s')
file=open('agonist-alldata.pickle', 'rb')
alldata=pickle.load(file)
file.close()
alldata.keys9)
alldata.keys()
newstates
newdbase
alldata.keys()
for x in alldata.keys():
    if x in newdbase:
        print x
for x in alldata.keys():
    location=numpy.where(newdbase==x)[0]
    matrix[location]
data.shape
location
print x
numpy.where(newdbase==x)[
)
newdbase=numpy.array(newdbase)
numpy.where(newdbase==x)[
)
numpy.where(newdbase==x)
matrix[197]
data[197]
data[197].shape
alldata
print x
alldata[x]
file=open('new-map-agonist.pickle', 'rb')
mapdata=pickle.load(file)
file.close()
mapdata
alldata[x]
import readline
readline.write_history_file('check-dbase.py')
import scipy
from scipy.sparse import lil_matrix
from scipy.io import *
T=mmread('tProb-19.mtx')
T.shape
T[2,2]
T
T=T.todense()
T
from scipy.io import *
import scipy
from scipy.sparse import lil_matrix
orig=mmread('Data/tProb.mtx')
T=mmread('tProb-18.mtx')
T.shape
T=T.todense()
orig=orig.todense()
T[3,3]
orig[3,3]
orig
T
orig.shape
T.shape
l
ls
from msmbuilder import io
asn=io.loadh('./Data/Assignments.h5)
asn=io.loadh('./Data/Assignments.h5')
asn['arr_0'].max()
import cipy
import scipy
from scipy.io import *
from scipy.sparse import lil_matrix
file='./bi_new_gen/tProb-2.mtx'
T=mmread(file)
from msmbuilder import MSMLib
rev_counts, t_matrix, Populations, Mapping = MSMLib.build_msm(T, symmetrize='MLE', ergodic_trimming=True)
data=dict()
data['bi']=numpy.loadtxt('./bi_new_gen/MSMl15/PopStdev_73iter_55000multi.dat')
import numpy
data['bi']=numpy.loadtxt('./bi_new_gen/MSMl15/PopStdev_73iter_55000multi.dat')
test=numpy.loadtxt('./bi_new_gen/MSMl15/PopStdev_8iter_55000multi.dat')
test=numpy.loadtxt('./bi_new_gen/MSMl15/PopStdev_8check_55000multi.dat')
test
data['bi']
import pylab
pylab.scatter(data['bi'], test)
pylab.show()
systems=['apo', 'bi', 'car']
data['car']=numpy.loadtxt('./car_new_gen/MSMl15/PopStdev_68iter_55000multi.dat')
data['apo']=numpy.loadtxt('./apo_new_gen/MSMl15/PopStdev_8iter_55000multi.dat')
data['apo']
max(data['apo'])
max(data['car'])
max(data['bi'])
import numpy, pylab
data=dict()
data['car']=numpy.loadtxt('./car_new_gen/MSMl15/PopStdev_68iter_55000multi.dat')
data['apo']=numpy.loadtxt('./apo_new_gen/MSMl15/PopStdev_8iter_55000multi.dat')
data['bi']=numpy.loadtxt('./bi_new_gen/MSMl15/PopStdev_73iter_55000multi.dat')
data['bi']
active=dict()
inactive=dict()
active['apo']=numpy.loadtxt('./apo_new_gen/Data_k_3000/TPT-strict/active_states.txt')
inactive['apo']=numpy.loadtxt('./apo_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
active['bi']=numpy.loadtxt('./bi_new_gen/Data_k_3000/TPT-strict/active_states.txt')
active['car']=numpy.loadtxt('./car_new_gen/Data_k_3000/TPT-strict/active_states.txt')
inactive['car']=numpy.loadtxt('./car_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
inactive['bi']=numpy.loadtxt('./bi_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
active
active['car']=numpy.loadtxt('./car_new_gen/Data_k_3000/TPT-strict/active_states.txt')
active
len(active['car'])
inactive['bi']=numpy.loadtxt('./bi_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
inactive['car']=numpy.loadtxt('./car_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
inactive['apo']=numpy.loadtxt('./apo_new_gen/Data_k_3000/TPT-strict/inactive_states.txt')
active['apo']=numpy.loadtxt('./apo_new_gen/Data_k_3000/TPT-strict/active_states.txt')
active['car']=numpy.loadtxt('./car_new_gen/Data_k_3000/TPT-strict/active_states.txt')
active['bi']=numpy.loadtxt('./bi_new_gen/Data_k_3000/TPT-strict/active_states.txt')
active
inactive
for sys in systems:
    active[sys]=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/active_states.txt' % sys, dtype=int)
    inactive[sys]=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/inactive_states.txt' % sys, dtype=int)
systems=['apo', 'bi', 'car']
for sys in systems:
    active[sys]=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/active_states.txt' % sys, dtype=int)
    inactive[sys]=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/inactive_states.txt' % sys, dtype=int)
active
inactive
for sys in systems:
    frames=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/active_states.txt' % sys, dtype=int)
    active[sys]=data[sys][frames]
    frames=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/inactive_states.txt' % sys, dtype=int)
    inactive[sys]=data[sys][frames]
active
len(data['bi'])
for sys in systems:
    frames=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/gen_active_states.txt' % sys, dtype=int)
    active[sys]=data[sys][frames]
    frames=numpy.loadtxt('./%s_new_gen/Data_k_3000/TPT-strict/gen_inactive_states.txt' % sys, dtype=int)
    inactive[sys]=data[sys][frames]
data
active['bi']
active.values()
for sys in systems:
    
    
    bla
populations=dict()
for sys in systems:
    populations[sys]=numpy.loadtxt('%s_new_gen/MSMl15/Populations.dat' % sys)
populations
populations['bi']
populations['bi']/active['bi']
section=dict()
import readline
readline.write_history_file('anapop.py')
from msmbuilder import io, Trajectory
t=Trajectory.load_from_lhdf('./Data_pairs/Gens.lh5')
t.save_to_xtc('./Data_pairs/Gens.xtc')
t=Trajectory.load_from_lhdf('./Data_tica/Gens.lh5')
t.save_to_xtc('./Data_tica/Gens.xtc')
60/60
20/60
20/60.0
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Data_wt_combo2/Gens.lh5')
t.save_to_xtc('./Data_wt_combo2/Gens.xtc')
from msmbuilder import metrics
ligand_inds=numpy.loadtxt('./AtomIndices-ligand.dat')
import numpy
ligand_inds=numpy.loadtxt('./AtomIndices-ligand.dat', dtype=int)
pairmetric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=atom_pairs)
atom_pairs=numpy.loadtxt('./atompairs.dat', dtype=int)
pairmetric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=atom_pairs)
from msmbuilder import Trajectory
pdb=Trajectory.load_from_pdb('fkbpLG2_xtal.pdb')
ppdb=pairmetric.prepare_trajectory(pdb)
from msmbuilder import Project
project = Project.load_from('ProjectInfo.yaml')
traj=project.load_traj(0)
ptraj=pairmetric.prepare_trajectory(traj)
output=metric.one_to_all_aligned(ppdb, ptraj, 0)
output=pairmetric.one_to_all_aligned(ppdb, ptraj, 0)
output=pairmetric.one_to_all(ppdb, ptraj, 0)
output
max(output)
min(output)
traj=project.load_traj(1)
ptraj=pairmetric.prepare_trajectory(traj)
output=pairmetric.one_to_all(ppdb, ptraj, 0)
max(output)
min(output)
rmsd=numpy.loadtxt('./trj-rmsd/run1-aln-dodec-all-prod-run1-LG2_3.rmsd.dat')
min(rmsd)
max(rmsd)
import readline
readline.write_history_file('checkpairs.py')
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
proj.keys()
proj
proj.traj_filename
proj.traj_filename()
proj.traj_filename(0)
proj.traj_filename(1)
proj.traj_filename(2)
proj.traj_filename(9985)
proj.traj_filename(9980)
import pylab. numpy
import pylab, numpy
pylab.figure()
import pylab, numpy
pylab.figure()
for x in range(0, 11467):
    data=numpy.loadtxt('Trajectories-metric/trj%s_pairs.dat' % x)
    data2=numpy.loadtxt('Trajectories-metric/trj%s_gmxhrmsd.dat' % x)
    pylab.scatter(data2, data)
    pylab.hold(True)
print x
data
len(data)
len(data2)
ls Trajectories-metric/trj0_gmxhrmsd.dat
data2=numpy.loadtxt('Trajectories-metric/trj0_gmxhrmsd.dat')
len(data2)
import pylab, numpy
import glob
dist=glob.glob('*.dist.dat')
rmsd=glob.glob('*.gmx*dat')
alldist=[]
allrmsd=[]
len(dist)
dist
dist=glob.glob('*.pairs.dat')
len(dist)
pwd
dist=glob.glob('*pairs.dat')
len(dist)
len(rmsd)
rmsd=glob.glob('*gmx*dat')
len(rmsd)
for x in range(0, 9985):
    data1=numpy.loadtxt('trj%s_pairs.dat' % x)
    data2=numpy.loadtxt('trj%s_gmxhrmsd.dat' % x)
    for i in data1:
        alldist.append(i)
    for j in data2:
        allrmsd.append(j)
import readline
readline.write_history_file('plot_pair_gmxrmsd.py')
import numpy, pylab
import glob
l
ls
ls Trajectories-metric/*pairs.dat | wc -l
rmsd=[]
pairs=[]
for i in range(0, 9985):
    data1=numpy.loadtxt('Trajectories-metric/trj%s_pairs.dat' % i)
    data2=numpy.loadtxt('Trajectories-metric/trj%s_lprmsd.dat' % i)
    for j in data1:
        pairs.append(j)
    for k in data2:
        rmsd.append(k)
len(rmsd)
len(pairs)
pylab.scatter(rmsd, pairs)
pylab.savefig('pair_vs_lprmsd.png', dpi=300)
import readline
readline.write_history_file('plot_pair_rmsd.py')
from msmbuilder import Project
project=Project.load_traj(1091)
project=Project.load_from('ProjectInfo.yaml')
traj=project.load_traj(1091)
traj['XYZList'].shape
frames=[45, 46, 47, 48, 49]
empty=project.empty_traj()
empty['XYZList'].shape
empty['XYZList']
empty['XYZList']=traj['XYZList'][0]
empty.save_to_pdb('test.pdb')
empty.keys()
empty['XYZList']
empty.save_to_xtc('test.xtc')
empty['XYZList'].shape
traj['XYZList'].shape
empty['XYZList']=numpy.zeros((1, 1712, 3))
import numpy
empty['XYZList']=numpy.zeros((1, 1712, 3))
empty['XYZList'][0]=traj['XYZList'][0]
empty.save_to_pdb('test.pdb')
for f in frames:
    empty=project.empty_traj()
    empty['XYZList']=numpy.zeros((1, 1712, 3))
    empty['XYZList'][0]=traj['XYZList'][frame]
    empty.save_to_pdb('weird_t1091_f%s.pdb' % f)
    print "done with %s" % f
for f in frames:
    empty=project.empty_traj()
    empty['XYZList']=numpy.zeros((1, 1712, 3))
    empty['XYZList'][0]=traj['XYZList'][f]
    empty.save_to_pdb('weird_t1091_f%s.pdb' % f)
    print "done with %s" % f
import glob
files=glob.glob('*pair*dat')
all=[]
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        all.append(i)
import numpy
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        all.append(i)
import pylab
pylab.plot(all)
pylab.show()
import glob, numpy, pylab
files=glob.glob('*pair*dat')
all=[]
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        all.append(i)
pylab.plot(all)
pylab.show()
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
proj.traj_lengths
min(proj.traj_lengths)
numpy.mean(proj.traj_lengths)
import numpy
numpy.mean(proj.traj_lengths)
from msmbuilder import io
test=io.loadh('./Data_pairs_d0.5_s100_hybrid/Assignments.h5')
test
from msmbuilder import Project
import numpy
rmsd=numpy.loadtxt('./Data_pairs_d0.5_s100_hybrid/Gens.rmsd.dat')
rmsd
import pylab
pylab.plot(sorted(rmsd))
pylab.show()
pops=dict()
lags=range(100,500,100)
lags
for lag in lags:
    pops[lag]=numpy.loadtxt('./Data_pairs_d0.5_s100_hybrid/msml%s/Populations.dat' % lag)
colors=['red', 'orange', 'green', 'blue']
import readline
readline.write_history_file('color.py')
from msmbuilder import io
test=io.loadh('./Data_lprmsd_d0.3_s100_hybrid/Assignments.h5')
test
from msmbuilder import io
cuts=[0.5 1.0 1.5 2.0]
from msmbuilder import io
cuts=[0.5, 1.0, 1.5, 2.0]
data=dict()
for x in cuts:
    data[x]=io.loadh('Data_pairs_d%s_s100_hybrid/Assignments.h5' % x)
data[0.5]
data[1.0]
data[1.5]
data[2.0]
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
proj.traj_lengths
min(proj.traj_lengths)
numpy.mean(proj.traj_lengths)
import numpy
numpy.mean(proj.traj_lengths)
numpy.std(proj.traj_lengths)
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
import numpy
numpy.where(proj.traj_lengths< 100)
numpy.where(proj.traj_lengths< 100)[0]
len(numpy.where(proj.traj_lengths< 100)[0])
len(numpy.where(proj.traj_lengths< 10)[0])
len(numpy.where(proj.traj_lengths< 20)[0])
len(numpy.where(proj.traj_lengths< 40)[0])
len(numpy.where(proj.traj_lengths< 60)[0])
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
proj.traj_lengths
sum(proj.traj_lengths)
sum(proj.traj_lengths)*100
(sum(proj.traj_lengths)*100)/1000
pwd
cd ../LG2-add/
proj=Project.load_from('ProjectInfo.yaml')
(sum(proj.traj_lengths)*100)/1000
pwd
cd ../AND-add/
proj=Project.load_from('ProjectInfo.yaml')
(sum(proj.traj_lengths)*100)/1000
cd ../LG9-add/
proj=Project.load_from('ProjectInfo.yaml')
(sum(proj.traj_lengths)*100)/1000
from msmbuilder import io
test=io.loadh('./Data_pairs_d0.5_s100_hybrid/Assignments.h5')_
test=io.loadh('./Data_pairs_d0.5_s100_hybrid/Assignments.h5')
test
import pickle
file=open('Data_pairs_d0.5_s100_hybrid/Gens.vmd_ligcoords.pickle', 'rb')
data=pickle.load(file)
file.close()
data
data.keys()
data[0]
data['0']
data['1']
import numpy, pickle
pwd
data=numpy.loadtxt('./Data_pairs_d0.5_s100_hybrid/msml400/prot_lig_distance.dat')
data
data=numpy.loadtxt('./Data_pairs_d0.5_s100_hybrid/msml400/prot_lig_distance.dat')
import numpy, pylab
from msmbuilder import metrics
pairs=numpy.loadtxt)'atompairs.dat')
pairs=numpy.loadtxt('atompairs.dat')
import numpy, pylab
pairs=numpy.loadtxt('atompairs.dat')
pairs
pairs=numpy.loadtxt('atompairs.dat', dtype=int)
metric=metrics.ContinuousContact(contacts=pairs)
from msmbuilder import io
from msmbuilder import Trajectory
gens=Trajectory.load_from_lhdf('./Data_pairs_d0.5_s100_hybrid/Gens.lh5')
distances=metric.prepare_trajectory(gens)
pairs
import readline
readline.write_history_file('check.py')
import check_gen_distance
distances=check_gen_distance.main()
distances
distances.shape
distances[,:]
distances[:,0]
more atompairs.dat
distances[0]
distances[1]
distances[2]
distances[3]
distances[4]
distances[5]
distances[6]
import check_gen_distance
distances=check_gen_distance.main()
distances[1]
distances[0]
min(distances[0])
min(distances[2])
min(distances[3])
min(distances[4])
min(distances[5])
distances.shape
for i in range(0, distances.shape[0]):
    
    bla
all=[]
for i in range(0, distances.shape[0]):
    all.append(numpy.sqrt(distances[0][i]**2+distances[1][i]**2+distances[2][i]**2+distances[3][i]**2+distances[4][i]**2+distances[5][i]**2))
import numpy
for i in range(0, distances.shape[0]):
    all.append(numpy.sqrt(distances[0][i]**2+distances[1][i]**2+distances[2][i]**2+distances[3][i]**2+distances[4][i]**2+distances[5][i]**2))
distances[0]
distances[5]
distances[6]
distances[7]
distances.shape
for i in range(0, distances.shape[0]):
    all.append(numpy.sqrt(distances[i][0]**2+distances[i][1]**2+distances[i][2]**2+distances[i][3]**2+distances[i][4]**2+distances[i][5]**2))
all
max(all)
min(all)
all=numpy.array(all)
numpy.where(all==0)
all[7]
distances[7,:]
distances.shape
distances[7,:]
sum(distances[7,:])
sum([i**2 for i in distances[7,:]])
numpy.sqrt(sum([i**2 for i in distances[7,:]]))
all=[]
for i in range(0, distances.shape[0]):
    all.append(numpy.sqrt(sum([j**2 for j in distances[i,:]]))
    
    )
all
min(all)
all=numpy.array(all)
numpy.where(all==0)
all[1]
distances[1,:]
distances
distances.shape
distances[7]
distances[20]
distances[20,:]
distances[1,:]
distances[0]
distances[1]
distances[2]
import readline
readline.write_history_file('check2')
import glob, numpy
files=glob.glob('./Trajectories-metric/*pair*dat')
problem=[]
for file in files:
    data=numpy.loadtxt(file)
    location=numpy.where(data==0)[0]
    if location.size:
        problem.append(file)
problem
from msmbuilder.metrics import core
from msmbuilder import Project
p=Project.load_from('ProjectInfo.yaml')
p.traj_lengths
min(p.traj_lengths)
sorted(p.traj_lengths)
sorted(p.traj_lengths)[:10]
sorted(p.traj_lengths)[:20]
sorted(p.traj_lengths)[:40]
sorted(p.traj_lengths)[:100]
numpy.where(p.traj_lengths < 50)
import numpy
numpy.where(p.traj_lengths < 50)
numpy.where(p.traj_lengths < 50)[0]
len(numpy.where(p.traj_lengths < 50)[0])
frames=numpy.where(p.traj_lengths < 50)[0]
p.traj_lengths[frames]
sum(p.traj_lengths[frames])*100
import pylab, numpy
pylab.plot(exp(-1))
pylab.plot(numpy.exp(-1))
pylab.show()
-0.6*numpy.log(0.90/0.01)
-0.6*numpy.log(0.99/0.001)
-0.6*numpy.log(0.999/0.0001)
-0.6*numpy.log(0.9999/0.00001)
-0.6*numpy.log(0.99999/0.000001)
-0.6*numpy.log(0.999999/0.0000001)
def testprior(x,a, max):
    return 1.0/(x+a)*numpy.log(a+max)/a
import random
sample=random_integers(2,40,1000)
sample=random.random_integers(2,40,1000)
sample=numpy.random.random_integers(2,40,1000)
import numpy
sample=numpy.random.random_integers(2,40,1000)
max(sample)
output=[testprior(x, 10, 40) for x in sample]
pylab.plot(output)
import pylab
pylab.plot(output)
pylab.show()
order=numpy.argsort(sample)
sorted_sample=sample[order]
sorted_prior=output[order]
output=numpt.array(output)
output=numpy.array(output)
sorted_prior=output[order]
pylab.plot(sorted_sample, sorted_prior)
pylab.show()
import readline
readline.write_history_file('jeff_prior.py')
systems='LG2 LG9 AND'
data=dict()
import glob, numpy, pylab
for sys in systems:
    data[sys]=[]
    files=glob.glob('%s-step/Trajectories-metric/*pair*dat')
    for file in files:
        tmp=numpy.loadtxt(file)
        for i in file:
            data[sys].append(i)
data.keys()
print file
files
pwd
cd ..
data=dict()
for sys in systems:
    data[sys]=[]
    files=glob.glob('%s-step/Trajectories-metric/*pair*dat')
    for file in files:
        tmp=numpy.loadtxt(file)
        for i in file:
            data[sys].append(i)
data
pwd
data=dict()
for sys in systems:
    data[sys]=[]
    files=glob.glob('%s-step/Trajectories-metric/*pair*dat' % sys)
    for file in files:
        tmp=numpy.loadtxt(file)
        for i in file:
            data[sys].append(i)
data
import readline
readline.write_history_file('check-metrics.py')
from msmbuilde import io
from msmbuilder import io
test=io.loadh('Data_lprmsd_d0.3_s100_hybrid/Assignments.h5')
test
import numpy
com=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/Gens.vmd_com.dat', usecols=(1,))
map=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/msml150/Mapping.dat')
map_com=[]
for i in range(0, len(com)):
    if map[i]!=-1:
        map_com.append(com[i])
com
len(com)
len(map)
com=com[1:]
map_com=[]
for i in range(0, len(com)):
    if map[i]!=-1:
        map_com.append(com[i])
map_com
binary=[]
for x in map_com:
    if x > 20:
        binary.append(0)
    else:
        binary.append(1)
binary
numpy.savetxt('Data_lprmsd_d0.3_s100_hybrid/msml150/mapped_binary.dat', binary)
numpy.where(binary==1)
binary=numpy.array(binary)
numpy.where(binary==1)
numpy.where(binary==1)[0]
len(numpy.where(binary==1)[0])
len(numpy.where(binary==0)[0])
sum(binary)
import readline
readline.write_history_file('test.py')
from msmbuilder import io
test=io.loadh('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_d1.2_Assignments_test.h5')
test
max(test)
max(test['arr_0'])
max(test['arr_0'].flatten())
names=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_gen/Coarsed1.2_r10_Gens_test.rmsd.dat', usecols=(0,))
import numpy
names=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_gen/Coarsed1.2_r10_Gens_test.rmsd.dat', usecols=(0,))
names=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_gen/Coarsed1.2_r10_Gens_test.rmsd.dat', usecols=(0,), dtype=str)
rmsd=numpy.loadtxt('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_gen/Coarsed1.2_r10_Gens_test.rmsd.dat', usecols=(2,))
rmsd
names
frames=numpy.where(names=='unbound')[0]
unbound=rmsd[frames]
frames=numpy.where(names=='bound')[0]
bound=rmsd[frames]
max(rmsd)
bins=range(0,55,5)
pylab.hist(bound, bins, label='bound')
import pylab
pylab.hist(bound, bins, label='bound')
pylab.hold(True)
pylab.hist(unbound, bins, label='unbound')
pylab.legend()
pylab.show()
ass
from msmbuilder import io
ass=io.loadh('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r10_d1.2_Assignments_test.h5')
ass
import numpy
data=numpy.loadtxt('Data_lprmsd_d1.2_s100_hybrid/Gens.rmsd.dat')
data
min(data)
pylab.plot(data)
import pylab
pylab.plot(data)
pylab.show()
pops=numpy.loadtxt('./Data_lprmsd_d1.2_s100_hybrid/msml100/Populations.dat')
map=numpy.loadtxt('./Data_lprmsd_d1.2_s100_hybrid/msml100/Mapping.dat')
new=[]
for i in range(0, len(data)):
    if map[i]!=-1:
        new.append(data[i])
new
new[0]
pylab.scatter(new, pops)
pylab.show()
new=numpy.array(new)
numpy.where(new < 10)
frames=numpy.where(new < 10)[0]
new[frames]
sum(pops[frames])
frames=numpy.where(new < 15)[0]
sum(pops[frames])
frames=numpy.where(new < 20)[0]
frames=numpy.where(new < 15)[0]
frames=numpy.where(new < 20)[0]
sum(pops[frames])
new
binary=[]
for i in new:
    if i > 15.0:
        binary.append(0)
    else:
        binary.append(1)
numpy.savetxt('./Data_lprmsd_d1.2_s100_hybrid/msml100/mapped_binary.dat', binary)
from msmbuilder import io
ass=io.loadh('Data_lprmsd_d0.3_s100_hybrid/Coarsed_r15_d1.2_Assignments.h5')
numpy.where(ass['arr_0'])==1389)
numpy.where(ass['arr_0']==1389)
import numpy
numpy.where(ass['arr_0']==1389)
numpy.where(ass['arr_0']==1458)
numpy.where(ass['arr_0']==1458)[0]
set(numpy.where(ass['arr_0']==1458)[0])
numpy.where(ass['arr_0']==1458)[0]
numpy.where(ass['arr_0']==1389)
import pylab, numpy
axis=numpy.loadtxt('Coarsed1.2_r10_Gens_k448_l100_intdistance_axis.dat')
volume=numpy.loadtxt('Coarsed1.2_r10_Gens_k448_l100_intdistance_volume.dat')
free=numpy.loadtxt('Coarsed1.2_r10_Gens_k448_l100_intdistance_free.dat')
pylab.plot(axis, free)
pylab.show()
import readline
readline.write_history_file('plot.py')
from msmbuilder import Projectt
from msmbuilder import Project
p=Project.load_from('ProjectInfo.yaml')
p.traj_lengths
sum(p.traj_lengths)
pwd
cd ../LG9-step/
p=Project.load_from('ProjectInfo.yaml')
sum(p.traj_lengths)
(2003592*100)/1000
(2003592*100)/(1000*1000)
(2003592*100)/(1000*1000.0)
(2054250*100)/(1000*1000.0)
pwd
cd ../AND-step/
p=Project.load_from('ProjectInfo.yaml')
(sum(p.traj_lengths)*100)/(1000*1000.0)
cd ../LG2-step
p=Project.load_from('ProjectInfo.yaml')
cd ../LG2-step-new/
p=Project.load_from('ProjectInfo.yaml')
(sum(p.traj_lengths)*100)/(1000*1000.0)
cd ../LG6-step/
p=Project.load_from('ProjectInfo.yaml')
(sum(p.traj_lengths)*100)/(1000*1000.0)
import readline
readline.write_history_file('check_projects.py')
import pylab
omega=[]
phi=[]
import numpy, lobg
import numpy, glob
pwd
cd ../LG3-step/
files=glob.glob('Trajectories-metric/*omega*dat')
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        omega.append(i)
files=glob.glob('Trajectories-metric/*phi*dat')
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        phi.append(i)
max(phi)
min(phi)
H,xedges,yedges=numpy.histogram2d(x_val,y_val,bins=360, range=((-180, 180), (-180, 180)), normed=True)
H,xedges,yedges=numpy.histogram2d(omega,phi,bins=360, range=((-180, 180), (-180, 180)), normed=True)
pylab.pcolor(H)
pylab.show()
pylab.pcolor(H.transpose())
pylab.show()
H2,xedges2,yedges2=numpy.histogram2d(omega,phi,bins=180, range=((0, 180), (0, 180)), normed=True)
pylab.pcolor(H.transpose())
pylab.show()
max(phi)
min(phi)
pylab.plot(hi)
pylab.plot(phi)
pylab.show()
len(phi)
pylab.hist(phi)
H,x=pylab.hist(phi)
(H,x)=pylab.hist(phi)
H=pylab.hist(phi)
H.shape
H[0]
H[1]
pylab.hist(phi)
pylab.xlabel(H[1])
pylab.show()
pylab.hist(phi)
pylab.show()
pylab.plot([i for i in phi if i % 2==0])
pylab.show()
len(phi)
pylab.plot([i for i in phi if i % 10==0])
pylab.show()
pylab.plot([i for i in phi if i % 20==0])
pylab.show()
import pylab
import readline
readline.write_history_file('check_phi_omega.py')
from msmbuilder import io
ass=io.loadh('./Data_lprmsd_d0.3_s100_hybrid/Coarsed_r15_d1.2_Assignments.h5')
ass['arr_0'].shape
ls Data_lprmsd_d0.3_s100_hybrid/msml200
ass['arr_0'].shape
from msmbuilder import Project
p=Project.load_from('ProjectInfo.yaml')
p.traj_lengths
min(p.traj_lengths)
numpy.where(p.traj_lengths< 20)
import numpy
numpy.where(p.traj_lengths< 20)
frames=numpy.where(p.traj_lengths< 20)[0]
len(frames)
p.traj_lengths
len(p.traj_lengths)
max(p.traj_lengths)
data=dict()
new=max(p.traj_lengths)-100
new
for i in ass['arr_0'].shape[0]:
    bla
data[new]=-numpy.ones(ass['arr_0'].shape[0], new)
for i in range(0, ass['arr_0'].shape[0]):
    data[n][i]=ass['arr_0'][i][:new]
for i in ass['arr_0'].shape[0]:
    data[i]=ass['arr_0'][i][:new]
for i in range(0, ass['arr_0'].shape[0]):
    data[i]=ass['arr_0'][i][:new]
data
data.shape
import readline
readline.write_history_file('subsample.py')
from msmbuider import io
from msmbuilder import io
test=io.loadh('./Coarsed_r15_d1.2_Assignments_sub-989.h5')
test.shape
test['arr_0'].shape
test
ls
rm *
from msmbuilder import io
test=io.loadh('./Data_lprmsd_d0.3_s100_hybrid/Coarsed_r15_d1.2_Assignments_sub771.h5')
test['arr_0'].shape
from msmbuilder import Project
proj=Project.load_from('../LG3-step/ProjectInfo.yaml')
proj.num_trajs
proj.num_trajs()
proj=Project.load_from('../LG3-step/ProjectInfo.yaml')
proj.n_trajs
proj.traj_lengths
sum(proj.traj_lengths)
import glob, numpy
files=glob.glob('*_omega*')
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        if i < 50 and i > 0: 
            print file
data=numpy.loadtxt('trj3305_omega.dat'
)
pylab.plot(data)
import pylab
pylab.plot(data)
pylab.show()
files=glob.glob('*_omega*')
for file in files:
    data=numpy.loadtxt(file)
    for i in data:
        if i < 50 and i > 0: 
            print file
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Trajectories/trj5366.lh5')
t.save_to_xtc('./trj5366.xtc')
t=Trajectory.load_from_lhdf('./Trajectories/trj3305.lh5')
t.save_to_xtc('./trj3305.xtc')
from msmbuilder import glob, numpy
import glob, numpy, pylab
pairs=numpy.loadtxt('Trajectories-metric/all-pairs-s1.dat')
lprmsd=numpy.loadtxt('Trajectories-metric/all-lprmsd-s1.dat')
pylab.scatter(lprmsd, pairs)
pylab.show()
files=glob.glob('./Trajectories-metric/trj*rmsd*dat')
lprmsd=[]
pairs=[]
for file in files:
    data=numpy.loadtxt(file)
    pairs=numpy.loadtxt('%s_pairs.dat' % file.split('_lp')[0])
    for (i,j) im zip(data, pairs):
        if i < 20 and j < 5.0
        :
            print file
import readline
readline.write_history_file('checkprob.py')
import numpy
x=numpy.loadtxt('../AND-step/Data_pairs_d0.5_s100_hybrid/Gens.vmd_ligcoords.dat', usecols=(1,))
y=numpy.loadtxt('../AND-step/Data_pairs_d0.5_s100_hybrid/Gens.vmd_ligcoords.dat', usecols=(2,))
z=numpy.loadtxt('../AND-step/Data_pairs_d0.5_s100_hybrid/Gens.vmd_ligcoords.dat', usecols=(3,))
min(x)
max(x)
from msmbuilder import Project
p=Project.load_from('./LG2-step-new/ProjectInfo.yaml')
p.traj_lengths
max(p.traj_lengths)
sorted(p.traj_lengths)
sorted(p.traj_lengths)[::-1]
sorted(p.traj_lengths)[::-1][:10]
list=['x1', 'y2', 'z3', 'y1', 'z1', 'x2', 'z2', 'x3', 'y3']
vector=[ 6.98745,   6.98745 ,  4.94088 ,  0.00000 ,  0.00000 ,  0.00000  , 0.00000 ,  3.49373  , 3.49373]
len(vector)
len(list)
for (a,b) in zip(list, vector):
    print a,b
import numpy
boxvolume=244.80*(10**3)
v0=1600
corr=-0.6*numpy.log(boxvolume/v0)
corr
from msmbuilder import io
test=io.loadh('../LG2-step/Data_pairs_d0.5_s100_hybrid/subsample/Assignments_sub611.h5')
test
from msmbuilder import Project
p=Project.load_from('../LG2-step/ProjectInfo.yaml')
p.traj_lengths
min(p.traj_lengths)
numpy.where(p.traj_lengths< 15)
import numpy
numpy.where(p.traj_lengths< 15)
frames=numpy.where(p.traj_lengths< 15)[0]
len(frames)
frames=numpy.where(p.traj_lengths > 100)[0]
len(frames)
p.traj_lengths
frames=numpy.where(p.traj_lengths > 100)[0]
len(frames)
frames=numpy.where(p.traj_lengths > 200)[0]
len(frames)
p.traj_lengths
numpy.mean(p.traj_lengths)
from msmbuilder import io
from scipy.io import *
from scipy.sparse import lil_matrix
modeldir='../LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/')
modeldir='../LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/'
T=mmread('%s/tProb.mtx' % modeldir)
from scipy.io import *
import numpy
modeldir='LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/'
T=mmread('%s/tProb.mtx' % modeldir)
from msmbuilder import io
ass=io.loadh('LG2-step/Data_pairs_d0.5_s100_hybrid/Assignments.h5')
T[3442,2683]
T=T.todense()
T[3442,2683]
T[3442,:]
max(T[3442,:])
T[3442,:]
row=T[3442,:]
row.max()
row=numpy.array(row)
numpy.where(row!=0)
frames=numpy.where(row!=0)[0]
row[frames]
frames=numpy.where(row!=0)
frames
frames=numpy.where(row!=0)[0]
row[frames]
state1=3442
state2=2683
state3=4619
T[state1,:]
T[state1,:].max()
numpy.where(ass['arr_0']==state1)
numpy.where(ass['arr_0']==state1).shape
numpy.where(ass['arr_0']==state1)
numpy.where(ass['arr_0']==state1)[0]
trajs=numpy.where(ass['arr_0']==state1)[0]
frames=numpy.where(ass['arr_0']==state1)[1]
for i in trajs:
    for f in frames:
        print ass['arr_0'][i][f-5:f+5]
state2
for i in trajs:
    for f in frames:
        list=ass['arr_0'][i][f-5:f+5]
        if state2 in list:
            print list
state1
state2
trajs
for i in trajs:
    frames=numpy.where(ass['arr_0'][i]==i)
    for f in frames:
        if state2 in ass['arr_0'][i][f-2:f+2]:
            print i, ass['arr_0'][i][f-2:f+2]
f-2
print f
for i in range(0, arr['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)
    for f in frames:
        if state2 in ass['arr_0'][i][f-2:f+2]:
            print i, ass['arr_0'][i][f-2:f+2]
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)
    for f in frames:
        if state2 in ass['arr_0'][i][f-2:f+2]:
            print i, ass['arr_0'][i][f-2:f+2]
print f
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)
    if frames.size:
        for f in frames:
            if state2 in ass['arr_0'][i][f-2:f+2]:
                print i, ass['arr_0'][i][f-2:f+2]
frames
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)[0]
    if frames.size:
        for f in frames:
            if state2 in ass['arr_0'][i][f-2:f+2]:
                print i, ass['arr_0'][i][f-2:f+2]
state2
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)[0]
    if frames.size:
         for f in frames:
            print i, ass['arr_0'][i][f-2:f+2]
state2
ass=io.loadh('LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/Assignments.Fixed.h5')
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)[0]
    if frames.size:
        for f in frames:
            if state2 in ass['arr_0'][i][f-2:f+2]:
                print i, ass['arr_0'][i][f-2:f+2]
state1
state2
adjacent=[]
for i in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][i]==state1)[0]
    if frames.size:
        for f in frames:
            select=ass['arr_0'][i][f-2:f+2]
            for j in select:
                adjacent.append(j)
adjacent
state2
numpy.where(adjacent==state2)
adjacent=numpy.array(adjacent)
numpy.where(adjacent==state2)
from msmbuilder import io
from scipy.io import io
from scipy.io import *
state1=3442
state2=2683
Tp=mmread('%s/tProb.mtx' % modeldir)
modeldir='LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/'
Tp=mmread('%s/tProb.mtx' % modeldir)
Tc=mmread('LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/tCounts.mtx')
Tc=Tc.todense()
Tp=Tp.todense()
Tp[state1, state2]
Tc[state1, state2]
Tc[state1]
Tc[state1].shape
Tc[state1].max()
ass=io.loadh('LG2-step/Data_pairs_d0.5_s100_hybrid/msml100/Assignments.Fixed.h5')
numpy.where(ass['arr_0']==state2)
import numpy
numpy.where(ass['arr_0']==state2)
numpy.where(ass['arr_0']==state2)[0]
frames1=numpy.where(ass['arr_0']==state2)[0]
frames2=numpy.where(ass['arr_0']==state2)[1]
for (i,j) in zip(frames1, frames2):
    print ass['arr_0'][i][j-2:j+2]
state2
state1
adjacent=[]
for (i,j) in zip(frames1, frames2):ass['arr_0'][i][j-2:j+2]
for (i,j) in zip(frames1, frames2):
    states=ass['arr_0'][i][j-2:j+2]
    for k in states:
        adjacent.append(k)
adjacent=numpy.array(adjacent)
numpy.where(adjacent==state1)
import numpy, glob
files=glob.glob('*oxo*dat')
for file in files:
    data=numpy.loadtxt(file)
    frames=numpy.where(data > 10.0)[0]
    if frames.size:
        print file
data=numpy.loadtxt('trj9647_oxos_dist.dat')
import pylab
pylab.plot(data)
pylab.show()
badfiles=[]
for file in files:
    data=numpy.loadtxt(file)
    frames=numpy.where(data > 10.0)[0]
    if frames.size:
        badfiles.append(file)
badfiles
pylab.figure()
for file in badfiles:
    data=numpy.loadtxt(file)
    pylab.plot(file)
    print file
    pylab.show()
ls
head trj2133_hydrophob_dist.dat
!head trj2133_hydrophob_dist.dat
badfiles
for file in files:
    data=numpy.loadtxt(file)
    pylab.figure()
    pylab.plot(data)
    print file
    pylab.show()
import readline
readline.write_history_file('checkoxo.py')
test=numpy.array([1,3,4,4,5,6,7])
import numpy
test=numpy.array([1,3,4,4,5,6,7])
len(test)
test=numpy.array([1,3,4,4,5,6,7,2])
[i for (n,i) in enumerate(test) if n % 2==0]
dir='/nobackup/mlawrenz/FKBP-FAH-results2/%s-step/Data_pairs_d0.5_s100_hybrid/subsample/' % sys
sys='LG3
sys='LG3;
sys='LG3'
dir='/nobackup/mlawrenz/FKBP-FAH-results2/%s-step/Data_pairs_d0.5_s100_hybrid/subsample/' % sys
import glob
glob.glob('%s/sub*/' % dir)
import numpy
-0.6*numpy.log(883/222.0)
-0.6*numpy.log(891/273.0)
-0.6*numpy.log(560/185)
-0.6*numpy.log(494/141)
-0.6*numpy.log(411/120)
import pylab,numpy
import random
random.random(0)
random.random()
A=[random.random() for i in range(0,1000)]
B=[random.random() for i in range(0,1000)]
pylab.plot(A,B)
pylab.show()
import numpy, scipy
n_states=10
C = scipy.sparse.lil_matrix((int(n_states), int(n_states)))
from scipy import *
C = scipy.sparse.lil_matrix((int(n_states), int(n_states)))
!more parallell-bootstrap_Tonly_2.6_slide.py
!more parallel-bootstrap_Tonly_2.6_slide.py
from scipy.sparse import lil_matrix
C = scipy.sparse.lil_matrix((int(n_states), int(n_states)))
C
C[0]
C=C.todense()
C
import random
random.random(10)
random.random(10,10)
random.random()
for i in range(0, n_states):
    for j in range(0, n_states):
        C[i,j]=random.random()
C
for i in range(0, n_states):
 bla
import readline
readline.write_history_file('test-boot2.py')
