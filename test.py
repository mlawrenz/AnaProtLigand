import numpy
dist=numpy.loadtxt('all_dist')
rmsd=numpy.loadtxt('all_rmsd')
import PRoject
from msmbuilder import Project
cd LG3-new/
proj=Project.load_from('ProjectInfo.yaml')
sum(proj.traj_lengths)
pwd
cd ../LG2-new/
dist=numpy.loadtxt('all_dist')
import numpy
dist=numpy.loadtxt('all_dist')
300/2000.0
-0.6*numpy.log(2000*10**(-9))
import numpy
-0.6*numpy.log(2000*10**(-9))
-0.6*numpy.log(660*10**(-9))
60/660.0
-0.6*numpy.log(12*10**(-9))
5/12.0
-0.6*numpy.log(7*10**(-9))
2/7.0
import numpy
dist=numpy.loadtxt('all_rmsd')
rmsd=dist
dist=numpy.loadtxt('all_dist')
import pylab
pylab.scatter(dist, rmsd)
exiy
import numpy, pylab
rmsd=numpy.loadtxt('all_rmsd')
dist=numpy.loadtxt('all_dist')
pylab.scatter(rmsd, dist)
pylab.show()
pylab.scatter(range(0, len(rmsd)), rmsd)
pylab.show()
import numpy, pylab
dist=numpy.loadtxt('all_dist')
rmsd=numpy.loadtxt('all_rmsd')
pylab.scatter(range(0, len(rmsd)), rmsd)
pylab.show()
rmsd[:10]
max(rmsd)
min(rmsd)
import numpy, pylab
rmsd=numpy.loadtxt('all_rmsd')
pylab.scatter(range(0, len(rmsd)), rmsd)
pylab.show()
pylab.scatter(range(0, len(dist)), dist)
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
pylab.show()
cd ../LG2-new/
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
dist=numpy.loadtxt('all_dist')
pylab.show()
cd ../LG6-new/
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
pylab.show()
cd ../LG9-new/
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
pylab.show9)
pylab.show()
import numpy, pylab
cd ../LG9-new/
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
import numpy, pylab
 dist=numpy.loadtxt('all_dist')
dist=numpy.loadtxt('all_dist')
pylab.scatter(range(0, len(dist)), dist)
pylab.show()
from msmbuilder import Project
cd LG3-new/
proj=Project.load_from('ProjectInfo.yaml')
len(proj.traj_lengths)
sum(proj.traj_lengths)
import numpy, pylab
import numpy, pylab
xmin=0
ymin=0
zmin=0
from test-orig import get_space
from test_orig import get_space
import numpy
GD=get_space()
from test_orig import get_space
GD=get_space()
from msmbuilder import Project
p=Project.load_from('ProjectInfo.yaml')
traj=p.load_traj(1)
from msmbuiler import mstrics
from msmbuilder import metrics
mymetric=metrics.Dihedral(metric='euclidean', p=2, angles=['chi'])
ptraj=mymetric.prepare_trajectory(traj)
ptraj=mymetric.prepare_trajectory(traj, degrees=True)
ptraj=mymetric.prepare_trajectory(traj)
traj['XYZList']
traj['XYZList'].shape
traj['ResidueName'].shape
import numpy
xyz=numpy.random.random((100,3))
xuz[:,:2]
xyz[:,:2]
xyz[:,:2].shape
xyz.shape
import numpy
test=numpy((2,3,4))
test=numpy.zeros((2,3,4))
test[0,2]=4.0
test
test[0,2,0]=4.0
test[0,2,2]=3.0
test[0,2,1]=3.0
test[2,2,1]=3.0
test[1,2,1]=3.0
test.shape
test=numpy.zeros((3,3,3))
test[1,1,1]=4.0
test
test[2,2,2]=1.0
test[1,1,1]=1.0
test
test.shape
test[1,0,0]=1.0
new=test.transpose()
test[0,0,1]
test[1,0,0]
test
new[0,0,1]
from msmbuilder import Project
48 % 24
48 % 23
49 % 24
import numpy
x_coor=numpy.loadtxt('./Data_plaink1000/Gens.vmdallcoords.dat', usecols=(1,))
y_coor=numpy.loadtxt('./Data_plaink1000/Gens.vmdallcoords.dat', usecols=(2,))
z_coor=numpy.loadtxt('./Data_plaink1000/Gens.vmdallcoords.dat', usecols=(3,))
min(x_coor)
max(x_coor)
min(y_coor)
max(y_coor)
min(z_coor)
max(z_coor)
import pylab, numpy
pylab.figure.func_defaults
fig=pylab.figure()
from msmbuilder import Project
cd LG2-new/
proj=Project.load_from('ProjectInfo.yaml')
sum(proj.traj_lengths)
cd ../AND-new/
proj=Project.load_from('ProjectInfo.yaml')
proj.traj_lengths
sum(proj.traj_lengths)
import numpy
test=numpy.zeros((5))
test
test[0]='test'
test=numpy.zeros((5), dtype='str')
test[0]='test'
test
test=numpy.zeros((5), dtype='a5')
test[0]='test'
test
import numpy
ass=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-assignments.dat')
ass.shape
ass[:20]
dist=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-distances.dat')
dist[:20]
import numpy
dist=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-distances.dat')
dist[:20]
dist=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-distances.dat')
ass=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-assignments.dat')
ass[:20]
import numpy
dist=numpy.loadtxt('../apo-2665/zinc/reports/dbase_1/dbase_1-distances.dat')
dist[:20]
import numpy
ls
ls re
ls reports/
ls reports/dbase_1/*ass*
numbers=range(1,99)
for num in numbers:
    ba
    
for num in numbers:
    data=numpy.loadtxt('./reports/dbase_%s/dbase_%s-assignments.dat' % (num, num), dtype=int)
    if num==1:
        alldata=data
    else:
        alldata=numpy.stack((alldata, data))
        
for num in numbers:
    data=numpy.loadtxt('./reports/dbase_%s/dbase_%s-assignments.dat' % (num, num))
    if num==1:
        alldata=data
    else:
        alldata=numpy.stack((alldata, data))
        
for num in numbers:
    data=numpy.loadtxt('./reports/dbase_%s/dbase_%s-assignments.dat' % (num, num))
    if num==1:
        alldata=data
    else:
        alldata=numpy.hstack((alldata, data))
        
import readline
readline.write_history_file('combine_dbases.py')
import glob
files=glob.glob('./reports/dbase_*/*ass*')
files=glob.glob('./reports/dbase_*/*ass*')
for file in files:
    data=numpy.loadtxt(file)
    data.shape
    
import numpy
for file in files:
    data=numpy.loadtxt(file)
    data.shape
    
for file in files:
    data=numpy.loadtxt(file)
    print data.shape
    
total=[]
for file in files:
    data=numpy.loadtxt(file)
    total.append(data.shape[0])
    
total
sum(total)
import numpy
file='reports/dbase_99/dbase_99-assignments.dat'
data=numpy.loadtxt(file)
data.shape
import numpy
gen_inds=[249, 95, 231, 161, 195]
gens=numpy.loadtxt('./reports/gens.list', dtype=str)
gens[gen_inds]
ls ../*gen*
bindb=numpy.loadtxt('../bindb_combo_1.5_gen_names.dat')
bindb=numpy.loadtxt('../bindb_combo_1.5_gen_names.dat', dtype=str)
gdd=numpy.loadtxt('../gdd_combo_1.5_gen_names.dat', dtype=str)
gdd_ind=numpy.loadtxt('../gdd_combo_1.5_gen_indices.dat', dtype=int)
bindb_ind=numpy.loadtxt('../bindb_combo_1.5_gen_indices.dat', dtype=int)
gens[gen_inds]
names=gens[gen_inds]
for x in names:
    location=numpy.where(bind==x)[0]
    if not location.size:
        location=numpy.where(gdd==x)[0]
        if not location.size:
            print "problem for ", x
        else:
            "in gdd ", x
    else:
        "in bindb ", x
        
import readline
readline.write_history_file('get_top_gens.py')
from msmbuilder import io
from msmbuilder import Project
p=Project.load_from('./LG2-new/ProjectInfo.yaml')
from msmbuilder import Trajectory
from msmbuilder import Trajectory
gens=Trajectory.load_from_xtc('./Data_plaink1000/Gens.xtc')
gens=Trajectory.load_from_xtc('./Data_plaink1000/Gens.xtc', conf='LG2_noh.pdb')
gens=Trajectory.load_from_xtc('./Data_plaink1000/Gens.xtc', Conf='LG2_noh.pdb')
gens=Trajectory.load_from_lhdf('./Data_plaink1000/Gens.lh5')
gens['XYZList'].shape
new=gens['XYZList'].copy()
new=gens.copy()
new['XYZList']=gens['XYZList'][1]
new['XYZList'].shape
from msmbuilder import Project
proj=Project.get_random_confs_from_states()
proj=Project.get_random_confs_from_states('./Data_plaink1000/Assignments.h5', 1, 1)
proj=Project.get_random_confs_from_states('./Data_plaink1000/Assignments.h5', 1, 1)
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
test=proj.get_random_confs_from_states('./Data_plaink1000/Assignments.h5', 1, 100)
test=proj.get_random_confs_from_states(ass=numpy.loadtxt('./Data_plaink1000/Assignments.h5'), 1, 100)
import io
test=proj.get_random_confs_from_states(ass=io.load_from('./Data_plaink1000/Assignments.h5'), 1, 100)
ass=io.open('./Data_plaink1000/Assignments.h5')
test=proj.get_random_confs_from_states(ass, 1, 100)
test=proj.get_random_confs_from_states(ass['arr_0], 1, 100)
test=proj.get_random_confs_from_states(ass['arr_0'], 1, 100)
ass
from msmbuilder import io
test=proj.get_random_confs_from_states(ass['arr_0'], 1, 100)
ass=io.loadh('./Data_plaink1000/Assignments.h5')
test=proj.get_random_confs_from_states(ass, 1, 100)
test=proj.get_random_confs_from_states(ass['arr_0'], 1, 100)
test
test.shpe
test.shape
test[0].shape
test[0]['XYZList'].sha[e
'
test[0]['XYZList'].shape
import readline
readline.write_history_file('ass.py')
import ass
test=numpy.array([2,2,2])
import numpy
test=numpy.array([2,2,2])
ohandle=open('test.out')
ohandle=open('test.out', 'w')
ohandle.write('%s\t%s' % (test, test))
val=1
ohandle.write('%s\t%s' % (val, test))
import numpy
rmsd=numpy.loadtxt('rmsd.dat')
for x in rmsd:
    bla
    
count=0
count=0
for x in rmsd:
    if x < 1.0:
        improt pdb
        pdb.set_trace*(
        
for x in rmsd:
    if x < 1.0:
        import pdb
        pdb.set_trace()
        
for x in rmsd:
    bla
    
count=0
frames=[]
for x in rmsd:
    if x < 1.0:
        frames.append(count)
        
frames
frames=[]
for x in rmsd:
    if x < 1.0:
        frames.append(count)
    count+=1
    
frames
frames[:10]
from msmbuilder import metrics
from msmbuilder.metrics import HYbrid
from msmbuilder.metrics import Hybrid
from metrics import hybrid
from msmbuilder import Trajectory
ref=Trajectory.load_from_pdb('fkbpLG2_xtal.pdb')
import numpy
prot=numpy.loadtxt('AtomIndices-min.dat', dtype=int)
lig=numpy.loadtxt('AtomIndices-ligand.dat', dtype=int)
ref
ref.keys()
ref['AtomID']
len(ref['AtomID'])
ref['XYZList'].shape
prot.shape
lig.shape
both=numpy.vstack((prot, lig))
both=numpy.hstack((prot, lig))
both.shape
new=Trajectory.load_from_pdb('fkbpLG2_xtal.pdb')
new=Trajectory.load_from_xtc('./Data/Gens.lh5')
new=Trajectory.load_from_lhdf('./Data/Gens.lh5')
new['XYZList'].shape
new['AtomID']
len(*new['AtomID'])
len(new['AtomID'])
ref['XYZList'].shape
from msmbuilder import Trajectory
gens=Trajectory.load_from_lhdf('./Data/Gens.lh5')
gens.save_to_xtc('./Data/Gens.xtc')
from msmbuilder import Trajectory
ef=Trajectory.load_from_lhdf('./Data/Gens.lh5')
ef.save_to_xtc('./Data/Gens.xtc')
import numpy
file=numpy.loadtxt('../LG2-new/all_rmsd')
max(file)
min(file)
len(file)
more ../LG2-new/get_all.sh
import numpy
bb=numpy.loadtxt('/home/mlawrenz/FKBP-FAH-new/LG2/trj-rmsd/all_bbrmsd')
lig=numpy.loadtxt('/home/mlawrenz/FKBP-FAH-new/LG2/trj-rmsd/all_rmsd')
len(lig(
)
bal
len(lig)
len(bb)
max(bb)
max(lig)
max(bb)/max(lig)
0.1*5+0.9*50
0.5*5+0.5*50
import numpy, os
data=numpy.loadtxt('./top_2percent.txt', dtype='str')
data[0]
for x in data:
    bla
    
ls ../../agonist/
for x in data:
    if 'ZINC' in data:
        os.system('../../decoys/agonist-decoys/%s .' % x)
    else:
        os.system('../../agonist/%s .' % x)
        
import readline
readline.write_history_file('copy.py')
import numpy
import numpy
improt os
import os
data=numpy.loadtxt('./top_2percent.txt', dtype='str')
for x in data:
    if 'ZINC' in data:
        os.system('cp ../../decoys/agonist-decoys/%s.mol2 .' % x)
    else:
        os.system('cp ../../agonist/%s.mol2 .' % x)
        
from msmbuilder import Trajectory
ref=Trajectory.load_from_lhdf('./Data/Gens.lh5')
ref.save_to_xtc('./Data/Gens.xtc')
import test
ls ../../decoys/agonist-decoys/
import test
import test
ls ../../agonist/*
print x
x=../../agonist/stereo-44216072-0.mol2
x='../../agonist/stereo-44216072-0.mol2'
x.split('-')
x.split('-')[-1]
x.split('-')[-1].split('.mol2')
x.split('.')
x.split('.')
x.split('-')[-1].split('.mol2')
num=x.split('-')[-1].split('.mol2')[0]
x.split(num)
import test
print num
from msmbuilder import Trajectory
g=Trajectory.load_from_lhdf('./Data_0.5/Gens.lh5')
g.save_to_xtc('./Data_0.5/Gens.xtc')
from msmbuilder import
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Data/Gens.lh5')
t.save_to_xtc('./Data/Gens.xtc')
from msmbuilder import Trajectory
re=Trajectory.load_from_lhdf('./Data_k100/Gens.lh5')
re.save_to_xtc('./Data_k100/Gens.xtc')
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Data_wt_k100/Gens.lh5')
t.save_to_xtc('./Data_wt_k100/Gens.xtc')
from msmbuilder import io
ass=io.loadh('./Data_wt_k100/Assignments.h5')
ass.shape
ass['arr_0'].shape
ls Trajectories/* | wc -l
for x in ass['arr_0']:
    for y in ass['arr_0'].shape[1]:
        bla
        
for x in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][x]==4)[0]
    if frames.size:
        print x, frames
        
import numpy
for x in range(0, ass['arr_0'].shape[0]):
    frames=numpy.where(ass['arr_0'][x]==4)[0]
    if frames.size:
        print x, frames
        
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Data_wt_k100/Gens.lh5')
t.save_to_xtc('./Data_wt_k100/Gens.xtc')
t=Trajectory.load_from_lhdf('./Data_eq_k100/Gens.lh5')
t.save_to_xtc('./Data_eq_k100/Gens.xtc')
from msmbuilder import Trajectory
g=Trajectory.load_from_lhdf('./Data_wt_k100/Gens.lh5')
g.save_to_xtc('./Data_wt_k100/Gens.xtc')
g=Trajectory.load_from_lhdf('./Data_eq_k100/Gens.lh5')
g.save_to_xtc('./Data_eq_k100/Gens.xtc')
from msmbuilder import Trajectory
g=Trajectory.load_from_lhdf('./Data_eq_k100/Gens.lh5')
g['AtomID']
g['AtomID'].shape
from scipy.io import *
from msmbuilder import MSMLib
from msmbuilder import io
 Assignments=io.loadh('bi_new_gen//Assignments.h5')
Assignments=io.loadh('bi_new_gen//Assignments.h5')
Counts = MSMLib.get_count_matrix_from_assignments(Assignments, lag_time=15, sliding_window=True)
Counts = MSMLib.get_count_matrix_from_assignments(Assignments['arr_0'], lag_time=15, sliding_window=True)
Counts = MSMLib.get_count_matrix_from_assignments(Assignments['Data'], lag_time=15, sliding_window=True)
scipy.io.mmwrite('%s/%s-%s' % (dir, str(tProb), 1), Counts)
from scipy.io import *
import scipy
scipy.io.mmwrite('%s/%s-%s' % (dir, str(tProb), 1), Counts)
iteration=1
scipy.io.mmwrite('bi_new_gen/tProb-%s/myx % iteration, Counts)
scipy.io.mmwrite('bi_new_gen/tProb-%s/mtx' % iteration, Counts)
scipy.io.mmwrite('bi_new_gen/tProb-%s.' % iteration, Counts)
test=range(0,6)
import pylab, numpy
pylab.plot(test)
pylab.show()
import glob
molecule=44216210
active=glob.glob('./active-examples/bound*%s*pdb' % molecule)
inactive=glob.glob('./inactive-examples/bound*%s*pdb' % molecule)
active
active[0].split('bound-')
active[0].split('bound-')[1]
active[0].split('bound-')[1].split('-')
active[0].split('bound-')[1].split('-')[:1]
active[0].split('bound-')[1].split('-')[:2]
'-'.join(active[0].split('bound-')[1].split('-')[:2])
newactive=[]
newinactive=[]
for file in active:
    newactive.append( '-'.join(file.split('bound-')[1].split('-')[:2])
    
    )
    
newactive
for file in inactive:
    newinactive.append( '-'.join(file.split('bound-')[1].split('-')[:2]))
    
newinactive
numpy.intersect1d(active,inactive)
import numpy
numpy.intersect1d(active,inactive)
for i in newactive:
    if i in newinactive:
        print i
        
newinactive
newactive
newactive=numpy.array(newactive)
newinactive=numpy.array(newinactive)
numpy.intersect1d(active,inactive)
numpy.where(newinactive=='bi-1782')
numpy.where(newactive=='bi-1782')
import numpy, pylab
molecules=[44216210, 44209282]
types=['inactive', 'active']
data=dict()
for m in molecules:
    data[m]=dict()
    data[m][type]=numpy.loadtxt('%s-%s-d192_k305.mindist.dat' % (m, type))
    
for m in molecules:
    data[m]=dict()
    for type in types:
        data[m][type]=numpy.loadtxt('%s-%s-d192_k305.mindist.dat' % (m, type))
        
data
for m in molecules:
    pylab.figure()
    pylab.hist(data[m]['active'], color='r', alpha=0.6)
    pylab.hist(data[m]['inactive'], color='b', alpha=0.6)
    
import readline
readline.write_history_file('plot_salt.py')
from msmbuilder import io
cd LG2-add/
project=io.loadh('ProjectInfo.yaml')
from msmbuilder import Project
project=Project.load_from('ProjectInfo.yaml')
project.n_trajs
project.traj_lengths
sum(project.traj_lengths)
ls *mdp
ls
ls ../*mdp
ls
from msmbuilder import Project
project=Project.load_from('ProjectInfo.yaml')
project.n_trajs
from msmbuilder import Project
project=Project.load_from('ProjectInfo.yaml')
project.n_trajs
from msmbuilder import Project
pwd
cd LG6-NEW
cd LG6-new/
from msmbuilder import Project
project=Project.load_from('ProjectInfo.yaml')
project.n_trajs
project.traj_lengths
min(project.traj_lengths)
max(project.traj_lengths)
sum(project.traj_lengths)
sum(project.traj_lengths)*0.1
sum(project.traj_lengths)*0.1*0.1
sum(project.traj_lengths)*0.1
pwd
cd ../LG6-add/
project=Project.load_from('ProjectInfo.yaml')
sum(project.traj_lengths)*0.1
198973.70000000001
project.n_trajs
min(project.traj_lengths)
project.traj_lengths
numpy.where(project.traj_lengths< 100)
import numpy
numpy.where(project.traj_lengths< 100)
frames=numpy.where(project.traj_lengths< 100)
frames=numpy.where(project.traj_lengths< 100)[0]
len(frames)
project.traj_lengths[frames]
len(project.traj_lengths)
import numpy
data=dict()
data['eq']=dict()
data['wt']=dict()
import numpy
data=dict()
data['eq']=dict()
data['wt']=dict()
data['wt']['bb']=numpy.loadtxt('Data_wt_k100/Gens.bb.dat')
data['eq']['bb']=numpy.loadtxt('Data_eq_k100/Gens.bb.dat')
data['eq']['lig']=numpy.loadtxt('Data_eq_k100/Gens.ligand.dat')
data['wt']['lig']=numpy.loadtxt('Data_wt_k100/Gens.ligand.dat')
pylab.hist(data['wt']['lig'])
import pylab
pylab.hist(data['wt']['lig'])
pylab.hist(data['wt']['bb'])
pylab.show()
pylab.hist(data['wt']['lig'], color='r')
pylab.hist(data['eq']['lig'], color='r')
pylab.show()
pylab.figure()
pylab.hist(data['eq']['lig'], color='r', alpha=0.6)
pylab.hist(data['wt']['lig'], color='b', alpha=0.6)
pylab.figure()
pylab.hist(data['eq']['bb'], color='r', alpha=0.6)
pylab.hist(data['wt']['bb'], color='b', alpha=0.6)
pylab.show()
import numpy, pylab
from msmbuilder import io
ass=io.loadh('./Assignments.h5)
ass=io.loadh('./Assignments.h5')
gens=numpy.loadtxt('Gens.ligand.dat')
numpy.where(gens>20)
frames=numpy.where(gens>20)[0]
len(gens)
len(frames)
bb=numpy.loadtxt('Gens.bb.dat')
bb[frames]
ass
ass['arr_0']
ass['arr_0'].shape
newass=ass
frames
frames=numpy.where(gens<20)[0]
frames
for i in range(0, ass['arr_0'].shape[0]):
    for j in ass['arr_0'][i]:
        if j in frames:
            pass
        else: 
        	newass[i][j]=3
        
newass
for i in range(0, ass['arr_0'].shape[0]):
    for j in ass['arr_0'][i]:
        if j in frames:
            pass
        else:
            newass['arr_0'][i][j]=3
            
import os
os.mkdir test/
os.mkdir('./test/')
from scipy.io import *
from scipy.sparse import lil_matrix
from scipy.sparse import lil_matrix
T=mmread('tCounts.mtx')
t.shape
T.shape
T
scipy.io.mmwrite('test/tProb-%s' % (dir, iteration), T)
iteration=4
importscipy
import scipy
scipy.io.mmwrite('test/tProb-%s' % (dir, iteration), T)
scipy.io.mmwrite('%s/tProb-%s' % (dir, iteration), T)
dir='./test/'
scipy.io.mmwrite('%s/tProb-%s' % (dir, iteration), T)
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Trajectories/trj1.lh5')
t.save_to_xtc('./trj1.xtc')
t=Trajectory.load_from_lhdf('./Trajectories/trj2.lh5')
t.save_to_xtc('./trj2.xtc')
ls
ls ../*pdb
from msmbuilder import Trajectory
ref=Trajectory.load_from_pdb('1dt7_A.pdb')
traj=Trajectory.load_from_lhdf('./Trajectories/trj0.lh5')
traj['XYZList']
traj['XYZList'].shape
traj['ResidueID']
set(traj['ResidueID'])
import glob
trajs=glob.glob('./run*/*xtc')
rmsd=glob.glob('../trj-rmsd/*.rmsd.dat')
rmsd[0]
len(trajs)
len(rmsd)
newrmsd=[]
newtrajs=[]
trajs[0]
for traj in sorted(trajs):
    newtrajs.append(traj.split('/')[-1].split('xtc')[0]
    
    )
    
len*newtrajs)
len(newtrajs)
newtrajs[0]
newtrajs=[]
for traj in sorted(trajs):
    newtrajs.append(traj.split('/')[-1].split('.xtc')[0])
    
for r in rmsd:
    bla
    
rmsd[0]
for r in rmsd:
    bla
    
for r in sorted(rmsd):
    newrmsd.append(r.split('/')[-1].split('.rmsd')[0]
    
    )
    
newrmsd[0]
trajs[0]
newtraj[0]
newtrajs[0]
len(newtrajs)
len(newrmsd)
set(newtrajs)
len(set(newtrajs))
duplicate=dict()
for traj in trajs:
    bla
    
trajs[0]
for traj in sorted(trajs):
    run=traj.split('/aln')[0])
    file=traj.split('/')[-1].split('.xtc')
    if file not in duplicate.keys():
        duplicate[file]=[]
    duplicate[file].append(run)
    
import readline
readline.write_history_file('run.py')
from msmbuilder import io
import numpy
systems=['LG2', 'LG3', 'LG6', 'LG9']
for sys in systems:
    bla
    
data=dict()
for sys in systems:
    data[sys]=numpy.loadtxt('./%s-add/trj-rmsd/all_rmsd.dat' % sys)
    
pylab.plot(data['LG2'])
import pylab
pylab.plot(data['LG2'])
pylab.show()
pylab.plot(data['LG2'], 'o')
pylab.show()
pylab.plot(data['LG3'], 'o')
pylab.show()
import pylab. numpy
import pylab, numpy
bbdata=dict()
ligdata=dict()
systems=['LG2', 'LG3', 'LG6', 'LG9']
for sys in systems:
    ligdata[sys]=numpy.loadtxt('./%s-add/trj-rmsd/all_rmsd.dat' % sys)
    bbdata[sys]=numpy.loadtxt('./%s-add/trj-rmsd/all_bbrmsd.dat' % sys)
    
for sys in systems:
    pylab.figure()
    pylab.plot(ligdata[sys], 'o')
    
pylab.show()
for sys in systems:
    pylab.figure()
    pylab.plot(bbdata[sys], 'o')
    pylab.title(sys)
    
pylab.show()
import pandas, pylab
names=pandas.read_csv('./all_bbrmsd.csv')
import pandas, pylab
test=pandas.read_csv('test.csv')
test
test.head()
test.yx[0:1]
test.ix[0:1]
test.ix[0]
test.ix[0]
test
test.ix[0, clone0]
test.ix[0, 'clone0']
test
test[1] <-0.5
test[1] <=0.5
test
test
test
names=pandas.read_csv('./all_bbrmsd.csv', names=['clone0', 'clone1', 'clone2'], header=0)
names
test
names=pandas.read_csv('./all_bbrmsd.csv', names=['index', 'clone0', 'clone1', 'clone2'], header=0)
names=pandas.read_csv('./test.csv', names=['clone0', 'clone1', 'clone2'], header=0)
names=pandas.read_csv('./test.csv', names=['index', 'clone0', 'clone1', 'clone2'], header=0)
head test.csv
mroe test.csv
more test.csv
names=pandas.read_csv('./test.csv', names=['clone0', 'clone1', 'clone2'], header=0)
names=pandas.read_csv('./test.csv', names=['clone0', 'clone1', 'clone2'], header=0, indexcol=0)
names=pandas.read_csv('./test.csv', names=['clone0', 'clone1', 'clone2'], header=0, index_col=0)
names=pandas.read_csv('./test.csv', names=['index', 'clone0', 'clone1', 'clone2'], header=0, index_col=0)
data=dict()
data['clone0']=numpy.array([0.5]*10)
import numpy
data['clone0']=numpy.array([0.5]*10)
data['clone0']=numpy.array([0.6, 0.5, 0.5, 0.5, 0.8, 1.1, 1.2, 1.5, 1.4, 1.2])
data['clone1']=numpy.array([0.6, 0.5, 0.5, 0.5, 0.8, 1.1, 1.2, 1.5, 1.4, 1.2])
data['clone0']=numpy.array([0.5]*10)
data['clone2']=numpy.array([0.7]*10)
s=pandas.Series(data, index=['clone0', 'clone1', clone2])
s=pandas.Series(data, index=['clone0', 'clone1', 'clone2'])
print s
s.index
s=pandas.Series(data)
s.index
s['clone0']
s[s >1.0]
s[s > 1.0]
print s
s0=pandas.Series(data[0])
data[0]
s=pandas.Series(data, index=['clone0', 'clone1', 'clone2'])
print s
data[0]
s=pandas.Series(data, index=['clone0', 'clone1', 'clone2'])
s0=pandas.Series(data['clone0'])
s0
s0[s0 == 0.5]
s0[s0 == 0.5]
s0=pandas.Series(data['clone1'])
s0[s0 > 1.0]
s0[[2,3]]
numpy.exp(s0)
s0=pandas.Series(data['clone1'], in)
df=pandase.DataFrame(data)
df=pandas.DataFrame(data)
df
df.index
df.columns
numbers=numpy.zeros((3,3))
numbers[1,2]=0.8
numbers[3,2]=0.8
numbers[2,1]=0.5
numbers[0,2]=0.8
pandas.DataFrame(numbers)
pandas.DataFrame(numbers, index=['y282', 't71', d34'])
pandas.DataFrame(numbers, index=['y282', 't71', 'd34'])
pandas.DataFrame(numbers, index=['y282', 't71', 'd34'], columns=['clone0', 'clone1', 'clone2'])
pandas.DataFrame(numbers, index=['y282', 't71', 'd34'], columns=['clone0', 'clone1', 'clone2'])
pandas.DataFrame.from_items('clone0',[0.5, 0.4], 'clone1', [0.6])
pandas.DataFrame.from_items('clone0',[0.5, 0.4], 'clone1', [0.6], orient=index, columns=['test1', 'test2'])
pandas.DataFrame.from_items('clone0',[0.5, 0.4], 'clone1', [0.6], orient='index', columns=['test1', 'test2'])
pandas.DataFrame.from_items(['clone0',[0.5, 0.4], 'clone1', [0.6]])
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6])])
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])])
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])], orient='column')
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])], orient='columns')
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])], orient='index')
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])], orient='index')
pandas.DataFrame.from_items([('clone0',[0.5, 0.4]), ('clone1', [0.6, 0.5])], orient='columns')
df
df[0]
df['clone0]
df['clone0']
df.columns
df.columns[0]
df[df.columns[0]]
df[df.columns[0]] >0.1
df > 0.1
df > 1.0
df['unbind']=df > 1.0
df['unbind']= df > 1.0
df['unbind']= (df > 1.0)
df
df > 1.0
new=df > 1.0
new
new[0]
new['clone0', 0]
new['clone0']
new['clone0'][0]
new['clone0',0]
new['clone0'][0]
df['clone0_trun']=df['clone0'][:5]
df
df.insert(1, df['clone2'])
df.insert(1, 'clone2_copy', df['clone2'])
df
df[2:3]
df.loc(2)
df.loc(2)
import glob
files=glob.glob('./*bbrmsd.dat')
import pandas
import numpy
from pandas import *
df=DataFrame(index=sorted(files))
df
for file in sorted(files):
    df[file]=numpy.loadtxt(file)
    
df
df.index
df.pop['all_bbrmsd.dat']
df['./aln-clone0-10.bbrmsd.dat]
df['./aln-clone0-10.bbrmsd.dat']
df['aln-clone0-10.bbrmsd.dat']
df.index
df.index[1]
df[df.index[1]]
df[0]
df[1]
df
index
df=DataFrame(index=sorted(files))
sorted(files)
sorted(files)[1:]
df=DataFrame(index=sorted(files)[1:])
for file in sorted(files)[1:]:
    df[file]=numpy.loadtxt(file)
    print file
    
print file
df[file]
df[0]
df=DataFrame(columns=sorted(files)[1:])
df[file]
for file in sorted(files)[1:]:
    df[file]=numpy.loadtxt(file)
    
df=DataFrame(columns=sorted(files)[1:])
df[df.index[0]]
df.index
df[df.columns[0]]
df.columns[0]
df[df.columns[0]]=Series(numpy.loadtxt(df.columns[0]))
 df[df.columns[0]]
df[df.columns[0]]
df.columns[0]
test=eries(numpy.loadtxt(df.columns[0]))
test=Series(numpy.loadtxt(df.columns[0]))
test
df[df.columns[0]]=test
df[df.columns[0]]
test[0]
test[:]
test
test=numpy.loadtxt(df.columns[0])
test
df[df.columns[0]]=test
stest=Series(test, column=0)
stest=Series(test, 'column'=0)
stest=Series(test, column=0)
stest=Series(test)
df[df.columns[0]]=stest
df[df.columns[0]]
df=DataFrame(numpy.zeros((len(files), 20000))columns=sorted(files)[1:])
df=DataFrame(numpy.zeros((len(files), 20000)), columns=sorted(files)[1:])
data=dict()
for file in files:
    data[file]=numpy.loadtxt(file)
    
df=DataFrame(data)
max(len(data[key]) for key in data]
max(len(data[key]) for key in data])
max(len(data[key]) for key in data.keys()])
max([len(data[key]) for key in data.keys()])
data.keys()
print files
data['./all_bbrmsd.dat']
data.pop('./all_bbrmsd.dat')
max([len(data[key]) for key in data.keys()])
data['./aln-clone353-6.bbrmsd.dat']
len(data['./aln-clone353-6.bbrmsd.dat'])
maxium=max([len(data[key]) for key in data.keys()])
test=numpy.hstack(data['./aln-clone353-6.bbrmsd.dat'], numpy.nan((maxim-311)))
test=numpy.hstack(data['./aln-clone353-6.bbrmsd.dat'], numpy.nan((maxium-311)))
test=numpy.hstack(data['./aln-clone353-6.bbrmsd.dat'], numpy.nan((maxium-311)))
filler=numpy.empty((maxium-311))
filler=numpy.nan
filler
filler=numpy.empty((maxium-311))
filler[:]=numpy.nan
filler
test=numpy.hstack(data['./aln-clone353-6.bbrmsd.dat'], filler)
test=numpy.hstack((data['./aln-clone353-6.bbrmsd.dat'], filler))
len(test)
maium
maxium
for file in files:
    filler=numpy.empty((maxium-len(data[file])))
    
newdaa=dict()
newdata=dict()
files.pop('./all_bbrmsd.dat')
numpy.where(files==0)
files.index('./all_bbrmsd.dat')
files.pop(948)
for file in files:
    filler=numpy.empty((maxium-len(data[file])))
    newdata[file]=numpy.hstack((data[file], filler))
    
newdata
df=DataFrame(newdata)
df
df.index
df.columns
df.columns[0]
df[df.columns[0]]
newdata=dict()
for file in files:
    filler=numpy.empty((maxium-len(data[file])))
    filler[:]=numpy.nan
    newdata[file]=numpy.hstack((data[file], filler))
    
df=DataFrame(newdata)
df[df.columns[0]]
df.loc['./aln-clone0-10.bbrmsd.dat']
df > 10
df.max()
import readline
readline.write_history_files('load_df.py')
import readline
readline.write_history_file('load_df.py')
from msmbuilder import io
import numpy
from msmbuilder import io
import pandas
from msmbuilder import io
import glob
import numpy
trajs=numpy.loadtxt('../trj/traj_list.txt', dtype=str)
rmsds=glob.glob('*bbrmsd.dat')
rmsds[0]
trajs[0]
new=[i.split('/')[-1].split('.xtc')[0] for i in trajs]
new
for file in rmsds:
    if file.split('.bb')[0] not in new:
        print file
        
for file in rmsds:
    test=file.split('aln-')[1].split('.bb')[0]
    if 'aln-%s' % test not in new:
        print file
        
for file in rmsds:
    test=file.split('ln-')[1].split('.bb')[0]
    if 'aln-%s' % test not in new:
        print file
        
for file in rmsds:
    test=file.split('ln-')[1].split('.bb')[0]
    print 'aln-%s' % test
    
for file in rmsds:
    print file
    print 'aln-%s' % test
    
for file in rmsds:
    print file
    test=file.split('ln-')[1].split('.bb')[0]
    print test
    
import readline
readline.write_history_file('remove.py')
from msmbuilder import PLAINRMSD
from msmbuilder import lprmsd
from lprmsd import LPRMSD
import LPRMSD
from LPRMSD import lprmsd
import plainrmsd
from msmbuilder import baseclasses
from msmbuilder import Trajectory
t=Trajectory.load_from_lhdf('./Trajectories/trj1.lh5')
t.save_to_xtc('trj1.xtc')
from msmbuilder import Trajectory
from msmbuilder import io
rmsd=io.loadh('./Data-rmsd/Assignments.h5.distances')
dist=io.loadh('./Data-dist/Assignments.h5.distances')
dist.shape
dist['arr_0'].shape
dist['arr_0'][0]
max(dist['arr_0'][0])
max(dist['arr_0'][1])
from msmbuilder import io, Trajectory
g=Trajectory.load_from_lhdf('./Data-rmsd/Gens.lh5')
g.save_to_xtc('./Data-rmsd/Gens.xtc')
g=Trajectory.load_from_lhdf('./Data-combo/Gens.lh5')
g.save_to_xtc('./Data-combo/Gens.xtc')
data=dict()
import numpy, pylab
data['rmsd']=numpy.loadtxt('./Data-combo/Gens.rmsd.dat')
data['dist']=numpy.loadtxt)('./Data-combo/Gens.rLG2n.rPHECA.dat')
data['dist']=numpy.loadtxt('./Data-combo/Gens.rLG2n.rPHECA.dat')
data
combodata=data.copy()
rmsddata=dict()
rmsddata['rmsd']=numpy.loadtxt('./Data-rmsd/Gens.rmsd.dat')
rmsddata['dist']=numpy.loadtxt('./Data-rmsd/Gens.rLG2n.rPHECA.dat'
)
pylab.scatter(rmsddata['dist'], rmsddata['rmsd'], color='r')
pylab.scatter(combodata['dist'], combodata['rmsd'], color='b')
pylab.show()
rmsdmap=numpy.loadtxt)
'-'.join(active[0].split('bound-')[1].split('-')[:2])
rmsdmap=numpy.loadtxt('./Data-rmsd/MSMl10/Mapping.dat', dtype=int)
combomap=numpy.loadtxt('./Data-combo/MSMl10/Mapping.dat', dtype=int)
numpy.where(combomap==-1)
numpy.where(combomap==-1)[0]
len(numpy.where(combomap==-1)[0])
len(numpy.where(rmsdmap==-1)[0])
len(numpy.where(combomap==-1)[0])
combodata['dist']
trimcombo=numpy.where(combomap==-1)[0]
rmsdcombo=numpy.where(rmsdmap==-1)[0]
trimrmsd=numpy.where(rmsdmap==-1)[0]
rmsddata['dist'][trimrmsd]
combodata['dist'][trimcombo]
rmsddata['dist'][trimcombo]
rmsddata['dist'][trimrmsd]
import readline
readline.write_history_file('format.py')
from msmbuilder import Project
proj=Project.load_from('ProjectInfo.yaml')
proj.traj_lengths
proj.traj_lengths
from msmbuilder import io
orig=io.loadh('./d6/Data/Assignments.h5')
coarse=io.loadh('./d12/Data/Assignments.h5')
orig['arr_0'].shape
orig['arr_0'][0]
orig_rmsd=numpy.loadtxt('./d6/Gens.rmsd.dat')
import numpy
orig_rmsd=numpy.loadtxt('./d6/Gens.rmsd.dat')
coarse_rmsd=numpy.loadtxt('./d12/Gens.rmsd.dat')
dist=io.loadh('d6/Data/Assignments.h5.distances')
dist
dsit['arr_0'][0]
dist['arr_0'][0]
orig_dist=io.loadh('d6/Data/Assignments.h5.distances')
coarse_dist=io.loadh('d12/Data/Assignments.h5.distances')
max(orig_dist['arr_0'].flatten())
import readline
readline.write_history_file('coarse_ass.py')
import coarse_ass
orig.keys()
ls
import coarse_ass
orig, coarse=coarse_ass()
orig, coarse=coarse_ass.main()
orig
import coarse_ass
coarse, orig=coarse_ass.main()
len(coarse['rmsd'])
for i in orig['ass'][0]:
    print rmsd[i]
    
for i in orig['ass'][0]:
    print orig['rmsd'][i]
    
for i in orig['ass'][0]:
    if orig['rmsd'][i] > 20.0:
        print i
        
for i in orig['ass'][0]:
for (n,i) in enumerate(orig['ass'][0]):
    if orig['rmsd'][i] > 20.0:
        state=coarse['ass'][n]
        print coarse['rmsd'][state]
        
coarse['ass'].shape
for (n,i) in enumerate(orig['ass'][0]):
    if orig['rmsd'][i] > 20.0:
        state=coarse['ass'][0][n]
        print coarse['rmsd'][state]
        
for (n,i) in enumerate(orig['ass'][0]):
    if orig['rmsd'][i] > 20.0:
        state=coarse['ass'][0][n]
        print coarse['dist'][0][n]
        
import readline
readline.write_history_file('test.py')
