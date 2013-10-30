from msmbuilder import metrics
import os
import lprmsd
import pylab
import numpy
from msmbuilder import Trajectory
from msmbuilder import Project

#load in atom pairs here
#atom_pairs=numpy.loadtxt('./atompairs.dat', dtype=int)
#pairmetric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=atom_pairs)

# load in residue pairs here
#residue_pairs=numpy.loadtxt('./residuepairs_0based_select_structure1.dat', dtype=int)
#pairmetric2=metrics.ContinuousContact(contacts=residue_pairs, scheme='closest-heavy')

# for computing lprmsd (the standard for comparison)
ligand_inds = numpy.loadtxt('p53_Calpha_indices.dat', dtype='int')
prot_inds = numpy.loadtxt('sirtuin_Calpha_indices.dat', dtype='int')
rmsdmetric = lprmsd.LPRMSD(atomindices=prot_inds, permuteindices=None, altindices=ligand_inds)

pdb=Trajectory.load_from_pdb('sir2_bound_reference.pdb')
project = Project.load_from('../sirtuin_round1/ProjectInfo.yaml')

if not os.path.exists('./Trajectories-metric'):
    os.mkdir('./Trajectories-metric')
# loop over 3 metrics
#metrics=[pairmetric, pairmetric2, rmsdmetric ] 
#names=['atompairs', 'residuepairs', 'lprmsd']
metrics=[ rmsdmetric, ] 
names=[  'lprmsd',]
for (metric, name) in zip(metrics, names):
    for n in range(0, project.n_trajs):
        print "on traj %s" % n
        traj=project.load_traj(n)
        if name=='lprmsd':
            ppdb=metric.prepare_trajectory(pdb)
            ptraj=metric.prepare_trajectory(traj)
            # for lprmsd, measure RMSD vs. Bound State
            output=metric.one_to_all(ppdb, ptraj, 0)
            numpy.savetxt('./Trajectories-metric/trj%s_%s.dat' % (n, name), [i*10 for i in output])
        else:
            distances=metric.prepare_trajectory(traj)
            # Multiple distances btw atoms or residues, combine all
            all=[]
            for i in range(0, distances.shape[0]):
                all.append(numpy.sqrt(sum([j**2 for j in distances[i,:]])))
            numpy.savetxt('./Trajectories-metric/trj%s_%s.dat' % (n, name), all)


