from msmbuilder import metrics
import multiprocessing
import os
import lprmsd
import pylab
import numpy
from msmbuilder import Trajectory
from msmbuilder import Project

ligand_inds = numpy.loadtxt('p53_Calpha_indices.dat', dtype='int')
prot_inds = numpy.loadtxt('s100b_Calpha_indices.dat', dtype='int')
rmsdmetric1 = lprmsd.LPRMSD(atomindices=prot_inds, permuteindices=None, altindices=ligand_inds)
pdb=Trajectory.load_from_pdb('structure_1.pdb')
project = Project.load_from('ProjectInfo.yaml')


data=dict()
if not os.path.exists('./Trajectories-metric/'):
    os.mkdir('./Trajectories-metric/')
for n in range(0, project.n_trajs):
    if os.path.exists('./Trajectories-metric/trj%s_lprmsd.dat' % n):
        pass
    else:
        print "on traj %s" % n
        traj=project.load_traj(n)
        ptraj=rmsdmetric1.prepare_trajectory(traj)
        ppdb=rmsdmetric1.prepare_trajectory(pdb)
        lprmsd, xout=rmsdmetric1.one_to_all_aligned(ppdb, ptraj, 0)
        numpy.savetxt('./Trajectories-metric/trj%s_lprmsd.dat' % n, [i*10 for i in lprmsd])

