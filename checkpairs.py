from msmbuilder import metrics
import lprmsd
import pylab
import numpy
from msmbuilder import Trajectory
from msmbuilder import Project

#residue_pairs=numpy.loadtxt('./residuepairs_select_0based.dat', dtype=int)

residue_pairs=numpy.loadtxt('./residuepairs_0based_select_structure1.dat', dtype=int)
pairmetric2=metrics.ContinuousContact(contacts=residue_pairs, scheme='closest-heavy')
#pairmetric min 0.00255089723477 max 10.8130541983

residue_pairs=numpy.loadtxt('./residuepairs_select_0based.dat', dtype=int)
pairmetric1=metrics.ContinuousContact(contacts=residue_pairs, scheme='closest-heavy')
#pairmetric1 min 0.00447584821743 max 22.0643260692 

#pairmetric=metrics.BooleanContact(metric='matching', contacts=residue_pairs, cutoff=0.5, scheme='closest-heavy')
ligand_inds = numpy.loadtxt('./p53_Calpha_indices.dat', dtype='int')
prot_inds = numpy.loadtxt('./s100b_Calpha_indices.dat', dtype='int')

rmsdmetric2=metrics.RMSD(ligand_inds)
rmsdmetric1 = lprmsd.LPRMSD(atomindices=prot_inds, permuteindices=None, altindices=ligand_inds)
pdb=Trajectory.load_from_pdb('structure_1.pdb')
#project = Project.load_from('/home/kkappel/p53/s100b_md/msm_80ns_p53CA/ProjectInfo.yaml')
project = Project.load_from('ProjectInfo.yaml')

hybrid1 = metrics.Hybrid([rmsdmetric2, pairmetric1], weights=[1/1.16, 1/22.1])
hybrid2 = metrics.Hybrid([rmsdmetric2, pairmetric2], weights=[1/1.16, 1/10.8])


data=dict()
# check pair metric for traj0
metrics=[pairmetric1, pairmetric2, hybrid1, hybrid2 ]
names=['pair-selectall', 'pair-selectstructure1', 'hybrid-selectall', 'hybrid-selectstructure1']
for (metric, name) in zip(metrics, names):
    data[name]=dict()
    data[name]['min']=1000000
    data[name]['max']=0
    pylab.figure()
    for n in range(0, project.n_trajs):
        print "on %s" % n
        traj=project.load_traj(n)
        ptraj=rmsdmetric1.prepare_trajectory(traj)
        ppdb=rmsdmetric1.prepare_trajectory(pdb)
        lprmsd, xout=rmsdmetric1.one_to_all_aligned(ppdb, ptraj, 0)
        ppdb=metric.prepare_trajectory(pdb)
        ptraj=metric.prepare_trajectory(traj)
        output=metric.one_to_all(ppdb, ptraj, 0)
        if max(output) > data[name]['max']:
            data[name]['max']=max(output)
        if min(output) < data[name]['min']:
            data[name]['min']=min(output)
        pylab.scatter([i*10 for i in lprmsd], output)
        #numpy.savetxt('./trj-rmsd/trj%s_lprmsd.dat' % n, lprmsd)
        #numpy.savetxt('./trj-structure1-pairs/trj%s_pairs.dat' % n, output)
        pylab.hold(True)
    pylab.title(name)
    pylab.savefig('checkdist-%s-select.png' % name, dpi=300)
    print "%smetric min %s max %s " % (name, data[name]['min'], data[name]['max'])

pylab.show()

