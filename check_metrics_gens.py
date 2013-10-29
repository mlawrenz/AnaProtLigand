import numpy, glob
import os
import optparse
from msmbuilder.geometry import dihedral as _dihedralcalc
from msmbuilder import Trajectory
from msmbuilder import metrics
import pylab


def check_metric_pairs(dir):
    atom_pairs=numpy.loadtxt('%s/atompairs.dat' % dir, dtype=int)
    names=numpy.loadtxt('%s/atompairs-map.txt' % dir, usecols=(2,),  dtype=str)
    indices=[i for i in atom_pairs]
    pairs=dict()
    for (ind, name) in zip(indices, names):
        index1=numpy.zeros((1,2)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
        index1[0]=ind
        pairmetric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=index1)
        pairs[name]=pairmetric.prepare_trajectory(traj)
    return pairs

def main(genfile):
    dir=os.path.dirname(genfile)
    traj=Trajectory.load_from_lhdf(genfile)
    rmsd=numpy.loadtxt('%s.rmsd.dat' % genfile.split('.lh5')[0])
    atom_pairs=numpy.loadtxt('./atompairs.dat', dtype=int)
    names=numpy.loadtxt('./atompairs-map.txt', usecols=(2,),  dtype=str)
    indices=[i for i in atom_pairs]
    for (ind, name) in zip(indices, names):
        index1=numpy.zeros((1,2)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
        index1[0]=ind
        pairmetric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=index1)
        distances=pairmetric.prepare_trajectory(traj)
        pylab.figure()
        pylab.scatter(rmsd, distances, label=name)
        pylab.legend()
        pylab.show()
        numpy.savetxt('%s_%spair.dat' % (genfile.split('.lh5')[0], name), distances)
    index1=numpy.zeros((1,4)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
    index1[0]=numpy.loadtxt('omega_indices.txt', dtype=int, ndmin=1)
    index2=numpy.zeros((1,4)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
    index2[0]=numpy.loadtxt('phi_indices.txt', dtype=int, ndmin=1)
    indices=[index1, index2]
    names=['omega', 'phi']
    for (ind, name) in zip(indices, names):
        dihed=_dihedralcalc.compute_dihedrals(traj['XYZList'], ind, degrees=True)
        dihed=[i[0] for i in dihed]
        numpy.savetxt('%s_%s.dat' % (genfile.split('.lh5')[0], name), dihed)
    atom_pairs=numpy.loadtxt('conformation_pairs.dat', dtype=int)
    names=['oxos_dist', 'hydrophob_dist']
    for (pair, name) in zip(atom_pairs, names):
        index=numpy.zeros((1,2), dtype=int)
        index[0]=pair
        metric= metrics.AtomPairs(metric='euclidean', p=1, atom_pairs=index)
        distances=metric.prepare_trajectory(traj)
        distances=[i[0]*10 for i in distances]
        pylab.figure()
        pylab.scatter(rmsd, distances, label=name)
        pylab.legend()
        pylab.show()
        numpy.savetxt('%s_pairs.dat' % genfile.split('.lh5')[0], distances)


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='input gens file xtc')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main( genfile=options.genfile)


