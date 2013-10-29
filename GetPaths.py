#!/bin/python
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt, metrics
from msmbuilder.geometry import dihedral as _dihedralcalc
import lprmsd
import glob
import random
from scipy.io import *
from msmbuilder import Conformation
from msmbuilder import MSMLib
import optparse
import numpy
import os
from numpy import linalg
import pylab


def get_dihedrals(dir, ind, traj):
    dihed=_dihedralcalc.compute_dihedrals(traj['XYZList'], ind, degrees=True)
    return numpy.array([i[0] for i in dihed])

def check_metric_pairs(dir, traj):
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

def check_rmsd(dir, pdb, ligand_inds, traj):
    prot_inds = numpy.loadtxt('%s/AtomIndices-ca.dat' % dir, dtype='int')
    rmsdmetric1 = lprmsd.LPRMSD(atomindices=prot_inds, permuteindices=None, altindices=ligand_inds)
    ptraj=rmsdmetric1.prepare_trajectory(traj)
    ppdb=rmsdmetric1.prepare_trajectory(pdb)
    output, xout=rmsdmetric1.one_to_all_aligned(ppdb, ptraj, 0)
    return numpy.array([i*10 for i in output])

def build_metric(dir, pdb, traj):
    data=dict()
    ligand_inds = numpy.loadtxt('%s/AtomIndices-ligand.dat' % dir, dtype=int)
    data['rmsd']=check_rmsd(dir, pdb, ligand_inds, traj)
    ligand_inds = numpy.loadtxt('/nobackup/mlawrenz/FKBP-FAH-results2/loop2.dat', dtype='int')
    data['loop2']=check_rmsd(dir, pdb, ligand_inds, traj)
    ligand_inds = numpy.loadtxt('/nobackup/mlawrenz/FKBP-FAH-results2/loop3.dat', dtype='int')
    data['loop3']=check_rmsd(dir, pdb, ligand_inds, traj)
    index1=numpy.zeros((1,4)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
    index1[0]=numpy.loadtxt('%s/omega_indices.txt' % dir, dtype=int, ndmin=1)
    data['omega']=get_dihedrals(dir, index1, traj)
    index2=numpy.zeros((1,4)) # here i just need to tell _dihedralcalc that we have one set of 4 coordinates
    index2[0]=numpy.loadtxt('%s/phi_indices.txt' % dir, dtype=int, ndmin=1)
    data['phi']=get_dihedrals(dir, index2, traj)
    pairs=check_metric_pairs(dir, traj)
    for key in pairs.keys():
        data[key]=[i*10 for i in pairs[key]]
    return data

def map_size(x):
    if x==0.5:
        size=500
    elif abs(x-0.5) < 0.1 and abs(x-0.5) > 0:
        size=300
    elif abs(x-0.5) < 0.2 and abs(x-0.5) > 0.1:
        size=200
    elif abs(x-0.5) < 0.3 and abs(x-0.5) > 0.2:
        size=100
    else:
        size=50
    return size

def main(modeldir, genfile,  type, write=False):
    data=dict()
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)
    frames=numpy.where(map!=-1)[0]

    unbound=numpy.loadtxt('%s/tpt-rmsd-%s/unbound_%s_states.txt' % (modeldir, type, type), dtype=int)
    bound=numpy.loadtxt('%s/tpt-rmsd-%s/bound_%s_states.txt' % (modeldir, type, type), dtype=int)

    dir=modeldir.split('Data')[0]
    name=glob.glob('%s/fkbp*xtal*pdb' % dir)
    pdb=Trajectory.load_from_pdb(name[0])
    paths=io.loadh('%s/tpt-rmsd-%s/Paths.h5' % (modeldir, type))

    committors=numpy.loadtxt('%s/commitor_states.txt' % modeldir, dtype=int)
    colors=['red', 'orange', 'green', 'cyan', 'blue', 'purple']
    colors=colors*40
    if type=='strict':
        ref=5
    elif type=='super-strict':
        ref=3
    elif type=='medium':
        ref=10
    elif type=='loose':
        ref=15
    #for p in range(0, 3):
    for p in range(0, 1):
        path=paths['Paths'][p]
        print "Bottleneck", paths['Bottlenecks'][p]
        flux=paths['fluxes'][p]/paths['fluxes'][0]
        if flux < 0.2:
            break
        print "flux %s" % flux
        frames=numpy.where(path!=-1)[0]
        path=numpy.array(path[frames], dtype=int)
        print path
        if write==True:
            size=(paths['fluxes'][p]/paths['fluxes'][0])*1000
            traj=Trajectory.load_from_xtc('%s/tpt-rmsd-%s/path%s_sample20.xtc' % (modeldir, type, p), Conf=pdb)
            data=build_metric(dir, pdb, traj)
            dir=modeldir.split('Data')[0]
            for op in sorted(data.keys()):
            #for op in residues:
                pylab.figure()
                pylab.scatter(data['rmsd'], data[op], c=colors[p], alpha=0.7) #, s=size)
                for j in paths['Bottlenecks'][p]:
                    frame=numpy.where(paths['Paths'][p]==j)[0]
                    pylab.scatter(data['rmsd'][frame*20], data[op][frame*20], marker='x', c='k', alpha=0.7, s=50)
                    location=numpy.where(committors==paths['Paths'][p][frame])[0]
                    if location.size:
                        print "path %s state %s bottleneck in committors" % (p, j)
                        print data['rmsd'][frame*20], data[op][frame*20]
                pylab.title('path %s' % p)
                pylab.xlabel('P-L RMSD')
                #pylab.xlabel('P-L COM')
                pylab.ylabel(op)
                pylab.xlim(0,max(data['rmsd'])+5)
                #pylab.ylim(0,max(data[op])+5)
                pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='genfile')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    parser.add_option('-w', action="store_true", dest="write")
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.write==True:
        main(modeldir=options.dir, genfile=options.genfile, type=options.type, write=True)
    else:
        main(modeldir=options.dir, genfile=options.genfile, type=options.type)

