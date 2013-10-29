#!/bin/python
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt
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

def map_size(x):
    if x==0.5:
        size=500
    elif abs(x-0.5) < 0.1 and abs(x-0.5) > 0:
        size=300
    elif abs(x-0.5) < 0.2 and abs(x-0.5) > 0.1:
        size=200
    elif abs(x-0.5) < 0.3 and abs(x-0.5) > 0.2:
        size=100
    elif x==0:
        size=10
    elif x==1:
        size=10
    else:
        size=50
    return size

def main(modeldir, genfile, type, write=False):
    proj=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)
    frames=numpy.where(map!=-1)[0]
    data=dict()
    data['rmsd']=numpy.loadtxt('%s.rmsd.dat' % genfile.split('.lh5')[0])
    data['rmsd']=data['rmsd'][frames]
    com=numpy.loadtxt('%s.vmd_com.dat' % genfile.split('.lh5')[0], usecols=(1,))
    refcom=com[0]
    data['com']=com[1:]
    data['com']=numpy.array(data['com'][frames])

    residues=['F36', 'H87', 'I56', 'I90', 'W59', 'Y82', 'hydrophob_dist', 'oxos_dist']
    loops=['loop1', 'loop2', 'loop3']
    for loop in loops:
        data[loop]=numpy.loadtxt('%s.%srmsd.dat' % (genfile.split('.lh5')[0], loop))
        data[loop]=data[loop][frames]
    for res in residues:
        file='%s_%spair.dat' % (genfile.split('.lh5')[0], res)
        if os.path.exists(file):
            data[res]=numpy.loadtxt(file)
            data[res]=data[res][frames]
    angles=['phi', 'omega']
    for ang in angles:
        file='%s_%s.dat' % (genfile.split('.lh5')[0], ang)
        if os.path.exists(file):
            data[ang]=numpy.loadtxt(file)
            data[ang]=data[ang][frames]
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    T=mmread('%s/tProb.mtx' % modeldir)
    unbound=numpy.loadtxt('%s/tpt-%s/unbound_%s_states.txt' % (modeldir, type, type), dtype=int)
    bound=numpy.loadtxt('%s/tpt-%s/bound_%s_states.txt' % (modeldir, type, type), dtype=int)

    Tdense=T.todense()
    Tdata=dict()
    for i in unbound:
        for j in unbound:
            if Tdense[i,j]!=0:
                if i not in Tdata.keys():
                    Tdata[i]=[]
                Tdata[i].append(j)
    #print Tdata
    cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive
    Q=tpt.calculate_committors(unbound, bound, T)
    ohandle=open('%s/commitor_states.txt' % modeldir, 'w')
    for i in range(0,len(Q)):
        if Q[i]>0.40 and Q[i]<0.6:
            ohandle.write('%s\n' % i)
            #t=project.get_random_confs_from_states(ass['arr_0'], [int(i),], 20)
            #t[0].save_to_xtc('%s/commottor_state%s.xtc' % (modeldir, i))
    if write==True:
        for op in sorted(data.keys()):
            pylab.figure()
            pylab.scatter(data['com'], data[op],  c=Q, cmap=cm, alpha=0.7, s=[map_size(i) for i in Q])
            pylab.xlabel('L RMSD')
            pylab.ylabel(op)
            pylab.colorbar()
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
        main(modeldir=options.dir,  genfile=options.genfile, type=options.type, write=True)
    else:
        main(modeldir=options.dir,  genfile=options.genfile, type=options.type)
