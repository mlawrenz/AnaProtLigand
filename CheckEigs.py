#!/bin/python
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt
import glob
import random
from scipy.io import *
from msmbuilder import Conformation
from msmbuilder import MSMLib, msm_analysis
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

def get_structure(modeldir, eig, gen_states, states, gens, project, ass, type=None):
    # for assignments, take mapped states
    # for gens, take gen states
    for (state, gstate) in zip(states, gen_states):
        sample=project.empty_traj()
        (a, b, c) =gens['XYZList'].shape
        sample['XYZList']=numpy.zeros((1, b, c), dtype=numpy.float32)
        sample['XYZList'][0]=gens['XYZList'][gstate]
        sample.save_to_pdb('%s/eig-states/eig%s-%s-state%s-centroid.pdb' % (modeldir, eig, type, int(state)))
        t=project.get_random_confs_from_states(ass['arr_0'], [int(state),], 5)
        for j in range(0, 5):
            sample=project.empty_traj()
            (a, b, c) =gens['XYZList'].shape
            sample['XYZList']=numpy.zeros((1, b, c), dtype=numpy.float32)
            sample['XYZList'][0]=t[0]['XYZList'][j]
            sample.save_to_pdb('%s/eig-states/eig%s-%s-state%s-%s.pdb' % (modeldir, eig, type, int(state), j))



def main(modeldir, gensfile, write=False):
    if not os.path.exists('%s/eig-states/' % modeldir):
        os.mkdir('%s/eig-states/' % modeldir)
    ohandle=open('%s/eig-states/eiginfo.txt' % modeldir, 'w')
    project=Project.load_from('../sirtuin_round1/ProjectInfo.yaml')
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)

    pdb=Trajectory.load_from_pdb('sir2_bound_xtal.pdb')
    gens=Trajectory.load_from_xtc(gensfile, Conf=pdb)
    T=mmread('%s/tProb.mtx' % modeldir)
    data=dict()
    data['rmsd']=numpy.loadtxt('%s.rmsd.dat' % gensfile.split('.xtc')[0])
    com=numpy.loadtxt('%s.vmd_com.dat' % gensfile.split('.xtc')[0], usecols=(1,))
    data['helixrmsd']=numpy.loadtxt('%s.helixrmsd.dat' % gensfile.split('.xtc')[0])
    com=numpy.loadtxt('%s.vmd_com.dat' % gensfile.split('.xtc')[0], usecols=(1,))
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    map_rmsd=[]
    map_helixrmsd=[]
    map_com=[]
    for x in range(0, len(data['rmsd'])):
        if map[x]!=-1:
            map_com.append(data['com'][x])
            map_rmsd.append(data['rmsd'][x])
            map_helixrmsd.append(data['helixrmsd'][x])
    
    map_com=numpy.array(map_com)
    map_rmsd=numpy.array(map_rmsd)
    map_helixrmsd=numpy.array(map_helixrmsd)

    T=mmread('%s/tProb.mtx' % modeldir)
    eigs_m=msm_analysis.get_eigenvectors(T, 10)

    cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive

    print numpy.shape(eigs_m[1][:,1])
    for i in range(0,3):
        order=numpy.argsort(eigs_m[1][:,i])
        if i==0:
            maxes=[]
            gen_maxes=[]
            values=[]
            ohandle.write('eig%s maxes\n' % i)
            ohandle.write('state\tgenstate\tmagnitude\thelixrmsd\tcom\n')
            for n in order[::-1][:5]:
                gen_maxes.append(numpy.where(map==n)[0])
                maxes.append(n)
                values.append(eigs_m[1][n,i])
                ohandle.write('%s\t%s\t%s\t%s\t%s\n' % (n, numpy.where(map==n)[0], eigs_m[1][n,i], map_helixrmsd[n], map_com[n]))
            print "maxes at ",  maxes, values
            maxes=numpy.array(maxes)
            if write==True:
                get_structure(modeldir, i, gen_maxes, maxes, gens, project, ass, type='max')
        else:
            maxes=[]
            gen_maxes=[]
            values=[]
            ohandle.write('eig%s maxes\n' % i)
            for n in order[::-1][:5]:
                gen_maxes.append(numpy.where(map==n)[0])
                maxes.append(n)
                values.append(eigs_m[1][n,i])
                ohandle.write('%s\t%s\t%s\t%s\t%s\n' % (n, numpy.where(map==n)[0], eigs_m[1][n,i], map_helixrmsd[n], map_com[n]))
            print "maxes at ",  maxes, values
            order=numpy.argsort(eigs_m[1][:,i])
            mins=[]
            gen_mins=[]
            values=[]
            ohandle.write('eig%s mins\n' % i)
            for n in order[:5]:
                gen_mins.append(numpy.where(map==n)[0])
                mins.append(n)
                values.append(eigs_m[1][n,i])
                ohandle.write('%s\t%s\t%s\t%s\t%s\n' % (n, numpy.where(map==n)[0], eigs_m[1][n,i], map_helixrmsd[n], map_com[n]))
            print "mins at ",  mins, values
            if write==True:
                get_structure(modeldir, i, gen_maxes,  maxes, gens, project, ass, type='max')
                get_structure(modeldir, i, gen_mins,  mins, gens, project, ass, type='min')
        pylab.scatter(map_helixrmsd[order], map_rmsd[order], c=eigs_m[1][order,i], cmap=cm, s=1000*abs(eigs_m[1][order,i]), alpha=0.5)
        print map_com[order][numpy.argmax(eigs_m[1][order,i])]
        print eigs_m[1][order,i][1]
        CB=pylab.colorbar()
        l,b,w,h=pylab.gca().get_position().bounds
        ll, bb, ww, hh=CB.ax.get_position().bounds
        CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
        CB.set_label('Eig%s Magnitudes' % i)
        ylabel=pylab.ylabel('Ligand RMSD to Xtal ($\AA$)')
        xlabel=pylab.xlabel('Ligand Beta RMSD')
        #xlabel=pylab.xlabel(r'P Active Site - L COM Distance ($\AA$)')
        pylab.legend(loc=8, frameon=False)
        pylab.savefig('%s/2deigs%i_com_prmsd.png' %(modeldir, i),dpi=300)
        pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--modeldir', dest='modeldir',
                      help='msm directory')
    parser.add_option('-g', '--gensfile', dest='gensfile',
                      help='gens files')
    parser.add_option('-w', action="store_true", dest="write")
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.write==True:
        main(modeldir=options.modeldir, gensfile=options.gensfile, write=True)
    else:
        main(modeldir=options.modeldir, gensfile=options.gensfile)
