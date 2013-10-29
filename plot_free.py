from msmbuilder import Trajectory, PMF3D
import multiprocessing
import optparse
import pylab
from numpy import * 
import glob
import os
import sys
import pickle
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from scipy import integrate


def main(dir):
    ref=loadtxt('%s/target.txt' % dir.split('Data')[0])
    upperstring='MSM'
    lowerstring=[x.lower() for x in upperstring.split()][0]
    if upperstring not in dir:
        if lowerstring in dir:
            if 'coarse' in dir:
                lag=dir.split('_coarse')[0].split('%sl' % lowerstring)[1]
            else:
                lag=dir.split('%sl' % lowerstring)[1].split('/')[0]
                string=lowerstring
    else:
        lag=dir.split('%sl' % upperstring)[1].split('/')[0]
        string=upperstring
    pops=loadtxt('%s/Populations.dat' % dir)
    k=len(pops)
    pylab.figure()
    bound_volumes=loadtxt('%s/Gens_k%s_l%s_pocketdistance_volume.dat' % (dir, k, lag))
    axis=loadtxt('%s/Gens_k%s_l%s_pocketdistance_axis.dat' % (dir,  k, lag))
    standard_frees=loadtxt('%s/Gens_k%s_l%s_pocketdistance_free.dat' % (dir, k, lag))
    pylab.plot(axis, [ref]*len(standard_frees), color='black', label='exp')
    pylab.plot(axis, standard_frees, color='red', label='standard free')
    pylab.ylim(-12,0)
    pylab.xlabel('P-L COM Distance')
    pylab.ylabel('standard free energy (kcal/mol)')
    pylab.legend()
    pylab.title('%s' % dir)
    pylab.savefig('%s/free.png' % dir)
    pylab.figure()
    pylab.plot(axis, bound_volumes, color='b', label='weighted bound volume')
    pylab.xlabel('P-L COM Distance')
    pylab.ylabel('volume (A^3)')
    pylab.title('%s' % dir)
    pylab.savefig('%s/volume.png' % dir)
    pylab.show()


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='input dir')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir)


