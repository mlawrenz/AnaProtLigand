#!/home/mlawrenz/Programs/epd-7.3-1-rh5-x86_64/bin/python
import pylab

from msmbuilder import Trajectory, Project, io
from msmbuilder import Project
import optparse
import numpy
import os
from numpy import linalg


def main(modeldir):
    proj=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    data=dict()
    data['dist']=numpy.loadtxt('%s/prot_lig_distance.dat' % modeldir, usecols=(1,))
    data['rmsd']=numpy.loadtxt('%s/Gens.rmsd.dat' % modeldir, usecols=(2,))
    com=numpy.loadtxt('%s/Gens.vmd_com.dat' % modeldir, usecols=(1,))
    refcom=com[0]
    data['com']=com[1:]
    data['com']=numpy.array(data['com'])
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)
    frames=numpy.where(map!=-1)[0]
    pylab.scatter(data['com'][frames], data['rmsd'][frames])
    pylab.scatter([refcom,], [0,], c='k', marker='x', s=100)
    pylab.xlabel('P-L COM')
    pylab.ylabel('P-L RMSD')
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.dir)

