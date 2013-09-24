from msmbuilder import Trajectory
import optparse
import numpy
import glob
import os
from numpy import linalg

def main(input, atoms):
    contacts=numpy.loadtxt(atoms, dtype=int, ndmin=2)
    print contacts.shape
    for n in range(0, contacts.shape[0]):
        atom1=int(contacts[n][0])+1
        atom2=int(contacts[n][1])+1
        t=Trajectory.load_trajectory_file(input)
        index1=numpy.where(t['AtomID']==atom1)[0]
        index2=numpy.where(t['AtomID']==atom2)[0]
        print t['ResidueNames'][index1], t['AtomNames'][index1]
        name1='r%s%s' % (t['ResidueNames'][index1][0], t['AtomNames'][index1][0])
        print t['ResidueNames'][index2], t['AtomNames'][index2]
        name2='r%s%s' % (t['ResidueNames'][index2][0], t['AtomNames'][index2][0])
        dist=[]
        for frame in range(0, t['XYZList'].shape[0]):
            diff=numpy.subtract(t['XYZList'][frame][index1], t['XYZList'][frame][index2])
            dist.append(linalg.norm(diff)*10)
        new=input.split('.lh5')[0]
        numpy.savetxt('%s.%s.%s.dat' % (new, name1, name2 ), dist)


def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input',
                      help='input file')
    parser.add_option('-p', '--pairs', dest='pairs',
                      help='atom pairs')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(input=options.input, atoms=options.pairs)



