#!/home/mlawrenz/Programs/epd-7.3-1-rh5-x86_64/bin/python

from msmbuilder import Trajectory, Project, io
from msmbuilder import Project
import optparse
import numpy
import os
from numpy import linalg


def main(dir, coarse , lag):
    proj=Project.load_from('ProjectInfo.yaml')
    ass=io.loadh('%s/msml%s_coarse_r10_d%s/Assignments.Fixed.h5' % (dir, lag, coarse))
    data=dict()
    rmsd=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.rmsd.dat' % (dir, coarse), usecols=(2,))
    data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.selfrmsd.dat' % (dir, coarse))
    com=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.vmd_com.dat' % (dir, coarse), usecols=(1,))
    com=[i/com[0] for i in com]
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/msml%s_coarse_r10_d%s/Populations.dat' % (dir, lag, coarse))
    map=numpy.loadtxt('%s/msml%s_coarse_r10_d%s/Mapping.dat' % (dir, lag, coarse))

    data['com']=numpy.array(data['com'])
    data['rmsd']=numpy.array(data['rmsd'])
    unbound1=numpy.where((data['rmsd']>7.0)&(data['com']>3.0))[0]
    print "%s strict unbound states" % len(data['rmsd'][unbound1])
    unbound2=numpy.where((data['rmsd']>6.0)&(data['com']>2.0))[0]
    print "%s loose unbound states" % len(data['rmsd'][unbound2])
    bound1=numpy.where((data['rmsd']<3.0)&(data['com']<1.1))[0]
    print "%s strict bound states" % len(data['rmsd'][bound1])
    bound2=numpy.where((data['rmsd']<4.0)&(data['com']<1.1))[0]
    print "%s loose bound states" % len(data['rmsd'][bound2])

    dirs=dict()
    dirs['strict']='%s/msml%s_coarse_r10_d%s/tpt-strict' % (dir, lag, coarse)
    dirs['loose']='%s/msml%s_coarse_r10_d%s/tpt-loose' % (dir, lag, coarse)

    unbinds=[unbound1, unbound2]
    binds=[bound1, bound2]
    names=['strict', 'loose']
    for (u, b, name) in zip(unbinds, binds, names):
        if not os.path.exists(dirs[name]):
            os.mkdir(dirs[name])
        ohandle=open('%s/unbound_%s_states.txt' % (dirs[name], name), 'w')
        ghandle=open('%s/gen_unbound_%s_states.txt' % (dirs[name], name), 'w')
        for i in u:
            ghandle.write('%s\n' % i)
            if map[i]!=-1:
                t=proj.get_random_confs_from_states(ass['arr_0'], [map[i],], 100)
                print "writing unbound state %s xtc" % map[i]
                t[0].save_to_xtc('%s/unbound_state%s.xtc' % (dirs[name], int(map[i])))
                ohandle.write('%s\n' % int(map[i]))
    
        ohandle=open('%s/bound_%s_states.txt' % (dirs[name], name), 'w')
        ghandle=open('%s/gen_bound_%s_states.txt' % (dirs[name], name), 'w')
        for i in b:
            ghandle.write('%s\n' % i)
            if map[i]!=-1:
                t=proj.get_random_confs_from_states(ass['arr_0'], [map[i],], 100)
                print "writing bound state %s xtc" % map[i]
                t[0].save_to_xtc('%s/bound_state%s.xtc' % (dirs[name], int(map[i])))
                ohandle.write('%s\n' % int(map[i]))

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-c', '--coarse', dest='coarse',
                      help='coarse grain cutoff')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, coarse=options.coarse, lag=options.lag)

