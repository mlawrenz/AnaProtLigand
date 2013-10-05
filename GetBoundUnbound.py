#!/home/mlawrenz/Programs/epd-7.3-1-rh5-x86_64/bin/python

from msmbuilder import Trajectory, Project, io
from msmbuilder import Project
import optparse
import numpy
import os
from numpy import linalg


def main(modeldir, write=False):
    proj=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    data=dict()
    data['dist']=numpy.loadtxt('%s/prot_lig_distance.dat' % modeldir, usecols=(1,))
    data['rmsd']=numpy.loadtxt('%s/Gens.rmsd.dat' % modeldir, usecols=(2,))
    com=numpy.loadtxt('%s/Gens.vmd_com.dat' % modeldir, usecols=(1,))
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    data['com']=numpy.array(data['com'])
    data['rmsd']=numpy.array(data['rmsd'])
    unbound_touch=numpy.where(data['dist'] > 4.0)[0]
    print "%s non touching unbound states" % len(data['dist'][unbound_touch])
    bound1=numpy.where(data['com']<5)[0]
    print "%s strict bound states" % len(data['dist'][bound1])
    bound2=numpy.where(data['com']<10)[0]
    print "%s medium bound states" % len(data['rmsd'][bound2])
    bound3=numpy.where(data['com']<15)[0]
    print "%s loose bound states" % len(data['rmsd'][bound3])

    dirs=dict()

    unbinds=[unbound_touch, unbound_touch, unbound_touch]
    binds=[bound1, bound2, bound3,]
    names=['strict', 'medium', 'loose']
    for (u, b, name) in zip(unbinds, binds, names):
        dirs[name]='%s/tpt-%s' % (modeldir, name)
        if not os.path.exists(dirs[name]):
            os.mkdir(dirs[name])
        ohandle=open('%s/unbound_%s_states.txt' % (dirs[name], name), 'w')
        ghandle=open('%s/gen_unbound_%s_states.txt' % (dirs[name], name), 'w')
        for i in u:
            ghandle.write('%s\n' % i)
            if map[i]!=-1:
                ohandle.write('%s\n' % int(map[i]))
                if write==True:
                    if not os.path.exists('%s/unbound_state%s.xtc' % (dirs[name], int(map[i]))):
                        t=proj.get_random_confs_from_states(ass['arr_0'], [map[i],], 100)
                        print "writing unbound state %s xtc" % map[i]
                        t[0].save_to_xtc('%s/unbound_state%s.xtc' % (dirs[name], int(map[i])))
    
        ohandle=open('%s/bound_%s_states.txt' % (dirs[name], name), 'w')
        ghandle=open('%s/gen_bound_%s_states.txt' % (dirs[name], name), 'w')
        for i in b:
            ghandle.write('%s\n' % i)
            if map[i]!=-1:
                ohandle.write('%s\n' % int(map[i]))
                if write==True:
                    if not os.path.exists('%s/bound_state%s.xtc' % (dirs[name], int(map[i]))):
                        t=proj.get_random_confs_from_states(ass['arr_0'], [map[i],], 100)
                        print "writing bound state %s xtc" % map[i]
                        t[0].save_to_xtc('%s/bound_state%s.xtc' % (dirs[name], int(map[i])))

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-w', action="store_true", dest="write")
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.write==True:
        main(modeldir=options.dir,  write=True)
    else:
        main(modeldir=options.dir,  write=False)

