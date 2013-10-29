#!/home/mlawrenz/Programs/epd-7.3-1-rh5-x86_64/bin/python

from msmbuilder import Trajectory, Project, io
from msmbuilder import Project
import optparse
import numpy
import os
from numpy import linalg


def main(modeldir, gensfile, write=False):
    proj=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    data=dict()
    data['dist']=numpy.loadtxt('%s.prot_lig_distance.dat' % gensfile.split('.lh5')[0], usecols=(1,))
    data['rmsd']=numpy.loadtxt('%s.rmsd.dat' % gensfile.split('.lh5')[0])
    com=numpy.loadtxt('%s.vmd_com.dat' % gensfile.split('.lh5')[0], usecols=(1,))
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    data['com']=numpy.array(data['com'])
    data['rmsd']=numpy.array(data['rmsd'])
    unbound_touch1=numpy.where(data['dist'] > 20.0)[0]
    bound1=numpy.where(data['rmsd']<3)[0]
    print "%s super strict bound states" % len(data['dist'][bound1])
    print "%s super strict unbound states" % len(data['dist'][unbound_touch1])
    unbound_touch2=numpy.where(data['dist'] > 12.0)[0]
    bound2=numpy.where(data['rmsd']<7)[0]
    print "%s strict unbound states" % len(data['dist'][unbound_touch2])
    unbound_touch3=numpy.where(data['dist'] > 8.0)[0]
    bound3=numpy.where(data['com']<10)[0]
    print "%s medium unbound states" % len(data['rmsd'][unbound_touch3])
    unbound_touch4=numpy.where(data['dist'] > 4.0)[0]
    bound4=numpy.where(data['com']<15)[0]
    print "%s loose unbound states" % len(data['rmsd'][unbound_touch4])

    dirs=dict()

    #unbinds=[unbound_touch, unbound_touch, unbound_touch, unbound_touch]
    #binds=[bound1, bound2, bound3, bound4]
    unbinds=[unbound_touch1, unbound_touch2, unbound_touch3, unbound_touch4]
    binds=[bound1, bound1, bound1, bound1]
    names=['super-strict', 'strict', 'medium', 'loose']
    for (u, b, name) in zip(unbinds, binds, names):
        dirs[name]='%s/tpt-rmsd-%s' % (modeldir, name)
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
    parser.add_option('-g', '--gensfile', dest='gensfile',
                          help='gens files')
    parser.add_option('-w', action="store_true", dest="write")
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.write==True:
        main(modeldir=options.dir, gensfile=options.gensfile,  write=True)
    else:
        main(modeldir=options.dir,gensfile=options.gensfile)

