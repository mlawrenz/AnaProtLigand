import numpy, glob
import os
import optparse
from msmbuilder.geometry import dihedral as _dihedralcalc
from msmbuilder import Trajectory, Project
from msmbuilder import metrics
import pylab



project = Project.load_from('ProjectInfo.yaml')
name='W59'
bound_cutoff=0.6
unbound_cutoff=1.5
#unbound_cutoff=2.0
bind_event=0
unbind_event=0
avg=0
count=0
total=[]
for n in range(0, project.n_trajs):
    unbound=False
    bound=False
    data=numpy.loadtxt('./Trajectories-metric/trj%s_%spair.dat' % (n, name))
    #t=project.load_traj(9982)
    for i in data:
        avg=((count*avg)+float(i))/(count+1)
        count+=1
        if i <= bound_cutoff:
            if avg >=unbound_cutoff:
                bind_event+=1
                print "bind %s to %s in traj %s" % (avg, i, n)
                avg=0
                count=0
            else:
                pass
        elif i >= unbound_cutoff:
            if avg <=bound_cutoff:
                unbind_event+=1
                print "unbind %s to %s in traj %s" % (avg, i, n)
                avg=0
                count=0
#pylab.figure()
#pylab.hist(total)
#pylab.savefig('W59dist.png')
print "unbind events: %s " % unbind_event
print "bind events: %s " % bind_event
ohandle=open('bind_events_avg.txt', 'w')
ohandle.write('bind\t%s\n' % bind_event)
ohandle.write('unbind\t%s\n' % unbind_event)
