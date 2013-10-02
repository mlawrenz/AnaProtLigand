import numpy, pylab
from msmbuilder import io
bw=dict()
bw_dcds=numpy.loadtxt('dimer_mdcrds/bw-trajs.txt', usecols=(0,), dtype=str)
bw_lengths=numpy.loadtxt('dimer_mdcrds/bw-trajs.txt', usecols=(1,))
for (i,j) in zip(bw_dcds, bw_lengths):
    name=i.split('bw-')[1].split('_nowat_')[0]
    bw[name]=j

totals=dict()
orig_dcds=numpy.loadtxt('dimer_mdcrds/orig-trajs.txt', usecols=(0,), dtype=str)
orig_lengths=numpy.loadtxt('dimer_mdcrds/orig-trajs.txt', usecols=(1,))
for (i,j) in zip(orig_dcds, orig_lengths):
    name=i.split('_nowat_')[0]
    totals[name]=j

dcds=numpy.loadtxt('mapped_trajs.txt', usecols=(0,), dtype=str)
xtcs=numpy.loadtxt('mapped_trajs.txt', usecols=(1,), dtype=str)
mapping=dict()
for (i,j) in zip(dcds, xtcs):
    name=i.split('_nowat_')[0]
    mapping[j]=name
    
#for j in mapping.keys():
#    total=mapping[j]+bw[j]
#    if total!=totals[j]:
#        print "problem"
#        import pdb
#        pdb.set_trace()

ohandle=open('d6/msml1000_coarse_r10_d20/traj_frames.txt', 'w')
ass=io.loadh('d6/msml1000_coarse_r10_d20/Assignments.Fixed.h5')
mapfile=numpy.loadtxt('d6/msml1000_coarse_r10_d20/Mapping.dat')
sample=False
for state in sorted(set(ass['arr_0'].flatten())):
    if state!=-1:
        traj=numpy.where(ass['arr_0']==state)[0]
        frames=numpy.where(ass['arr_0']==state)[1]
        indices=numpy.random.random_integers(0, len(traj)-1, len(traj))
        for ind in indices:
            traj_ind=traj[ind]
            mapped_traj=mapping['trj%s' % traj_ind]
            if mapped_traj in bw.keys():
                minval=bw[mapped_traj]
            else:
                minval=0
            location=numpy.where((traj==traj_ind)&(frames>minval))[0]
            if not location.size:
                sample=False
                continue
            else:
                sample=True
                break
        if sample==False:
            print "no sample for state %s" % state
            continue
        else:
            frame_ind=numpy.random.random_integers(0, len(location)-1, 1)
            frame_number=frames[location][frame_ind]
            if state!=ass['arr_0'][traj_ind,frame_number]:
                print "problem"
                import pdb
                pdb.set_trace()
            else:
                name=mapping['trj%s' % traj_ind]
                if name not in bw.keys():
                    print mapping['trj%s' % traj_ind], "no bw", int(frame_number[0]-minval)
                else:
                    print mapping['trj%s' % traj_ind], int(bw[name]), int(frame_number[0]-minval)
                ohandle.write('%s\t%s\n' % (mapping['trj%s' % traj_ind], int(frame_number[0]-minval)))
