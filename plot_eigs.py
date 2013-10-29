import numpy, pylab
import os
import glob
from parse_eig import parse_eiginfo

def get_ref_sys(sys):
    if sys=='LG2':
        value=2.77
    elif sys=='LG3': 
        value=2.43
    elif sys=='LG6':
        value=3.40
    elif sys=='LG9':
        value=3.30
    elif sys=='AND':
        value=2.43 # same as LG3 cuz smallest
    return value

systems=['LG2', 'LG9', 'LG3',  'LG6', 'AND']
colors=['red', 'blue', 'pink', 'purple', 'grey']
lags=[10,20]
eig='0'
for lag in lags:
    for (color, sys) in zip(colors, systems):
        dir='/nobackup/mlawrenz/FKBP-FAH-results2/%s-step/Data_pairs_d0.5_s100_hybrid/subsample/' % sys
        ref_sub=numpy.loadtxt('%s/times.h5' % dir, usecols=(1,)) 
        ref_total=numpy.loadtxt('%s/times.h5' % dir, usecols=(2,))
        times=dict()
        for (i,j) in zip(ref_sub, ref_total):
            times[i]=(j*100.0)/(1000*1000)
        count=0
        totals=[]
        coms=[]
        for subdir in glob.glob('%s/sub*/' % dir):
            print subdir
            num=int(subdir.split('sample/sub')[1].split('/')[0])
            eigfile='%s/msml%s/eig-states/eiginfo.txt' % (subdir, lag)
            if not os.path.exists(eigfile):
                continue
            else:
                eigs=parse_eiginfo(eigfile)  
                if len(eigs.keys())==0:
                    continue
                else:
                    ref_com=get_ref_sys(sys)
                    free=-0.6*numpy.log(eigs[eig]['maxes']['mag'][0])
                    totals.append(times[num])
                    coms.append(eigs[eig]['maxes']['com'][0]-ref_com)
        pylab.scatter(totals, coms, label='%s' % sys, c=color)
        pylab.hold(True)
    pylab.title('Eig0 Highest Probability State from MSMl%s' % lag)
    pylab.xlabel('Total Simulation Time (microsec)')
    pylab.ylabel('P-L COM Distance from Crystal')
    pylab.ylim(0,10)
    pylab.legend(loc=9, frameon=False)
    pylab.savefig('eig0_p_l_msml%s.png' % lag, dpi=300)
    pylab.show()
