from msmbuilder import io
import os
import optparse
from msmbuilder import Project
import numpy

def main(file):
    ass=io.loadh(file)
    dir=os.path.dirname(file)
    base=os.path.basename(file)
    newdir='%s/subsample' % dir
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    p=Project.load_from('%s/ProjectInfo.yaml' % dir.split('Data')[0])
    data=dict()
    totals=dict()
    iterations=int(ass['arr_0'].shape[1]/10.0)
    start=max(p.traj_lengths)
    for iter in range(0, iterations):
        new=start-10
        if new < 10:
            break
        totals[new]=0
        data[new]=-numpy.ones((ass['arr_0'].shape[0], new), dtype=int)
        for i in range(0, ass['arr_0'].shape[0]):
            data[new][i]=ass['arr_0'][i][:new]
            frames=numpy.where(data[new][i]!=-1)[0]
            totals[new]+=len(frames)
        start=new

    ohandle=open('%s/times.h5' % (newdir), 'w')
    for key in sorted(data.keys()):
        print data[key].shape
        print "total time is %s" % totals[key]
        ohandle.write('%s\t%s\t%s\n' % (data[key].shape[0], data[key].shape[1], totals[key]))
        #io.saveh('%s/%s_sub%s.h5' % (newdir, base.split('.h5')[0], key), data[key])

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='file',
                      help='assignments file')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(file=options.file)

