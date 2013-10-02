#!/bin/python
import pickle
import optparse


def coors_to_dict(file):
    data=dict()
    fhandle=open(file)
    for line in fhandle.readlines():
        frame=int(line.split()[0])
        state=frame-1
        if state not in data.keys():
            print "on state %s" % state
            data[state]=[]
        data[state].append([float (i) for i in line.split()[1:]])
    ofile=open('%s.pickle' % file.split('.dat')[0], 'w')
    pickle.dump(data, ofile)
    ofile.close()
    return data

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='file',
                      help='file')
    (options, args) = parser.parse_args()
    return (options, args)



if __name__ == "__main__":
    (options, args) = parse_commandline()
    coors_to_dict(file=options.file)

