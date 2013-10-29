#!/bin/python
import numpy
import pylab
import optparse

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)


def make(input):
    colors=['#0000ff','#00ccff','#00ffcc','#ff0000','#ff00cc','#ff66ff','#ff6600','#ff9900','#6600cc','#330099','#00ff00']
    # 1 lag=20 steps= 2 ns
    lags=numpy.loadtxt(input, usecols=(0,))
    times=numpy.loadtxt(input, usecols=(1,))
    data=dict()
    for (x,y) in zip(lags, times):
        if x in data.keys():
            data[x].append(y)
        else:
            data[x]=[]
            data[x].append(y)
    
    return data, colors

def main(input, sys):
    sys=input.split('-step')[0].split('/')[1]
    print sys
    data, colors=make(input)
    length=len(data[data.keys()[0]])
    for (count, color) in zip(range(0, length), colors):
        ydata=[]
        for key in sorted(data.keys()):
            ydata.append(data[key][count])
        pylab.scatter([i*100/1000.0 for i in sorted(data.keys())], [numpy.log10(i*100/1000.0) for i in ydata], c=color, edgecolor=color, marker='o')
        pylab.hold(True)
        #pylab.plot([i*100/1000.0 for i in sorted(data.keys())], [numpy.log10(i*100/1000.0) for i in ydata], c=color)
    pylab.plot([10]*8, numpy.arange(0.5,4.5,0.5), 'k--')
    xlabel=pylab.xlabel('Lag Time $\\tau$ (ns)',size='large')
    pylab.ylim(0.5, 4.0, 0.5)
    #pylab.xlim(5, 30)
    #loc,labels=pylab.xticks(range(50, 5000, 200) , [round(i*10.0/1000,2) for i in range(50, 5000, 200)])
    ylabel=pylab.ylabel('log(ns)',size='large')
    #pylab.xlim(-5,1020)
    #lg=pylab.legend(loc=9)
    #lg.get_frame().set_linewidth(0)
    pylab.savefig('%s.%s.png' % (input.split('.dat')[0], sys), dpi=300)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input',
                      help='input file')
    parser.add_option('-s', '--system', dest='sys',
                      help='system name')
    parser.add_option('-k', '--kclusters', dest='k',
                      help='number of clusters')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(input=options.input, sys=options.sys)


