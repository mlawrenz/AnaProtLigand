#!/bin/python
import numpy
import pylab
import optparse

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

def main(input, sys, k, pnorm=False):
    data, colors=make(input)
    length=len(data[data.keys()[0]])
    for (count, color) in zip(range(0, length), colors):
        ydata=[]
        for key in sorted(data.keys()):
            ydata.append(data[key][count])
        pylab.scatter(sorted(data.keys()), [numpy.log10(i*10/1000.0) for i in ydata], c=color, edgecolor=color, marker='o')
        pylab.hold(True)

    xlabel=pylab.xlabel('Lag Time $\\tau$ (ns)',size='large')
    #pylab.ylim(0.5, 3.5, 0.5)
    #pylab.xlim(50, 4500)
    #loc,labels=pylab.xticks(range(50, 5000, 200) , [round(i*10.0/1000,2) for i in range(50, 5000, 200)])
    #loc,labels=pylab.xticks(range(50, 5000, 200) , range(50, 5000, 200))
    #loc,labels=pylab.xticks(sorted(data.keys()),[int((x*2)/1000) for x in sorted(data.keys())])
    ylabel=pylab.ylabel('log(ns)',size='large')
    #pylab.xlim(-5,1020)
    #lg=pylab.legend(loc=9)
    #lg.get_frame().set_linewidth(0)
    pylab.title('%s %s State MSM' % (sys, k) )
    if pnorm==True:
        pylab.savefig('%s_Data_HybridPnorm_%s_Implied.png' % ( sys, k), dpi=300)
    else:
        pylab.savefig('%s_Data_Hybrid_%s_Implied.png' % ( sys, k), dpi=300)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input',
                      help='input file')
    parser.add_option('-s', '--system', dest='sys',
                      help='system name')
    parser.add_option('-k', '--kclusters', dest='k',
                      help='number of clusters')
    parser.add_option('-p', action="store_true", dest="pnorm")
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    if options.pnorm==True:
        main(input=options.input, sys=options.sys,k=options.k, pnorm=True)
    else:
        main(input=options.input, sys=options.sys,k=options.k )


