#!/bin/python
from scipy import stats
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt, io
import glob
import matplotlib.gridspec as gridspec
import random
from scipy.io import *
from msmbuilder import Conformation
from msmbuilder import MSMLib
import optparse
import numpy
import os
from numpy import linalg
import pylab

def map_size(x):
    if x==0.5:
        size=500
    elif abs(x-0.5) < 0.1 and abs(x-0.5) > 0:
        size=300
    elif abs(x-0.5) < 0.2 and abs(x-0.5) > 0.1:
        size=200
    elif abs(x-0.5) < 0.3 and abs(x-0.5) > 0.2:
        size=100
    else:
        size=50
    return size

def get_ss_matrix(residues, frames, struct, stride=0):
    unique_res=sorted(set(residues))
    unique_frames=sorted(set(frames))
    matrix=numpy.zeros((len(unique_res), len(unique_frames)-1))
    bound_matrix=numpy.zeros((len(unique_res), 1))
    res_num=dict()
    helix=True
    for (n, r) in enumerate(unique_res):
        res_num[r]=n
    check_matrix=numpy.zeros((len(unique_res), 20), dtype=numpy.str)
    refframe=0
    for n in range(0, len(residues)):
        res=int(residues[n])
        res_map=res_num[res]
        frame=int(frames[n])
        str=struct[n]
        if frame==0:
            bound_matrix[res_map][frame]=number_code(str)
            continue
        else:
            frame=frame-1
            matrix[res_map][frame]=number_code(str)
    return bound_matrix, matrix

def number_code(input):
    data=dict()
    data['E']=0
    data['B']=1
    data['C']=2
    data['T']=3
    data['I']=4
    data['G']=5
    data['H']=6
    return data[input]

def structure_code(input):
    data=dict()
    data[0]='Beta Extend' 
    data[1]='Bridge'        
    data[2]='Coil'
    data[3]='Helix Turn'
    data[4]='Pi Helix'
    data[5]='3-10 Helix'
    data[6]='Alpha Helix'   
    return data[input]

def main(file):
    data=dict()
    data['selfrmsd']=numpy.loadtxt('%s.selfrmsd.dat' % file.split('.xtc')[0])
    #data['selfhelix']=numpy.loadtxt('%s.selfhelixrmsd.dat' % file.split('.xtc')[0])
    data['helix']=numpy.loadtxt('%s.helixrmsd.dat' % file.split('.xtc')[0])
    #data['rmsd']=numpy.loadtxt('%s.rmsd.dat' % file.split('.xtc')[0])
    #com=numpy.loadtxt('%s.vmd_com.dat' % file.split('.xtc')[0], usecols=(1,))
    #com=[i/com[0] for i in com]
    #data['com']=com[1:]

    residues=numpy.loadtxt('%s.ss.dat' % file, usecols=(0,), dtype=int)
    frames=numpy.loadtxt('%s.ss.dat' % file, usecols=(1,), dtype=int)
    struct=numpy.loadtxt('%s.ss.dat' % file, usecols=(2,), dtype=str)
    ref_matrix, ss_matrix=get_ss_matrix(residues, frames, struct, stride=20)
    cm=pylab.cm.get_cmap('RdYlBu_r', 7)
    fig=pylab.figure(figsize=(8,8))
    gs=gridspec.GridSpec(1, 2, width_ratios=[8,1])
    ax1=fig.add_subplot(gs[0] )
    ax2=fig.add_subplot(gs[1])
    fig1=ax1.pcolor(ss_matrix, cmap=cm, vmax=7, vmin=0)
    ax1.yaxis.set_ticks(range(0, 22))
    ax1.yaxis.set_ticklabels(range(94,116))
    ax1.set_ylim([0,22])
    if 'path' in file:
        ax1.xaxis.set_ticks(numpy.arange(0, max(set(frames))+1, 20))
        ax1.xaxis.set_ticklabels([int(l/20.0) for l in numpy.arange(0, max(set(frames))+1, 20)])
    else:
        ax1.xaxis.set_ticks(numpy.arange(0, max(set(frames))+1, 100))
        ax1.xaxis.set_ticklabels([round(((l*10.0)/1000.0), 2) for l in numpy.arange(0, max(set(frames))+1, 100)])
    #cbar = ax1.colorbar(fig, ticks=range(0,8))
    #cbar.ax.set_yticklabels([structure_code(i) for i in range(0,8)]) # vertical
    fig2=ax2.pcolor(ref_matrix, cmap=cm, vmax=7, vmin=0)
    ax2.set_ylim([0,22])
    ax2.yaxis.set_ticks(range(0, 22))
    ax2.yaxis.set_ticklabels([' ']*len(range(0, 22)))
    ax2.xaxis.set_ticks(numpy.arange(0,1.1, 1))
    ax2.xaxis.set_ticklabels([' ', ' '])
    cbar = pylab.colorbar(fig2, ticks=range(0,7))
    cbar.ax.set_yticklabels([structure_code(i) for i in range(0,7)]) # vertical
    pylab.subplots_adjust(left = 0.10, right = 0.8, bottom = 0.10, top = 0.85, wspace = 0.1, hspace = 0.1)
    l,b,w,h=pylab.gca().get_position().bounds
    ll, bb, ww, hh=cbar.ax.get_position().bounds
    #cbar.ax.set_position([ll, b+0.1*h, ww, h*0.8])
    cbar.ax.set_position([ll-0.6*w, b+0.1*h, ww, h])
    ax1.set_ylabel('P53 Resid')
    if 'path' in file:
        ax1.set_xlabel('Pathway State')
    else:
        ax1.set_xlabel('MicroSeconds')
    ax2.set_xlabel('Bound P53')
    pylab.savefig('%s_ss.png' % (file.split('.xtc')[0]), dpi=300)

    count=0
    colors=[]
    for n in range(0, len(data['helix'])):
        colors.append(count/(len(data['helix'])/20.0))
        if n % 20 ==0:
            count+=1
    for op in data.keys():
        if op=='helix':
            continue
        else:
            pylab.figure()
            pylab.scatter(data['helix'], data[op], c=colors, alpha=0.4)
            pylab.subplots_adjust(left = 0.1, right = 1.02, bottom = 0.10, top = 0.85, wspace = 0, hspace = 0)
            CB=pylab.colorbar()
            #l,b,w,h=pylab.gca().get_position().bounds
            #ll, bb, ww, hh=CB.ax.get_position().bounds
            #CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
            CB.set_ticks(numpy.arange(0,1.1,1))
            CB.set_ticklabels(['Init', 'Final'])
            CB.set_label('Trajectory Progress')
            pylab.xlabel('P53-Binding Site COM')
            pylab.ylabel('P53 RMSD')
            pylab.xlim(0,50)
            pylab.ylim(0,12)
            if 'micro' in file:
                pylab.title('%s $\mu$s' % file.split('_')[-1].split('.xtc')[0])
            elif 'path' in file:
                dir=os.path.dirname(file)
                pylab.title('%s' % file.split('_sample')[0])
                paths=io.loadh('%s/Paths.h5' % dir)
                map=numpy.loadtxt('%s/Mapping.dat' % dir.split('tpt')[0])
                p=file.split('path')[1].split('_sample')[0]
                frames=numpy.where(paths['Paths'][int(p)]!=-1)[0]
                gen_path=[]
                for state in paths['Paths'][int(p)][frames]:
                    gen_path.append(numpy.where(map==state)[0])
                print gen_path
                print paths['fluxes'][int(p)]/paths['fluxes'][0]
            else:
                pylab.title('%s' % file.split('.xtc')[0])
        pylab.savefig('%s_%s_helix.png' % (file.split('.xtc')[0], op), dpi=300)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='file',
                      help='trajectory file')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(file=options.file)

