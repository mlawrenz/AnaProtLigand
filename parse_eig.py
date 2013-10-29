import numpy, pylab

def parse_eiginfo(file):
    ohandle=open(file)
    eigs=dict()
    for line in ohandle.readlines():
        if 'maxes' in line.split():
            num=line.split()[0].split('eig')[1]
            eigs[num]=dict()
            eigs[num]['maxes']=dict()
            eigs[num]['maxes']['rmsd']=[]
            eigs[num]['maxes']['mag']=[]
            eigs[num]['maxes']['com']=[]
            type='maxes'
        elif "mins" in line.split():
            num=line.split()[0].split('eig')[1]
            if num not in eigs.keys():
                eigs[num]=dict()
            eigs[num]['mins']=dict()
            eigs[num]['mins']['rmsd']=[]
            eigs[num]['mins']['mag']=[]
            eigs[num]['mins']['com']=[]
            type='mins'
        try:
            test=int(line.split()[0])
            eigs[num][type]['mag'].append(float(line.split()[2]))
            eigs[num][type]['rmsd'].append(float(line.split()[3]))
            eigs[num][type]['com'].append(float(line.split()[4]))
        except ValueError:
            pass
    return eigs
