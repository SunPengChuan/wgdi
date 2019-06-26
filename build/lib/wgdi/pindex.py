import os
import random
import re
import sys

import numpy as np
import pandas as pd


class dotplot():
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def run():
        alignment = pd.read_csv(self.alignment,header=None,sep='\t',index_col=0)
        alignment.replace('\w+',1,regex=True,inplace=True)
        alignment.replace('.',0,inplace=True)
        alignment.fillna(0,inplace=True)
        gff = pd.read_csv(self.gff, sep="\t", header=None)
        
        result = []
        for i in range(10):
            newdata = genome(data)
            x = times(newdata)
            result.append(x)
        print(sum(result)/len(result))

def genome(data):
    chr = [k[0] for k in data]
    d = {k: chr.count(k) for k in set(chr)}
    newdata = []
    for i in range(1, len(data[0])):
        c = []
        for j in sorted(d.keys()):
            b = [k[i] for k in data if k[0] == j]
            c.append(subgenome(b))
        newdata.append(c)
    return newdata


def times(newdata):
    a = []



    
    for i in range(len(newdata)-1):
        for j in range(i+1, len(newdata)):
            b = []
            for k in range(len(newdata[i])):
                chr1, chr2 = newdata[i][k], newdata[j][k]
                b.append(Pindex(chr1, chr2))
            print(b)
            c = sum([abs(k[0]) for k in b])/sum([k[1] for k in b])
            a.append(c)
    print(a)
    return sum(a)/len(a)

def Pindex(chr1, chr2):
    r1 = retain(chr1)
    r2 = retain(chr2)
    r = []
    for i in range(len(r2)):
        if(r1[i] == 0 and r2[i] == 0):
            r.append(0)
            continue
        if((r1[i]-r2[i])/(r1[i]+r2[i])*0.5 > 0.05):
            r.append(1)
        elif((r2[i]-r1[i])/(r1[i]+r2[i])*0.5 > 0.05):
            r.append(-1)
            #print(len(r))
        else:
            r.append(0)
    #r=[i for i in r if i in [-1,1]]
    #print(chr1,len(r))
    r1 = [sum(r), len(r)]
    #if len(r)>0:
    #    r1=sum(r)/len(r)
    return r1

    def retain(arr):
        a = []
        for i in range(0, len(arr), 100):
            start, end = i-50, i+50
            genenum, retainnum = 0, 0
            for j in range(start, end):
                if((j >= int(len(arr))) or (j < 0)):
                    continue
                else:
                    retainnum += arr[j]
                    genenum += 1
            a.append(float(retainnum/genenum))
        return a



