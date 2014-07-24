import numpy as np
import math


def _calc_powersum(nsib,nlevels):
    if nlevels < 1:
        print 'n levels cannot be < 1'
        return 0
    elif nlevels == 1:
        return nsib
    else:
        return nsib*(1 + _calc_powersum(nsib,nlevels-1))


def calc_num_segments(nb,nsib,nlevels):
    nseg_branch = _calc_powersum(nsib,nlevels-1)+1
    ntotal= nseg_branch*nb
    print 'Number of segments: per branch = ', nseg_branch
    print '                    total      = ', ntotal
    return nseg_branch, ntotal

def find_trees(num_segments,tolerance=3):
    nb = 1
    ntot = 0
    collection = []
    while ntot <= num_segments+tolerance: #loop for nb
        print 'loop A - nb = ',nb
        ns = 1
        while ntot <= num_segments+tolerance:#loop for ns
            print 'loop B - ns = ',ns
            nl = 2
            while ntot <= num_segments+tolerance:#loop for nl
                print 'loop C - nl = ',nl
                nseg_b,ntot = calc_num_segments(nb,ns,nl)
            
                if ntot >  num_segments+tolerance:
                    break
                elif ntot >=  num_segments-tolerance:
                    collection.append((nb,ns,nl))
                nl += 1
            ns += 1
        nb +=1
        
    return collection

def find_trees2(num_segments,tolerance=3,nb=1,ns=1,nl=2,collection={}):
    if collection.has_key((nb,ns,nl)):
        return
    
    nseg_b,ntot = calc_num_segments(nb,ns,nl)
    
    if ntot >=  num_segments-tolerance and ntot <=  num_segments+tolerance :
        collection[(nb,ns,nl)] = (nb,ns,nl)
    
    #tmp_col = []
    if ntot <=  num_segments+tolerance:
        find_trees2(num_segments,tolerance,nb+1,ns,nl,collection)
        find_trees2(num_segments,tolerance,nb,ns+1,nl,collection)
        find_trees2(num_segments,tolerance,nb,ns,nl+1,collection)
    #if ntot >=  num_segments-tolerance:
    #    collection.append((nb,ns,nl))
    if ntot > num_segments+tolerance:
        return collection
    return collection




"""
tmp =  find_trees2(124,tolerance=3).keys()
tmp.sort()
print tmp
"""