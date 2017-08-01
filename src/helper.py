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


[(1, 1, 121),
 (1, 1, 122),
 (1, 1, 123),
 (1, 1, 124),
 (1, 1, 125),
 (1, 1, 126),
 (1, 1, 127),
 (1, 2, 7),
 (1, 3, 5),
 (1, 120, 2),
 (1, 121, 2),
 (1, 122, 2),
 (1, 123, 2),
 (1, 124, 2),
 (1, 125, 2),
 (1, 126, 2),
 (2, 1, 61),
 (2, 1, 62),
 (2, 1, 63),
 (2, 2, 6),
 (2, 60, 2),
 (2, 61, 2),
 (2, 62, 2),
 (3, 1, 41),
 (3, 1, 42),
 (3, 40, 2),
 (3, 41, 2),
 (4, 1, 31),
 (4, 2, 5),
 (4, 5, 3),
 (4, 30, 2),
 (5, 1, 25),
 (5, 24, 2),
 (6, 1, 21),
 (6, 4, 3),
 (6, 20, 2),
 (7, 1, 18),
 (7, 17, 2),
 (9, 1, 14),
 (9, 13, 2),
 (11, 1, 11),
 (11, 10, 2),
 (14, 1, 9),
 (14, 8, 2),
 (18, 1, 7),
 (18, 2, 3),
 (18, 6, 2),
 (21, 1, 6),
 (21, 5, 2),
 (25, 1, 5),
 (25, 4, 2),
 (31, 1, 4),
 (31, 3, 2),
 (41, 1, 3),
 (41, 2, 2),
 (42, 1, 3),
 (42, 2, 2),
 (61, 1, 2),
 (62, 1, 2),
 (63, 1, 2)]
 
 
 near misses removed:
 [(1, 1, 124),
 (1, 2, 7),
 (1, 3, 5),
 (1, 123, 2),
 
 (2, 1, 62),
 (2, 2, 6),
 (2, 61, 2),
 (3, 1, 41),
 (3, 1, 42),
 (3, 40, 2),
 (3, 41, 2),
 (4, 1, 31),
 (4, 2, 5),
 (4, 5, 3),
 (4, 30, 2),
 (5, 1, 25),
 (5, 24, 2),
 (6, 1, 21),
 (6, 4, 3),
 (6, 20, 2),
 (7, 1, 18),
 (7, 17, 2),
 (9, 1, 14),
 (9, 13, 2),
 (11, 1, 11),
 (11, 10, 2),
 (14, 1, 9),
 (14, 8, 2),
 (18, 1, 7),
 (18, 2, 3),
 (18, 6, 2),
 (21, 1, 6),
 (21, 5, 2),
 (25, 1, 5),
 (25, 4, 2),
 (31, 1, 4),
 (31, 3, 2),
 (41, 1, 3),
 (41, 2, 2),
 (42, 1, 3),
 (42, 2, 2),
 (62, 1, 2)]

 

"""
"""

tmp =  find_trees2(124,tolerance=1).keys()

full: 
{(1, 1, 123): (1, 1, 123),
 (1, 1, 124): (1, 1, 124),
 (1, 1, 125): (1, 1, 125),
 (1, 122, 2): (1, 122, 2),
 (1, 123, 2): (1, 123, 2),
 (1, 124, 2): (1, 124, 2),
 (2, 1, 62): (2, 1, 62),
 (2, 61, 2): (2, 61, 2),
 (3, 1, 41): (3, 1, 41),
 (3, 40, 2): (3, 40, 2),
 (4, 1, 31): (4, 1, 31),
 (4, 2, 5): (4, 2, 5),
 (4, 5, 3): (4, 5, 3),
 (4, 30, 2): (4, 30, 2),
 (5, 1, 25): (5, 1, 25),
 (5, 24, 2): (5, 24, 2),
 (25, 1, 5): (25, 1, 5),
 (25, 4, 2): (25, 4, 2),
 (31, 1, 4): (31, 1, 4),
 (31, 3, 2): (31, 3, 2),
 (41, 1, 3): (41, 1, 3),
 (41, 2, 2): (41, 2, 2),
 (62, 1, 2): (62, 1, 2)}
 
 
 
 near misses removed:
 {(1, 1, 124): (1, 1, 124),
 (1, 123, 2): (1, 123, 2),
 (2, 1, 62): (2, 1, 62),
 (2, 61, 2): (2, 61, 2),
 (3, 1, 41): (3, 1, 41),
 (3, 40, 2): (3, 40, 2),
 (4, 1, 31): (4, 1, 31),
 (4, 2, 5): (4, 2, 5),
 (4, 5, 3): (4, 5, 3),
 (4, 30, 2): (4, 30, 2),
 (5, 1, 25): (5, 1, 25),
 (5, 24, 2): (5, 24, 2),
 (25, 1, 5): (25, 1, 5),
 (25, 4, 2): (25, 4, 2),
 (31, 1, 4): (31, 1, 4),
 (31, 3, 2): (31, 3, 2),
 (41, 1, 3): (41, 1, 3),
 (41, 2, 2): (41, 2, 2),
 (62, 1, 2): (62, 1, 2),
 (124,1, 1): (124,1, 1)}
 
 
"""