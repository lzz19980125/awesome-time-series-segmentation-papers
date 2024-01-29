___author__='Sorour Ekhtiari Amiri & Liangzhe Chen'
__email__='esorour@vt.edu'

import numpy as np
import math
from collections import defaultdict
import sys
import os
import ALP as falp
import time
import subprocess


def aib(data_dir, MinW, N_mdl, CalcMode, _dir, matlab_path):
    matlab_function = "AIB('" + _dir + "', '" + data_dir + "', '" + str(MinW) + "', '" + str(N_mdl) + "', '" + \
                      str(CalcMode) + "')"
    subprocess.call([matlab_path, "-nosplash", "-nodisplay", "-r", matlab_function])


def toposort(graph):
    print('toposort...')
    """http://code.activestate.com/recipes/578272-topological-sort/

    Dependencies are expressed as a dictionary whose keys are items
and whose values are a set of dependent items. Output is a list of
sets in topological order. The first set consists of items with no
dependences, each subsequent set consists of items that depend upon
items in the preceeding sets.

>>> print '\\n'.join(repr(sorted(x)) for x in toposort2({
...     2: set([11]),
...     9: set([11,8]),
...     10: set([11,3]),
...     11: set([7,5]),
...     8: set([7,3]),
...     }) )
[3, 5, 7]
[8, 11]
[2, 9, 10]

"""
    from functools import reduce
    data = defaultdict(set)
    for x, y in graph.items():
        for z in y:
            data[z[0]].add(x)

    # Ignore self dependencies.
    for k, v in data.items():
        v.discard(k)
    # Find all items that don't depend on anything.
    extra_items_in_deps = reduce(set.union, data.values()) - set(data.keys())
    # Add empty dependences where needed
    data.update({item:set() for item in extra_items_in_deps})
    while True:
        ordered = set(item for item, dep in data.items() if not dep)
        if not ordered:
            break
        yield ordered
        data = {item: (dep - ordered)
                for item, dep in data.items()
                    if item not in ordered}
    assert not data, "Cyclic dependencies exist among these items:\n%s" % '\n'.join(repr(x) for x in data.items())


def longestpathDAG(graph, startnode, endnode):
    print('longsetpath...')
    """http://www.geeksforgeeks.org/find-longest-path-directed-acyclic-graph/"""
    ### TOPOLOGICALLY SORT THE VERTICES
    order = []
    for part in toposort(graph):
        order.extend(list(part))
    # order.reverse()

    ### INITIALIZE DISTANCE MATRIX
    LOWDIST=-99999999999999999
    dist = dict((x, LOWDIST) for x in graph.keys())
    dist[startnode] = 0

    ### MAIN PART
    comesfrom = dict()
    for node in order: # u
        for nbr, t1, t2, t3, nbrdist in graph[node]: # v
            if dist[nbr] < dist[node] + nbrdist :
                dist[nbr] = dist[node] + nbrdist
                comesfrom[nbr] = node

    ### BACKTRACKING FOR MAXPATH
    maxpath = [endnode]
    while maxpath[-1] != startnode:
        maxpath.append(comesfrom[maxpath[-1]])
    maxpath.reverse()
    Length = len(maxpath)
    return dist[endnode], maxpath[1:Length-1]


def Convert(ALP_arr, Y, dir_, FName):
    End = len(Y[0]) - 1
    F = dir_ + FName
    with open(F,'w')as f:
        # f.write('Segmentation :' + '\n')
        # f.write("['Source',")
        for i in ALP_arr:
            f.write("'" + str(float(Y[int(i)][End - 1])) + '-' + str(float(Y[int(i)][End])) + "'" + ',')
        # f.write("'Target']" + '\n')
    f.close()


def Finding_Neigbours(YID, Y, node):
    Neigbours = []
    # YID.remove(node)
    End = len(Y[0]) -1
    E1 = Y[node][End]
    for item in YID:
        B2 = Y[item][End-1]
        if (E1 == B2):
            Neigbours.append(item)
    return Neigbours


def get_distance(node1, node2, P_Xtilde_Y):
    V1 = node1
    V2 = node2
    if len(P_Xtilde_Y) > 1:
        a = P_Xtilde_Y[:, V1]
        b = P_Xtilde_Y[:, V2]
    else:
        a = P_Xtilde_Y[0][V1]
        b = P_Xtilde_Y[0][V2]
    dist = 0
    for i in range(len(a)):
        dist = (a[i] - b[i])**2

    dist **= 0.5
    return dist


def Datapointnumber(End, Begin, Timestamp, EEnd):
    if End != EEnd:
        members = Timestamp[np.where((Timestamp >= Begin) & (Timestamp < End))]
    else:
        members = Timestamp[np.where((Timestamp >= Begin) & (Timestamp <= End))]
    number = len(members)
    return number


def get_quality(Distance, Min_n, Std, Mean, Total, method, fun, Frac):
    epsilon = 0.0001
    if method == 'Fraction':
        if fun == 'Sharp':
            if Min_n< (Frac* Total): #Mean - (2*Std):
                quality = min(Distance,epsilon)
            else:
                quality = Distance
        else:
            #print Min_n,Frac,Total
            try:
                quality = Distance / (1+math.exp(-Min_n+(Frac * Total)))
            except OverflowError:
                quality = 0.0
    else: # if method = 'Sigma':
        if fun == 'Sharp':
            if (Min_n< (Mean - (2*Std))): #:
                quality = min(Distance,epsilon)
            else:
                quality = Distance
        else:
            quality = (Distance / (1+math.exp(-Min_n+(Mean - (2*Std)))))
    return quality


def generate_graph(P_Xtilde_Y,Timestamp,Y,BeginTime,EndTime,method,fun,Frac):
    print('generate_graph...')
    G={}
    N_points_arr = {}
    End = len(Y[0])-1
    YID = range(0,len(Y))
    Source = 'Source' #len(Y)
    Target = 'Target' #len(Y) + 1
    G[Source] = []
    G[Target] = []
    TotalW = 0
    Totale = 0
    for node1 in YID:
        size1 = Y[node1][End]- Y[node1][End-1]
        number1 = Datapointnumber(Y[node1][End],Y[node1][End-1],Timestamp,EndTime)

        try:
            N_points_arr[size1].append(number1)
        except KeyError:
            N_points_arr[size1] = [number1]

        Neigbours = Finding_Neigbours(YID,Y,node1)
        ## Adding 2 extra nodes
        if (Y[node1][End] == EndTime) & (Y[node1][End-1] == BeginTime):
            continue
        if Y[node1][End] == EndTime:
            try:
                G[node1].append((Target,0,0,0,0))
            except KeyError:
                G[node1] = [(Target,0,0,0,0)]
        if Y[node1][End-1] == BeginTime:
            try:
                G[Source].append((node1,0,0,0,0))
            except KeyError:
                G[Source] = [(node1,0,0,0,0)]
        #################################################
        for node2 in Neigbours:
            size2 = Y[node2][End]- Y[node2][End-1]
            number2 = Datapointnumber(Y[node2][End],Y[node2][End-1],Timestamp,EndTime)


            if number1<=number2:
                Min_n = number1
                Min_s = size1
            else:
                Min_n = number2
                Min_s = size2
            Distance = get_distance(node1,node2,P_Xtilde_Y)
            quality = 0

            Weight = Distance
            TotalW += Weight
            Totale +=1
            try:
                G[node1].append((node2,Weight,Min_n,Min_s,quality))
            except KeyError:
                G[node1] = [(node2,Weight,Min_n,Min_s,quality)]
    Mean = {}
    Std = {}
    Graph = {}
    for Key in N_points_arr.keys():
        Mean[Key] = np.array(N_points_arr[Key]).mean()
        Std[Key] = np.array(N_points_arr[Key]).std()
    Total = Datapointnumber(EndTime,BeginTime,Timestamp,EndTime)
    for Key in G.keys():
        for i in range(0,len(G[Key])):
            node2,Weight,Min_n,Min_s,quality = G[Key][i]
            if Min_s> 0:
                quality = get_quality(Weight,Min_n,Std[Min_s],Mean[Min_s],Total,method,fun,Frac)
                G[Key][i] = (node2,Weight,Min_n,Min_s,quality)
    for Key in G.keys():
        for node2,Weight,Min_n,Min_s,quality in G[Key]:
            try:
                Graph[Key][node2] =  quality
            except KeyError:
                Graph[Key] = {}
                Graph[Key][node2] =  quality
    Graph['Target']={}
    return Graph,G,(float(TotalW)/Totale)


def Calc_Weights(G,maxpath):
    Weight = []
    for i in range(0,len(maxpath)-1):
        for (x,t1,t2,t3,y) in G[maxpath[i]]:
            if (x == maxpath[i+1]):
                Weight.append(y)
    return Weight


def ReadData(data_dir):
    print('ReadData')
    PFile = data_dir + 'P_Xtilde_Y.txt'
    P_Xtilde_Y = []
    with open(PFile) as f:
        Rows=f.read().splitlines()
    for row in Rows:
        #items=row.strip().split('\t')
        items=row.strip().split()
        tmp = [float(item) for item in items]
        P_Xtilde_Y.append(tmp)
    P_Xtilde_Y = np.array(P_Xtilde_Y)
    f.close()
    FYtilde =  data_dir + 'Ytilde.txt'
    Ytilde = []
    with open(FYtilde) as f:
        Rows=f.read().splitlines()
    for row in Rows:
        Ytilde.append(float(row))
    Ytilde = np.array(Ytilde)
    f.close()
    FXtilde = data_dir + 'Xtilde.txt'
    Xtilde = []
    with open(FXtilde) as f:
        Rows=f.read().splitlines()
    for row in Rows:
        Xtilde.append(float(row))
    f.close()
    FY = data_dir+'Y.txt'
    Y = []
    with open(FY) as f:
        Rows=f.read().splitlines()
    for row in Rows:
        items=row.strip().split('\t')
        tmp = [float(item) for item in items]
        Y.append(tmp)
    f.close()
    FTime = data_dir + 'BeginEndTime.txt'
    with open(FTime) as f:
        Rows = f.read().splitlines()
    for row in Rows:
        items=row.strip().split('\t')
        BeginTime = float(items[0])
        EndTime = float(items[1])
    f.close()

    FX = data_dir+'X.txt'
    X = []
    with open(FX) as f:
        Rows=f.read().splitlines()
    for row in Rows:
        items=row.strip().split('\t')
        tmp = [float(item) for item in items]
        X.append(tmp)
    X = np.array(X)
    f.close()
    End = len(X[0])-1
    if len(X)>1:
        Timestamp=X[:,End]
    else:
        Timestamp = X[0][End]

    return P_Xtilde_Y,Ytilde,Xtilde,Y,Timestamp,BeginTime,EndTime


def main(data_dir, ff, MinW, N_mdl, CalcMode, _dir, matlab_path, PathMode):
    aib(data_dir, MinW, N_mdl, CalcMode, _dir, matlab_path)
    start = time.clock()
    print('start time recording')
    s_lz=time.time()
    P_Xtilde_Y,Ytilde,Xtilde,Y,Timestamp,BeginTime,EndTime = ReadData(data_dir)
    l = int((EndTime - BeginTime) / MinW + 3)
    G, G1, TotalAVG = generate_graph(P_Xtilde_Y, Timestamp, Y, BeginTime, EndTime, 'Fraction', 'Sigmoid', ff)
    dir_ = data_dir
    directory = os.path.dirname(dir_)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if PathMode == 'ALP':
        sf = open(data_dir + 'ALP_runnimgtime.txt', 'w')
        start1 = time.clock()
        start2 = time.time()
        ALP_arr, avg, ALP_length = falp.main(G, 'Source', 'Target', l)
        end2 = time.time()
        end1 = time.clock()
        e_lz=time.time()

        sf.write(data_dir + '\n')
        sf.write(str(e_lz-s_lz)+'\n')
        sf.close()
        # print('fast alp cpu time: ' + str(end1 - start1))
        # print('fast alp time: ' + str(end2 - start2))
        # print('ALP_arr:' + str(ALP_arr))
        # print('avg: ' + str(avg))
        # print('ALP_length' + str(ALP_length))

        # FName = 'Segmentation.txt'
        # Name = 'ALP.txt'
        # CostName = 'ALP-AVG.txt'

    # elif PathMode == 'LP':
    #     avg, ALP_arr = longestpathDAG(G1, 'Source', 'Target')
    #     print str(ALP_arr)
    #     FName = 'LP2.txt'
    #     Name = 'LP.txt'
    #     CostName = 'LP-AVG.txt'
    #     avg = float(avg) / len(ALP_arr)

    FName = 'Segmentation.txt'
    Convert(ALP_arr, Y, dir_, FName)

    # with open(dir_ + Name, 'wb')as f:
    #     for j in ALP_arr:
    #         f.write(str(j) + '\t')
    # with open(dir_ + CostName, 'wb') as g:
    #     g.write('AVG of path: ' + str(avg) + '\n')
    #     g.write('Totla AVG: ' + str(TotalAVG))
    elapsed = (time.clock() - start)
    print (data_dir + ':' + str(elapsed))


if __name__ == '__main__':
    data_dir = sys.argv[1] #'./data/test/' # Dirrectory of the data
    matlab_path = sys.argv[2] #'/Applications/MATLAB_R2016a.app/bin/matlab'
    data_dir += '/'
    data_dir = data_dir.replace('//', '/')
    _dir = data_dir + 'input.txt'
    ff = 0.001
    MinW = int(sys.argv[3]) #2*86400 #604800
    N_mdl = int(sys.argv[4]) #5
    CalcMode = int(sys.argv[5]) #0
    PathMode = 'ALP'
    main(data_dir, ff, MinW, N_mdl, CalcMode, _dir, matlab_path, PathMode)
