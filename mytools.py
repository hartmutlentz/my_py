#!/usr/bin/env python3.4
#
#  
#
#
import math, random, sys
import numpy, scipy, scipy.sparse
import networkx as nx

from statsmodels.distributions.empirical_distribution import ECDF

assert sys.version_info >= (3,4)



def hitcode_bundesland(revert=False):
    """ dict with HIT codes for Bundesland
        HIT Code: 276+01+Nummer=D+SH+Nummer
    """
    d={}
    d["01"]="SH" # Schleswig-H.
    d["02"]="HH"
    d["03"]="NI"
    d["04"]="HB"
    d["05"]="NW"
    d["06"]="HE"
    d["07"]="RP"
    d["08"]="BW"
    d["09"]="BY"
    d["10"]="SL"
    d["11"]="BE"
    d["12"]="BB"
    d["13"]="MV"
    d["14"]="SN" # Sachsen
    d["15"]="ST" # Sachsen-Anhalt
    d["16"]="TH"
    
    if revert: return revert_dictionary(d,True)
    else: return d
    

def histogram(seq):
    dmax=max(seq)+1
    freq=[0 for d in range(dmax)]
    for d in seq:
        freq[d] += 1
    return freq


def giant_component(G, strongly=True):
    """ returns the giant component of a network as a (Di)Graph.
        If network is directed:
            strongly=True, returns GSCC
            GWCC otherwise.
    """
    if G.is_directed():
        if strongly:
            components = nx.strongly_connected_component_subgraphs
        else:
            components = nx.weakly_connected_component_subgraphs
    else:
        components = nx.connected_component_subgraphs
    
    ccs = components(G)
    ccs = sorted(ccs, key=len, reverse=True)
    return ccs[0]
    

def giant_out_component(G, return_giant_component=True):
    """ Returns the giant out component of a directed network.
        If return_giant_component == False, only GOC nodes without GC nodes
        are returned.
    """
    GC = giant_component(G, strongly=True)
    
    # start at a random node and do BFS
    node = GC.nodes()[0]
    GOC_nodes = list(nx.dfs_preorder_nodes(G, node))
    
    if not return_giant_component:
        GOC_nodes = set(GOC_nodes) - set(GC.nodes())

    return G.subgraph(GOC_nodes)

def giant_in_component(G, return_giant_component=True):
    """ Returns the giant in component of a directed network.
        If return_giant_component == False, only GIC nodes without GC nodes
        are returned.
    """
    return giant_out_component(G.reverse(), return_giant_component)


def ranges_percolating_system(G):
    """ computes Ranges and uses giant connected component
        to save cpu time.
    
    """
    # print 'Ranges for percolating network is not tested yet!'
    
    if G.is_directed():
        #comps = nx.strongly_connected_component_subgraphs(G)
        #comps = list(comps)
        #comps.sort(compare, reverse=True)
        comps = sorted(nx.strongly_connected_component_subgraphs(G),\
            key=len, reverse=True)
        LCC = comps[0]
    else:
        comps = sorted(nx.connected_component_subgraphs(G),\
            key=len, reverse=True)
        LCC = comps[0]

    rang = {}

    # an arbitrary LCC node
    laenge = nx.single_source_shortest_path_length(G,LCC.nodes()[0])
    lcc_range = len(laenge)-1
    # = all LCC nodes
    for node in LCC.nodes():
        rang[node] = lcc_range
    print("ranges: LCC done.")

    # remaining nodes
    nodes = set(G.nodes()) - set(LCC.nodes())
    remaining_nodes = len(nodes)
    for (i, start) in enumerate(nodes):
        print("Node ", i, " of ", remaining_nodes)
        laenge = nx.single_source_shortest_path_length(G,start)
        rang[start] = len(laenge) - 1
    return rang


def ranges(G,nodes=False,display=False):
    return ranges_single_sources(G,nodes,display)

def ranges_single_sources(G,nodes=False,display=False):
    """
    Space-saving version of ranges_small_graphs.
    Returns dict.
    """
    if nodes==False:
        nodes=G.nodes()
    if display:
        i=0
        no=G.number_of_nodes()
        rang={}
        for start in nodes:
            laenge=nx.single_source_shortest_path_length(G,start)
            rang[start]=len(laenge)-1
            print(i,' von ',no)
            i=i+1
    else:
        rang={}
        for start in nodes:
            laenge=nx.single_source_shortest_path_length(G,start)
            rang[start]=len(laenge)-1
    return rang

def node_range(G,node=False):
    if node:
        return range_single_node(G,node)
    else:
        return ranges_single_sources(G)

def range_single_node(G,node):
    """
    Space-saving version of ranges_small_graphs.
    Returns dict. *Use rang.get(i,' ') to avoid key errors*
    """
    laenge=single_source_shortest_path_length(G,node)
    rang=len(laenge)-1

    return rang    

def reachabilities(G,nodes):
    return reachabilities_single_sources(G,nodes)

def reachabilities_single_sources(G,nodes=None):
    """
    Number of nodes that can reach a node i
    Returns dict. 
    """
    if nodes==None:
        the_nodes=G.nodes()
    G1=G.reverse()
    return ranges_single_sources(G1,the_nodes)

def revert_dictionary(di,bijection=False):
    """
    Reverts dictionary. Output as dict of lists
    Input: {1:2, 5:2, 3:8,...}
    Output: {2:[1,5], 8:[3]}
    """
    rev={}
    if bijection:
        for x in di:
            rev[di[x]]=x
        if len(rev)!=len(di): 
            print('Dictionary not bijective! Cannot revert.')
            return None
    else:
        for x in di:
            if di[x] in rev: rev[di[x]].append(x)
            else: rev[di[x]]=[x]
    
    return rev
    
def sort_list_by_length(lis):
    """
    returns list sorted by the lengths of entries
    [[1,2],[1],[3,4,5,6]] -> [[3,4,5,6],[1,2],[1]]
    """
    def cmp(a,b):
        return len(a)-len(b)
    lis.sort(cmp,reverse=True)

def dict2file(d,nameoffile='dict.txt',sorted=True):
    """ Writes dictionary (or list or tuple) to a textfile
        Sorted by keys, if sorted=True.
    """
    def list2dict(li):
        x={}
        for (i,el) in enumerate(li):
            x[i]=el
        return x
    
    if not isinstance(d,dict): d=list2dict(d)
    
    dk=list(d.keys())
    if sorted: dk.sort()
    
    # if d={ 1: [a,b,c,...], 2:[d,e,f,...],... }
    s=list(d.values())[0]
    if isinstance(s,dict) or isinstance(s,list) or isinstance(s,tuple) or isinstance(s,numpy.ndarray):
        laenge=len(list(d.values())[0])
        g=file(nameoffile,'w+')
        for k in dk:
            wstring=''
            for l in range(laenge): wstring += '\t'+str(d[k][l])
            g.writelines(( str(k)+wstring+'\n' ))
        g.close
        return
    
    g=file(nameoffile,'w+')
    for k in dk:
        g.writelines((str(k),'\t',str(d[k]),'\n'))
    g.close
    
    return

def cdf(seq):
	"""
	The CDF of a sequence.
	"""
	ecdf = ECDF(seq)
	return 1 - ecdf(range(int(max(seq))))
	


if __name__=="__main__":
    c = hitcode_bundesland(True)
    print(c)

