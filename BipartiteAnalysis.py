#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:18:27 2019

@author: x
"""

# -*- coding: utf-8 -*-
"""
Script to take an R igraph-encoded graphml of a bipartite graph, calculate 
"""

import networkx as nx
import argparse
parser =  argparse.ArgumentParser()
parser.add_argument("graph")
parser.add_argument("out")
args = parser.parse_args()
path = args.graph
path_out = args.out

B = nx.read_graphml(path)

nx.set_node_attributes(B,
                       values = dict(nx.degree(B)),
                       name = "grado" 
                       )

#F = B.subgraph((n for n in B if B.node[n]['grado']>1))
F = B
F = F.to_undirected()

print(F)
print(nx.is_bipartite(F))
#clustering coefficients

nx.set_node_attributes(F, 
                       values = nx.bipartite.clustering(F, mode = "dot"),
                       name = "ClusteringDot"
                       )


nx.set_node_attributes(F, 
                       values = nx.bipartite.robins_alexander_clustering(F),
                       name   = "ClusteringRobbinsAlexander"
                       )


#if redundancy fails with  Cannot compute redundancy coefficient for a node that has fewer than two neighbors
redundancia= dict()
for nodo in F.nodes():
    if F.nodes[nodo]["grado"]<3:
        redundancia[nodo] = 0
    else:
        redundancia[nodo] = nx.bipartite.node_redundancy(F, nodes = [nodo])[nodo]

        
nx.set_node_attributes(F, 
                       values = redundancia,
                       name ="Redundancy" 
                       )

nx.write_graphml(G= F, infer_numeric_types= True, path= path_out, prettyprint=True)
#nx.write_graphml(F, path_out)
print("victoria!")
