import os
import sys
from igraph import *

graphml_path = "/graphml"
graphmat_path = "/mat"


for graph_file in os.listdir(graphml_path):
    if graph_file[-8:] == ".graphml":
        graph = Graph()
        graph = graph.Read_GraphML(graphml_path + "/" + graph_file)
        
        graph.write_adjacency(graphmat_path + "/" + graph_file[:-8] + ".mat")

