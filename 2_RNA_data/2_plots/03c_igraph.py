#!/home/lauber/anaconda3/bin/python

# load modules
from igraph import *


# parameters
file_edges    = "edges.csv"
file_vertices = "nodes.csv"
n_vertices    = 70534
out_fig_name  = "igraph.eps"
layout_iter   = 10000
dimensions    = (10000,10000)

# Create graph
g = Graph()

# Add vertices
g.add_vertices(n_vertices)

edges = []
weights = []

with open("edges.csv", "r") as edges_file:
    line = edges_file.readline()
    
    while line != "":
        
        strings = line.rstrip().split(",")
        
        # Add edge to edge list
        edges.append(((int(strings[0])-1), (int(strings[1])-1)))
        
        # Add weight to weight list
        weights.append(float(strings[2]))
        
        
        line = edges_file.readline()

# Add edges to the graph
g.add_edges(edges)

# Add weights to edges in the graph
g.es['weight'] = weights



# Visualize the graph according to labels

n_classes = 6

bins = [[] for x in range(n_classes)]

with open("nodes.csv", "r") as labels_file:
    line = labels_file.readline()
    
    while line != "":
        
        strings = line.rstrip().split(",")
        
        vertex_id = int(strings[0])-1
        bin_id = int(strings[1])-1
        bins[bin_id].append(vertex_id)
        
        line = labels_file.readline()


# style the graph

node_colours = []

for i in range(n_vertices):
    if   i in bins[0]:
        node_colours.append("skyblue")
    elif i in bins[1]:
        node_colours.append("hotpink")
    elif i in bins[2]:
        node_colours.append("red")
    elif i in bins[3]:
        node_colours.append("blue")
    elif i in bins[4]:
        node_colours.append("green")
    elif i in bins[5]:
        node_colours.append("orange")
    else:
        node_colours.append("grey")

g.vs["color"]       = node_colours
g.vs["frame_color"] = node_colours


# layout and style
visual_style = {}

# the layout
#my_layout = g.layout_fruchterman_reingold( niter=1000 ) #, start_temp=LFR_stemp )
#my_layout = g.layout_kamada_kawai()
my_layout = g.layout_graphopt( niter=layout_iter )
visual_style["layout"] = my_layout

# Define colors used for outdegree visualization
colours = ['#fecc5c', '#a31a1c']

# Set bbox and margin
visual_style["bbox"] = dimensions #(5000,5000) #(3000,3000)
visual_style["margin"] = 17

# Set vertex size
visual_style["vertex_size"] = 10 #20

# Set vertex lable size
visual_style["vertex_label_size"] = 8

# Don't curve the edges
visual_style["edge_curved"] = False

# Plot the graph
plot(g, out_fig_name, **visual_style)

