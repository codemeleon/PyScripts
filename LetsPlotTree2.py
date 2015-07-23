# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 21:20:04 2015

@author: devil
"""
import matplotlib
matplotlib.use("Agg")
import colorsys
import pandas as pd
import numpy as np
from ete2 import Tree
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches


rec_path = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
            ]
def rectangle(x,y,w,h,col):
   #print x,y,w,h,col
    verts = [
            (x, y), # left, bottom
            (x, y+h), # left, top
            (x+w, y+h), # right, top
            (x+w, y), # right, bottom
            (x, y), # ignored
            ]
    path = Path(verts, rec_path)
    patch = patches.PathPatch(path, facecolor=col,lw=0.00,alpha=0.25)# edgecolor='orange' if col=='red' else 'white',
    return patch



def create_colors(N):
    HSV_tuples = [( x*1.0/N, 0.5, 0.5) for x in range( N)]
    RGB_tuples = map( lambda x: colorsys.hsv_to_rgb( *x), HSV_tuples)


allcolor = 'rgbcmyk'*2

def tree_node_numbering(tree):
    """Numbering the node which doesn't have a name"""
    j = 0
    for node in tree.traverse("preorder"):
        if node.name == 'NoName' or not node.is_leaf():
            node.name = str(j)
            j +=1
    return tree


def node_arrangements(tree):
    leafs = tree.get_leaf_names()
    nleafs = len(leafs)
    nodes = list(leafs)
    leaf_y_pos = np.linspace(1.0,0.0,nleafs)
    leaf_half_dist = (leaf_y_pos[0]-leaf_y_pos[1])/2.0
    Y = {leafs[j]:leaf_y_pos[j]+leaf_half_dist for j in xrange(nleafs)}
    max_Y = max(Y.values( ))+leaf_half_dist
    Y = {k:Y[k]/float(max_Y) for k in Y.keys( )}

    def ycoors(tree):   # Calculate Y_distances of nodes from the bottom
        for child in tree.get_children():
            if child.name not in nodes:
                ycoors(child)
        Y[tree.name] = (Y[tree.get_children()[0].name] + Y[tree.get_children()[1].name])/2.0
        nodes.append(tree.name)

    ycoors(tree)


    X = {}  # X_distances of nodes from root
    for node in tree.traverse("preorder"):
        X[node.name] = tree.get_distance(node.name)
    max_X = max(X.values( ))
    X = {k:X[k]/float(max_X) for k in X.keys( )}

    return X,Y

def plot_tree(tree,X,Y,plt):
    for node in tree.traverse("preorder"):
        root = 0
        if root and node.is_root( ):
            plt.plot([X[node.name],X[node.name]-10],[Y[node.name],Y[node.name]],color='k')
        if not node.is_leaf( ):
            for child in node.get_children( ):
                plt.plot([X[node.name],X[child.name]],[Y[child.name],Y[child.name]],color='k')
                plt.plot([X[node.name],X[node.name]],[Y[node.name],Y[child.name]],color="k")
                node_num = 0
                if node_num:
                    plt.text(X[node.name],Y[node.name],node.name,fontsize=5,color="k",ha='left',va='center')

        leaf_names = 0
        if leaf_names:
            alt_name = 0
            if alt_name:
                altname = pd.read_csv("")
            for leaf in tree.get_leaf_names( ):
                plt.text(X[leaf],Y[leaf],leaf,fontsize=5,color="k",ha='left',va='center')
        if not leaf_names:
            max_X = max(X.values( ))
            Y_leaf_dist = (max(Y.values( ))-min(Y.values( )))/float(len(tree.get_leaf_names( )))
            for leaf in tree.get_leaf_names( ):
                line, = plt.plot([X[leaf],max_X],[Y[leaf],Y[leaf]],color='k')
                line.set_dashes([1,1,1,1])
                continue
                rec = rectangle(max_X,Y[leaf]-Y_leaf_dist/2.0,10,Y_leaf_dist,col='r')##Correct width with pixel ammount
                plt.gca( ).add_patch(rec)
            ""
    return plt

def blocks(tree,data,plt,next_X):
    data['Names'] = map(str,data['Names'])
    index = list(data['Names'])
    data = data.set_index(data['Names'])
    del data['Names']
    ""
def blocks_with_names(tree,data,plt,next_X):
    data["Names"] = map(str,data["Names"])
    index = list(data["Names"])
    #print data
    data = data.set_index(data["Names"])
    del data["Names"]
    #print data
    Y_leaf_dist = (max(Y.values( ))-min(Y.values( )))/float(len(tree.get_leaf_names( )))
    X_w = 30
    col_comb = {}
    for col in data.columns:
        ## find a way to add more color
        for i,d in  enumerate(set(data[col])):
     #       print i,d
            col_comb[d] = allcolor[i]

    for leaf in tree.get_leaf_names( ):
        if leaf not in index:continue
        t_X = next_X
        for column in data.ix[leaf]:
            #plt.text(t_X,Y[leaf],col,fontsize=5,color='k',ha='left',va='center')
            rec = rectangle(t_X,Y[leaf]-Y_leaf_dist/2.0,X_w,Y_leaf_dist,col=col_comb[column])
            plt.gca().add_patch(rec)
            t_X += X_w
    t_X = next_X
    t_Y = max(Y.values( ))+Y_leaf_dist/2.0
    for col in data.columns:
        plt.text(t_X+X_w/2.0,t_Y,col,fontsize=5,color='k',ha='center',va='bottom',rotation='vertical')
        t_X +=X_w
    return plt


if __name__=='__main__':
    height = 8
    width = 11
    fig = plt.Figure(figsize =(11,8),dpi=300)
    height,width = height*300,width*300
    input_files = ["/home/devil/Documents/Jen/GAS/Final Data/Vir_Emm.csv","/home/devil/Documents/Jen/GAS/Final Data/test_antibiotics.csv"]
    tree_file ="../RAxML_bestTree.hypervirulent"
    tree_file = "/home/devil/Documents/Jen/GAS/Final Data/RAxML_bestTree.GAS_L.rax"
    tree =Tree(tree_file,format=1)
    tree.set_outgroup(tree&'NP')
    tree.ladderize(direction=0)
    tree = (tree)
    tree = tree_node_numbering(tree)
    X,Y = node_arrangements(tree)
    X = {k:X[k]*width for k in X.keys( )}
    Y = {k:Y[k]*height for k in Y.keys( )}
    print max(X.values( ))
    plt = plot_tree(tree,X,Y,plt)
    next_X  = max(X.values( )) + 20
    #top margin of 10 char
    for fl in input_files:
        #continue
        data = pd.read_csv(fl)
        if len(data.columns) < 10:
            plt = blocks_with_names(tree,data,plt,next_X)
            next_X += 30*len(data.columns) + 20
        else:
            plt = blocks(tree,data,plt,next_X)
            next_X += 5*len(data.columns)+20
    plt.xlim(0,3400)

    plt.axis("off")
    plt.savefig("Anmol.pdf",bbox_inches="tight",dpi=300)



    #print "Not yet Done",a

