import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
import pandas as pd
import os
from IPython.display import display
from IPython.display import HTML
import dateutil;
import json;
import urllib;

import matplotlib
matplotlib.style.use('ggplot')

def computePCA(d, components = 2, columnNames=None):

    if (d.shape[1] < 2) :
        print ("Number of components are already less than 2")

    pca = PCA(n_components= components)
    pca.fit(d)

    #usorted=  zip(df.columns.values, pca.components_[0])
    #sort = sorted (usorted, key = lambda x:x[1])
    #print (sort)

    return pca
#======================================================================
'''
Evaluating the cost of the clusters
'''
'''
The costOfCluster is the cost assigned to a k-means cluster.
K-means algorithm optimizes on the cost function SUM( ||xi - Mj||^2) - where Mj is the
centroid of of cluster. In this case, if xi ( i is from 1..m) is assigned to cluster j ( j from 1..k)
then the ||xi - Mj||^2) is the square of the distaance
'''
#from scipy.spatial.distance import *
#from sklearn import metrics
#from sklearn.metrics import pairwise_distances

def costOfCluster(kmeans, data) :
    l = kmeans.labels_
    c = kmeans.cluster_centers_
    s = 0.0
    for i, x in enumerate(data):
        c1 = euclidean (x, c[l[i]]);
        s =  s + c1 ** 2
    return s;

# x must be an np array. Thi
def elbowIndex(x):
    if ( x.shape[0] <= 2):
        return 0
    d = 0
    d = x[:-1] - x[1:]
    second_d = np.abs(d[:-1] - d[1:])
    if (second_d.shape[0] <=0):
        return 0;
    print ("+++", second_d.shape)
    return 1 + np.argmax(second_d)
