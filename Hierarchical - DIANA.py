
# coding: utf-8

# In[258]:


#IMPORT REQUIRED LIBRARIES
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import random
import copy
from copy import deepcopy
import numpy as np
import pandas as pd


# In[259]:


#READ RECORDS FROM AMINO.FASTA
records = list(SeqIO.parse("amino.fasta", "fasta"))


# In[260]:


#READ DISTANCES from CSV
data = pd.read_csv('distances.csv')
#Add all data in a list of lists
dist = []  
for i in range(0, len(data)):  
    dist.append([data.values[i,j] for j in range(0, len(data.values[0]))])


# In[261]:


#INITIALISE LIST OF CLUSTER WITH A SINGLE CLUSTERING CONTAINING ALL POINTS
cluster=[]
a=set()
for i in range(len(records)):
    a.add(i)
cluster.append(a)


# In[262]:


#TO STORE POINTS FOR DENDOGRAM
xlist=[]


# In[263]:


#CREATING COPY OF CLUSTER
newclust=(copy.deepcopy(cluster))


# In[264]:


#FUNCTION TO BREAK EACH CLUSTER INTO TWO
x=[]
y=[]

pointx=[i for i in range(0,len(dist))]
pointy=[0 for i in range(0,len(dist))]
def breakclust(clust):
    if(len(clust)==1):
        return
    child1=set()
    child2=set()
    #IDENTIFY MOST DISSIMILAR POINTS
    max=-99999
    for i in clust:
        for j in clust:
            if(i!=j):
                if(dist[i][j]!=0):
                    if(dist[i][j]>max):
                        if(i!=j):
                            max=dist[i][j]
                            maxi=i
                            maxj=j
    #FORM NEW CLUSTERS WITH MEDIODS AS THE ABOVE POINTS
    child1.add(maxi)
    child2.add(maxj)
    clust.remove(maxi)
    clust.remove(maxj)
    if(len(clust)!=0):
        for i in clust:
            min1=1
            for j in child1:
                if(dist[i][j]<min1):
                    min1=dist[i][j]
            min2=1
            for j in child2:
                if(dist[i][j]<min2):
                    min2=dist[i][j]
            if(min1<min2):
                child1.add(i)
            else:
                child2.add(i)
    
    cluster.append(child1)
    cluster.append(child2)  
    #BREAK CLUSTER INTO TWO BY CALLING BREAKCLUST FOR EACH CHILD RECURSIVELY
    breakclust(copy.deepcopy(child1))   
    breakclust(copy.deepcopy(child2))

    #STORING COORDINATES FOR DENDOGRAMS
    midpoint=(pointx[maxi]+pointx[maxj])/2
    x.append(pointx[maxi])
    x.append(pointx[maxi])
    x.append(pointx[maxj])
    x.append(pointx[maxj])

    y.append(pointy[maxi])
    y.append(1/abs(max))
    y.append(1/abs(max))
    y.append(pointy[maxj])

    pointx[maxi]=midpoint
    pointx[maxj]=midpoint
    pointy[maxi]=1/abs(max)
    pointy[maxj]=1/abs(max)


# In[265]:


#CALL FUNCTION FOR FIRST CLUSTER
breakclust(copy.deepcopy(cluster[0]))


# In[266]:


#PLOTTING THE DENDOGRAM
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure
figure(num=None, figsize=(17, 17), dpi=300)
count=350
loop=int(len(x)/4)
x[len(x)-1]=(x[len(x)-7]+x[len(x)-6])/2
x[len(x)-2]=(x[len(x)-7]+x[len(x)-6])/2
y[len(y)-1]=(y[len(y)-6])
for i in range(loop):
    xlist=[]
    ylist=[]
    for j in range(0,4):
        xlist.append(x[j+i*4])
        ylist.append(y[j+i*4]) 
    plt.plot(xlist, ylist)

