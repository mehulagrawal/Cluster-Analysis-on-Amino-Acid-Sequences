
# coding: utf-8

# In[15]:


#IMPORT REQUIRED LIBRARIES
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import random
import copy
from copy import deepcopy
import pandas as pd


# In[16]:


#READ RECORDS FROM AMINO.FASTA
records = list(SeqIO.parse("amino.fasta", "fasta"))


# In[17]:


#Read DISTANCES from CSV
data = pd.read_csv('distances.csv')
#Add all data in a list of lists
dist = []  
for i in range(0, len(data)):  
    dist.append([data.values[i,j] for j in range(0, len(data.values[0]))])


# In[18]:


#FUNCTION TO UPDATE MEDIAN
def updatemed(i,clust):
    min=99999999
    pos=0
    cost=0
    for k in clust:
        if(clust[k]==i):
            cost=0
            for j in clust:
                if(clust[j]==i):
                    cost=cost+dist[j][k]
        if(cost<min):
            min=cost
            pos=k
    return pos


# In[19]:


#CALCULATE TOTAL SSE
def totcost(k,clust):
    cost=0
    for i in range(0,len(k)):
        for j in clust:
            if(clust[j]==i):
                cost=cost+dist[j][i]
    return cost


# In[26]:


#MAIN FUNCTION TO PERFORM K-MEANS
mincost=99999999
clustfinal={}
for i in range(0,10):
    k=random.sample(range(0, 236),5)
    clust={}
    p=1
    while(p==1):
        for i in range(0,236):
            min=999999
            for j in range(0,len(k)):
                d=dist[i][j]
                if(d<min):
                    min=d
                    clust[i]=j
        k2=copy.deepcopy(k)
        for i in range(0,len(k)):
            k[i]=updatemed(i,clust)
        if(k==k2):
            break
    costf=totcost(k,clust)
    if(mincost>costf):
        mincost=costf
        clustfinal=copy.deepcopy(clust)


# In[27]:


#PRINT CLUSTER
for j in range(0,len(k)):
    print(j+1,"th Cluster")
    for i in clust:
        if(clustfinal[i]==j):
            print(i,end=", ")
    print("\n")

