
# coding: utf-8

# In[163]:


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


# In[164]:


#CHOOSE LINKAGE
choice=3
#1 for min
#2 for max
#3 for avg
#READ RECORDS FROM AMINO.FASTA
records = list(SeqIO.parse("amino.fasta", "fasta"))


# In[165]:


#Read DISTANCES from CSV
data = pd.read_csv('distances.csv')
#Add all data in a list of lists
dist = []  
for i in range(0, len(data)):  
    dist.append([data.values[i,j] for j in range(0, len(data.values[0]))])
distc=copy.deepcopy(dist)


# In[166]:


#INITIALLY TAKE EACH SINGLE DATA POINT IN A CLUSTER
cluster=[set() for i in range(len(records))]
dist2=[i for i in range(len(records))]
for i in range(0,len(records)):
    cluster[i].add(i)


# In[167]:


#MIN-LINK
def mergemin(minj,index):
#INDEX IS POSITION OF NEW CLUSTER IN CLUST
    dist[minj][minj]=-99999
    for i in range(0,len(records)):
        min=9999
        for k in cluster[index]:
            if(dist[i][k]<min):
                if(i!=k):
                    min=dist[i][k]
        dist[minj][i]=dist[i][minj]=min

#MAX-LINK        
def mergemax(minj,index):
#INDEX IS POSITION OF NEW CLUSTER IN CLUST
    dist[minj][minj]=-99999
    max=-9999
    for i in range(0,len(records)):
        max=-9999
        for k in cluster[index]:
            if(dist[i][k]>max):
                if(i!=k):
                    if(dist[i][k]!=1):
                        max=dist[i][k]
        if(max!=-9999):
            dist[minj][i]=dist[i][minj]=max

# #AVERAGE-LINK            
# def mergeavg(minj,index,l1,l2):
#     #INDEX IS POSITION OF NEW CLUSTER IN CLUST
#     dist[minj][minj]=-99999
#     sum=0
#     for i in range(0,len(records)):
#         max=-9999
#         for k in cluster[index]:
#             if(i!=k):
#                 if(dist[i][k]!=1):
#                     sum=sum+dist[i][k]
#         if(sum!=0):
#             dist[minj][i]=dist[i][minj]=(sum/(l1*l2))
def mergeavg(minj,index,l1,l2):
    #INDEX IS POSITION OF NEW CLUSTER IN CLUST
    dist[minj][minj]=-99999
#     sum=0
    for i in range(0,len(records)):
        sum=0
        min=9999999999
        for k in cluster[index]:
            if(i!=k):
                if(dist[i][k]!=1):
                    sum=sum+dist[i][k]
        if(sum/(l1*l2)<min):
            min=sum/(l1*l2)
        dist[minj][i]=dist[i][minj]=sum/(l1*l2)
    print(min)


# In[168]:


#CALCULATE DISTANCES FOR DENDOGRAM
b=1
x=[]
y=[]
pointx=[i for i in range(0,len(dist))]
pointy=[0 for i in range(0,len(dist))]
while(b==1):
    #COMPUTE MINIMUM VALUE IN MATRIX
    min=9999
    for i in range(0,len(dist[0])):
        for j in range(0,len(dist[0])):
            if(dist[i][j]<min):
                if(i!=j):
                    min=dist[i][j]
                    mini=i
                    minj=j
    #CREATE CLUSTER ITEM
    
    a=set()
    l1=len(cluster[dist2[minj]])
    l2=len(cluster[dist2[mini]])
    for i in cluster[dist2[minj]]:
        a.add(i)
    for i in cluster[dist2[mini]]:
        a.add(i)
    if (min==0):
        break
    cluster.append(a)
#     print(len(a))
    h=(pointx[mini]+pointx[minj])/2
    x.append(pointx[minj])
    x.append(pointx[minj])
    x.append(pointx[mini])
    x.append(pointx[mini])

    y.append(pointy[minj])
    y.append(1/-min)
    y.append(1/-min)
    y.append(pointy[mini])

    pointx[mini]=h
    pointx[minj]=h
    pointy[mini]=1/-min
    pointy[minj]=1/-min
    
    print(len(a)," ")
    if(len(a)==len(records)):
        break
    #SET CLUSTER POINTER
    dist2[minj]=len(cluster)-1
    if(choice==1):
        mergemin(minj,dist2[minj])
    if(choice==2):
        mergemax(minj,dist2[minj])
    if(choice==3):
        mergeavg(minj,dist2[minj],l1,l2)
    # SET COLUMN AND ROW =1
    for i in range(0, len(records)):
        dist[mini][i]=1
        dist[i][mini]=1


# In[169]:


#COPY DISTANCE MATRIX AGAIN BECASUE ORIGINAL IS MODIFIED
dist=copy.deepcopy(distc)


# In[170]:


#PLOTTING THE DENDOGRAM
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure
figure(num=None, figsize=(15, 15), dpi=100)
# x=[0,0,1,1,0.5,0.5,2,2]
# y=[0,1,1,0,1,2,2,0]
count=0
loop=int(len(x)/4)
for i in range(0,loop):
    xlist=[]
    ylist=[]
    for j in range(0,4):
        xlist.append(x[j+i*4])
        ylist.append(y[j+i*4])
        count=count+1
    plt.plot(xlist, ylist) 


# In[ ]:



            

