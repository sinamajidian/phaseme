#!/usr/bin/env python3
#import networkx as nx
#from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
#import ast
#from sys import argv
import itertools




#if __name__ == "__main__":
hap_name = 'haplotype_samples.txt'
hap_all=np.zeros(89) # initaial row
#with open(hap_name, 'r') as file_haps:
file_haps=open(hap_name, 'r')
for line in file_haps:
    line_np=np.array(line.split())
    hap_all=np.vstack([hap_all, line_np.astype(int)])
hap_all=hap_all[:,1:] # removing initaial row
haplolength, nsample=np.shape(hap_all)  
hap_first= hap_all[:,0] #np.zeros(haplolength)
number_flipped=0
for i in range(1,nsample): # expecpt first one
    hap_i=hap_all[:,i]
    hap_all_ordered=hap_all;
    index_best_match=np.argmin([np.sum(hap_i==hap_first),np.sum(hap_i!=hap_first)])
    if index_best_match: # change it to its flip
        hap_i_flip=~hap_i.astype(bool)
        hap_all_ordered[:,i]= hap_i_flip.astype(int)
        number_flipped+=1

print('From',nsample+1,'samples, ',number_flipped,'were flipped!')    
    
    
hap_all_ordered[1:5,1:5]




k=8;
#start=0
#hap_first=hap_all[start:k+start,0] #np.zeros(haplolength)
#possible=np.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]])

possible = np.array(list(itertools.product([0, 1], repeat=k)))

freq=np.zeros([haplolength-k,2**k])
for start in range(1,3000): # haplolength-k

    for i in range(1,nsample): # expecpt first one
        for ii in range(2**k):
            if np.array_equal(possible[ii,:],hap_all_ordered[start:k+start,i]):
                freq[start,ii]+=1
#print(freq) 
freq_max=np.max(freq,axis=1)/nsample

plt.plot(freq_max[1:3000])
plt.ylabel('Max frequency of haplotype of length eight')
plt.xlabel('SNP index')
plt.show()

