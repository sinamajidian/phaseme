



"""
Credit: https://github.com/vibansal/IntegratedPhasing/blob/master/encodereads.py 

Vikas Bansal, Integrating read-based and population-based phasing for dense and accurate haplotyping of individual genomes, Bioinformatics 2019.

Edited by Sina Majidian


"""



import os,sys,subprocess,random,math
#from mstkruskal import Graph








def pairwiseLD(haplist,nsamples,i,j):
	counts = [0,0];
	for r in range(nsamples): 
		if haplist[r][i] == '0' and haplist[r][j] == '0': counts[0] +=1; 
		elif haplist[r][i] == '1' and haplist[r][j] == '1': counts[0] +=1; 
		elif haplist[r][i] == '0' and haplist[r][j] == '1': counts[1] +=1; 
		elif haplist[r][i] == '1' and haplist[r][j] == '0': counts[1] +=1; 
	return counts;



HAPfile = sys.argv[1] # input 'ont_500/ont.hapsamples' #




File = open(HAPfile,'r'); 
nsamples=0; haplist = []
for line in File: 
    if nsamples ==0: ## first line in new format has variant indexes
        varindex = list(map(int, line.strip().split()));  
    else: 
        haplist.append(list(line.strip())); 
    nsamples +=1; 
File.close(); 

nsamples -=1; 



file_report= open('pairs.txt','w');  # output



THRESH=0.98

nsamples = len(haplist); 
v = len(haplist[0]); 
edges =[];
W = 10; #20; 5
size = v;
discarded=0;
for i in range(size):
    for j in range(i+1,min(v,i+W)):
        counts = pairwiseLD(haplist,nsamples,i,j);
        f = float(counts[0]+0.5)/(counts[0]+counts[1]+1);
        if f > THRESH: 
            edges.append([i,j,f,'00']);  
        elif 1.0-f > THRESH: 
            edges.append([i,j,1.0-f,'01']); 
print(len(edges))

edges_all_sorted=sorted(edges,key=lambda item: item[0])
fid = 1;
for edge in edges_all_sorted:#mst_edges:  edges_all_sorted
    blocks =2;

    out1=str(varindex[edge[0]])+'\t'+str(edge[3][0])+'\t'+str(varindex[edge[1]])+'\t'+str(edge[3][1]+'\n')
    file_report.write(out1)
    fid +=1;

