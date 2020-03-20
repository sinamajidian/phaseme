
from sys import argv


vcf_file_address =  argv[1] # "sample2.vcf"

lines_het=[]
if 1:
#def read_vcf_file_with_truth(vcf_file_address):    
  
    vcf_file = open(vcf_file_address,'r')
  
    for line in vcf_file:

        line_strip = line.strip() 
              
        if line_strip.startswith('#'):
            lines_het.append(line_strip)
        else:
    
            line_parts=line_strip.split('\t') 
            if len(line_strip.split('\t')[3]+line_strip.split('\t')[4]) ==2: # for removing non bi alleilc 
                val=float(line_parts[7].split(';')[8].split('=')[1])
                if (val !=0) : #and val !=1 
                    # 'AC=1;AF=0.000199681;AN=5008;NS=2504;DP=21490;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0;EUR_AF=0;SAS_AF=0;AA=.|||;VT=SNP'
                    line_parts[-1]='0|1' # [5] EAS 

                    lines_het.append('\t'.join(line_parts))


        
    vcf_file_fake_name= vcf_file_address[:-19]+"_fake_n0.vcf"
    vcf_file_fake = open(vcf_file_fake_name,'w');  # phased_vcf_dic



    for header_line in lines_het:
        vcf_file_fake.write(header_line+'\n')
       
    vcf_file_fake.close()
                
             
