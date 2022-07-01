import re
import numpy as np
import requests, sys, json
input_file = "/Users/qzhang/Downloads/test_vcf_subset.txt"
output_file = "/Users/qzhang/Downloads/test_vcf_subset_output.txt"
outHandle = open(output_file,"w")
info_flag = 0 #Used for printing new header lines
with open(input_file,"r") as inFile:
    for aline in inFile.readlines():
        aline = aline.rstrip()
        if re.match("^\#.*",aline) is None:
            #print(aline)
            cols = aline.split("\t")
            ref = cols[3]
            alt = cols[4] 
            info = cols[7]
            form = cols[8]
            samp_info = cols[9]
            #Get info
            depth = samp_info.split(":")[form.split(":").index("NR")].split(",")
            var_supp = samp_info.split(":")[form.split(":").index("NV")].split(",")
            var_percent = [round((int(x)/int(depth[var_supp.index(x)]))*100,1) for x in var_supp]
            #Add info
            form = form+":DP:VS:VP"
            samp_info = samp_info+":"+",".join(depth)+":"+",".join(var_supp)+":"+",".join(map(str,var_percent))
            if len(depth) == 1:
                af = int(var_supp[0])/int(depth[0])
                if af > 0.5:
                    maf = 1-af
                else:
                    maf = af
                maf = round(maf,3)
                if len(alt) == len(ref):
                    var_type = "substitution"
                elif len(alt) > len(ref):
                    var_type = "insertion"
                else:
                    var_type = "deletion"
                info = info+";VTYPE="+var_type
                form = form+":MAF"
                samp_info = samp_info+":"+str(maf)
            cols[8] = form
            cols[9] = samp_info
            #Reconstruct line
            raline = "\t".join(cols)
            #print(raline)
            #VEP - getting only the first hit from VEP. 
            bline = " ".join(raline.split("\t")[0:5])
            r = requests.post(server+ext,headers=headers,data='{ "variants": ["'+bline+'"]}')
            if not r.ok:
                r.raise_for_status()
            decoded = r.json()
            if 'transcript_consequences' in decoded[0].keys():
                #print(decoded[0]['transcript_consequences'][0])
                gene_name = decoded[0]['transcript_consequences'][0]['gene_id']
                gene_cons = decoded[0]['transcript_consequences'][0]['consequence_terms'][0]
                gene_impact = decoded[0]['transcript_consequences'][0]['impact']
            #Add VEP to info
                info = info+";AGN="+gene_name+";AGC="+gene_cons+";AGI="+gene_impact
            cols[7] = info
            #reconstruct
            raline = "\t".join(cols)
            outHandle.write(raline+"\n")
        #Adding header lines using flag
        elif re.match("^\##FORMAT.*",aline) is not None:
            if info_flag == 0:
                outHandle.write("##INFO=<ID=VTYPE,Number=1,Type=String,Description=\"Type of variant.Only applicable for bi-allelic sites\">\n")
                outHandle.write("##INFO=<ID=AGN,Number=1,Type=String,Description=\"Affected gene Ensembl ID. Ensembl ID of top vep hit\">\n")
                outHandle.write("##INFO=<ID=AGC,Number=1,Type=String,Description=\"Affected gene consequence. Gene consequence of top vep hit\">\n")
                outHandle.write("##INFO=<ID=AGI,Number=1,Type=String,Description=\"Affected gene Impact. Gene impact of top vep hit\">\n")
                outHandle.write(aline+"\n")
                info_flag = 1
            else:
                outHandle.write(aline+"\n")
        elif re.match("^\#CHROM.*",aline) is not None:
                outHandle.write("##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Depth of sequence coverage at the site\">\n")
                outHandle.write("##FORMAT=<ID=VS,Number=.,Type=Integer,Description=\"Number of reads supporting the alt variant(s)\">\n")
                outHandle.write("##FORMAT=<ID=VP,Number=.,Type=Integer,Description=\"Percentage of reads supporting the alt variant(s)\">\n")
                outHandle.write("##FORMAT=<ID=MAF,Number=.,Type=Integer,Description=\"Minor (second most common allele) allele frequency. Only applicable for bi-allelic sites\">\n")
                outHandle.write(aline+"\n")
        else:
            outHandle.write(aline+"\n")
            #print (maf)
            #print(var_percent)