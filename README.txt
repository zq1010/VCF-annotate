#This is the README for tempus coding challenge. 
#All relevant files and code can be found at : https://github.com/zq1010/VCF-annotate.git
#Input data provided : test_vcf_data.txt
#Script name: vcf_process.py
#Test dataset used : test_vcf_subset.txt - Smaller sample was used due to run time constraints arising from VEP API.
#Test dataset output : test_vcf_subset_output.txt

#Description of the script

1. Add depth of sequence coverage at variant sites.
2. Add number of reads supporting the alternate variant(s).
3. Add percentage of reads supporting the alternate variant(s).
4. Add minor allele frequency for bi-allelic sites. 
5. Add type of variant - substitution, insertion or deletion.
6. Add impacted gene id from VEP - Only the top hit from VEP API was processed.
7. Add gene consequence from VEP - Only the top hit from VEP API was processed.
8. Add gene impact from VEP - Only the top hit from VEP API was processed.
