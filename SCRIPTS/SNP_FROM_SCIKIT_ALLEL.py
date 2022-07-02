import allel

df = allel.vcf_to_dataframe('SNP/baseplate/VCF_baseplate.vcf', fields='*', alt_number=5)
print(df)

df.to_csv("SNP/baseplate/scikit_allel/VCF_baseplate_5.csv")