from Defs_Auxiliares import write_id_file
import vcf


vcf_reader = vcf.Reader(open('SNP/baseplate/VCF_baseplate.vcf', 'r'))

record = next(vcf_reader)

snp_info = {}

mutation_type = [f"Is SNP,Is Indel,Is Transition (ts),Is Deletion,Is Monomorphic,Variant type,Variant subtype"]


for record in vcf_reader:
    mutation_type.append(f"{record.is_snp},{record.is_indel},{record.is_transition},{record.is_deletion},{record.is_monomorphic},{record.var_type},{record.var_subtype}")

print(len(mutation_type))

write_id_file(file_name="pyvcf_snp_baseplate", info=mutation_type, path="SNP/baseplate/PyVCF")
# print(len(snp_info))