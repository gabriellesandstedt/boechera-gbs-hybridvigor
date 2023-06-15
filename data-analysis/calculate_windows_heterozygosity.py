import vcf

# Create a VCF reader object
vcf_reader = vcf.Reader(open('boech_gbs_allsamples_combined_final.vcf', 'r'))

# Set the window size and starting position
window_size = 100000
current_pos = {chromosome: 1 for chromosome in ["CM031545.1", "CM031546.1", "CM031547.1", "CM031548.1", "CM031549.1", "CM031550.1", "CM031551.1"]}

# Initialize counters for homozygous and heterozygous genotypes for all samples
hom_counts = {chromosome: [0] * 156 for chromosome in ["CM031545.1", "CM031546.1", "CM031547.1", "CM031548.1", "CM031549.1", "CM031550.1", "CM031551.1"]}
het_counts = {chromosome: [0] * 156 for chromosome in ["CM031545.1", "CM031546.1", "CM031547.1", "CM031548.1", "CM031549.1", "CM031550.1", "CM031551.1"]}

# Loop through each record in the VCF file
for record in vcf_reader:
    chromosome = record.CHROM

    # If the current record is beyond the current window, print the counts and reset them
    if record.POS > current_pos[chromosome] + window_size:
        print(f"{chromosome}:{current_pos[chromosome]}-{current_pos[chromosome]+window_size-1}: Homozygous: {hom_counts[chromosome]}, Heterozygous: {het_counts[chromosome]}")
        hom_counts[chromosome] = [0] * 156
        het_counts[chromosome] = [0] * 156
        current_pos[chromosome] += window_size

    # Increment the appropriate genotype counter based on the genotype of each sample
    for i, sample in enumerate(record.samples):
        genotype = sample['GT']
        if genotype == '0/0' or genotype == '1/1':
            hom_counts[chromosome][i] += 1
        elif genotype == '0/1':
            het_counts[chromosome][i] += 1

# Print the counts for the final window
print(f"{chromosome}:{current_pos[chromosome]}-{current_pos[chromosome]+window_size-1}: Homozygous: {hom_counts[chromosome]}, Heterozygous: {het_counts[chromosome]}")
