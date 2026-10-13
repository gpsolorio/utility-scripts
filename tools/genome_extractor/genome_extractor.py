import argparse
from collections import defaultdict

# Read in a BED file and extract chromosomal regions. Returns output in the following format:
# chr1 → [(start1, end1), (start2, end2), ...]
# chr2 → [(start3, end3), ...]

parser = argparse.ArgumentParser(description='Extract genome regions from FASTA using BED')
parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
parser.add_argument('-b', '--bed', required=True, help='Input BED file')
args = parser.parse_args()

coord_dict = defaultdict(list)

# Read the whole bed file
with open(args.bed) as bf:
    for line in bf:
        if not line.strip():
            continue  # skip empty lines
        bf_split = line.split()
        chrom, start, end = bf_split[0:3]
        coord_dict[chrom].append((int(start), int(end)))


# Read one chromosome at a time and then cross reference it with the previous bed file

# can be neater by using a parser that already exists

with open(args.fasta) as ff:
    seq_lst = []
    extracted_seqs = []
    chr_extractions = {}
    uid = None

    for line in ff:
        if line.startswith('>'):
            if uid is not None:
                seq = ''.join(seq_lst)

                for start, end in coord_dict[uid]:
                    extracted_seqs.append(seq[start:end])
                chr_extractions[uid] = extracted_seqs
                seq_lst = []
                extracted_seqs = []
            
            uid = line[1:].split()[0]
            # coords_lst = coord_dict[uid]
            continue
        else:
            seq_lst.append(line.strip())
    
    seq = ''.join(seq_lst)
    for start, end in coord_dict[uid]:
        extracted_seqs.append(seq[start:end])
    chr_extractions[uid] = extracted_seqs

print(chr_extractions)