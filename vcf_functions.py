import sqlite3
import time
import os
import random
from collections import defaultdict as dd

# For each chosen chromosome, takes a random position on it
def select_regions(chr_positions, read_len):
    pos_dict = dd(list)
    for entry in chr_positions:
        chr_name = entry[1]
        start_pos = random.randint(1, entry[2]-read_len)
        pos_dict[chr_name].append(start_pos)
    for chr, positions in pos_dict.items():
        pos_dict[chr] = sorted(positions)
    return(pos_dict)

# Loads vcf index in memory
def load_index(index_file):
    index = {}
    with open(index_file, 'r') as f:
        for line in f:
            chrom, pos, offset = line.strip().split('\t')
            pos = int(pos)
            offset = int(offset)
            if chrom not in index:
                index[chrom] = []
            index[chrom].append((pos, offset))
    return index


# With random chromosomes and regions, obtains VCF entries with positions
# in given regions
def query_vcf(vcf_file, index_file, chrom, start, end):
    index = load_index(index_file)
    result = []
    with open(vcf_file, 'r') as f:
        for pos, offset in index[chrom]:
            if start <= pos <= end:
                f.seek(offset)
                line = f.readline()
                result.append(line.strip())
    return result

# Creates vcf index file
def index_vcf(vcf_file, index_file, db_name):
    index = {}
    db = sqlite3.connect(db_name)
    cursor = db.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS vcf_genomes (
        id INTEGER PRIMARY KEY,
        chr TEXT NOT NULL
        
        )""")
    with open(vcf_file, 'r') as f:
        offset = 0
        for line in f:
            if line.startswith('##'):
                offset += len(line)
                continue
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                offset += len(line)
                continue
            parts = line.split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            if chrom not in index:
                index[chrom] = []
            index[chrom].append((pos, offset))
            offset += len(line)
    with open(index_file, 'w') as f:
        ind_entries = []
        for chrom in index:
            for pos, offset in index[chrom]:
                ind_entries.append((chrom,))
                f.write(f"{chrom}\t{pos}\t{offset}\n")
    cursor.executemany("INSERT INTO vcf_genomes (chr) VALUES (?)", ind_entries)
    db.commit()
    db.close()
    return(samples)


# Creates sqlite table for genome
def parse_fasta(fasta_file, db_name):
    now = time.time()

    db = sqlite3.connect(db_name)
    cursor = db.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS genome (
            id INTEGER PRIMARY KEY,
            chr TEXT NOT NULL,
            seq TEXT NOT NULL,
            len INTEGER NOT NULL
            
    )""")

    with open(fasta_file, 'r') as fa:
        chr_list = []
        lines = iter(fa)
        chr = next(lines)[1:].strip().split()[0]
        chr_list.append(chr)
        seq = []
        entries = []
        num_entries=0
        for line in lines:
            if line[0] == '>':
                num_entries+=1
                seq = ''.join(seq)
                entries.append((chr, seq, len(seq)))
                if num_entries%500==0:
                    print(num_entries, 'chromosomes')
                    cursor.executemany("INSERT INTO genome (chr, seq, len) VALUES (?,?,?)", entries)
                    entries = []
                chr = line[1:].strip().split()[0]
                chr_list.append(chr)
                seq = []
            else:
                seq.append(line.strip())
        seq = ''.join(seq)
        entries.append((chr, seq, len(seq)))
    fa.close()

    cursor.executemany("INSERT INTO genome (chr, seq, len) VALUES (?,?,?)", entries)
    db.commit()
    db.close()
    print(time.time()-now)
    return chr_list


# Fills out tables for variants and samples
def process_variants(samples, lines, af_threshold, db_name):
    num_entries = 0
    sample_id = 0
    var_id = 0
    db = sqlite3.connect(db_name)
    cursor = db.cursor()
    entries = []
    sample_entries = []
    variant_sample_entries = []
    for line in lines:
        var_id += 1
        values = line.strip().split('\t')
        chr = values[0]
        pos = values[1]
        ref = values[3]
        alt = values[4].split(',')
        alt_order = f'{ref},{values[4].strip()}'
        info = values[7].split(';')
        for i in info:
            if i[:2] == 'AF':
                af = i.split('=')[1].split(',')

        for i in range(len(alt)):
            local_num_entries = num_entries
            if float(af[i]) < af_threshold:
                continue
            num_entries+=1
            entries.append((var_id, chr, pos, ref, alt[i], alt_order, i+1, af[i]))
            local_sample_id = sample_id
            for j in samples:
                local_sample_id+=1
                variant_sample_entries.append((num_entries, local_sample_id))
        if num_entries==local_num_entries:
            continue
        sample_id = local_sample_id
        format =  values[8].split(':')
        gt_id = format.index('GT')
        for i, sample in enumerate(samples):
            sample_gt = values[9+i].split(':')[gt_id]
            sample_entries.append((sample, sample_gt))

        if num_entries%50000 == 0: # Insertions are split into batches to save memory
            print('variants processed', num_entries)
            cursor.executemany("INSERT INTO variants (var_id, chr, pos, ref, alt, alt_order, gt_idx, af) VALUES (?,?,?,?,?,?,?,?)", entries)
            cursor.executemany("INSERT INTO samples  (sample, gt) VALUES (?,?)", sample_entries)
            cursor.executemany("INSERT INTO variants_samples (variant_id, sample_id) VALUES (?,?)", variant_sample_entries)
            entries = []
            sample_entries = []
            variant_sample_entries = []
    cursor.executemany("INSERT INTO samples  (sample, gt) VALUES (?,?)", sample_entries)        
    cursor.executemany("INSERT INTO variants (var_id, chr, pos, ref, alt, alt_order, gt_idx, af) VALUES (?,?,?,?,?,?,?,?)", entries)
    cursor.executemany("INSERT INTO variants_samples (variant_id, sample_id) VALUES (?,?)", variant_sample_entries)
    db.commit()
    db.close()