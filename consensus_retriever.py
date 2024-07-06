import vcf_functions as vf
import sqlite3
import random
from collections import defaultdict as dd
import os
import time
from argparse import ArgumentParser


def main(fasta_file, vcf_file, out_fasta, num_cons, read_len, af_threshold, pos_dict):
    global_now = time.time()
    db_name = f'{fasta_file}.db'
    db = sqlite3.connect(db_name)


    cursor = db.cursor()
    tables_to_drop = ['vcf_genomes', 'variants', 'samples', 'variants_samples']
    for table in tables_to_drop:
        cursor.execute(f"""DROP TABLE IF EXISTS {table}""")

    # Indexing VCF
    print('Indexing VCF')
    index_file = f'{vcf_file}.idx'
    if os.path.isfile(index_file):
        os.remove(index_file)
    samples = vf.index_vcf(vcf_file, index_file, db_name) # saving samples from vcf
    print('VCF index finished', time.time()-global_now)

    # Parsing genome fasta
    gen_tab = False
    cursor.execute("""SELECT name FROM sqlite_master WHERE type='table'""")
    data = cursor.fetchall()
    for entry in data:
        if 'genome' in entry:
            gen_tab = True
    if not gen_tab:
        print('Parsing genome fasta')
        chr_list = vf.parse_fasta(fasta_file, db_name)
    else:
        print('Genome table exists')

    print('Retrieving variants')

    if pos_dict == False: # if chromosomes and positions were not provided
        # Choosing random chromosomes
        cursor.execute(f"""SELECT g.id, g.chr, g.len FROM (
                    SELECT id, chr FROM vcf_genomes ORDER BY RANDOM() LIMIT {num_cons}
                    ) AS v 
                    JOIN genome g ON v.chr = g.chr """)

        chr_positions = cursor.fetchall()

        # Choosing random substrings from chromosomes
        pos_dict = vf.select_regions(chr_positions, read_len)

    

    all_vcf_entries = []
    for chr, positions in pos_dict.items():
        if chr not in chr_list:
            continue
        positions = list(set(positions))
        pos_dict[chr] = positions
        for position in positions:
            vcf_entries = vf.query_vcf(vcf_file, index_file, chr, position, position+read_len)
            all_vcf_entries.extend(vcf_entries)
    print(time.time()-global_now)

    # Creating database for variants
    print('Inserting variants in db')
    cursor.execute("""CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY,
            var_id INTEGER NOT NULL,
            chr TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            alt_order TEXT NOT NULL,
            gt_idx INTEGER NOT NULL,
            af REAL NOT NULL
            
    )""")

    # Creating database for samples
    cursor.execute("""CREATE TABLE IF NOT EXISTS samples (
            id INTEGER PRIMARY KEY,
            sample TEXT NOT NULL,
            gt TEXT NOT NULL
            
    )""")

    # Creating database to link variants and samples
    cursor.execute("""CREATE TABLE IF NOT EXISTS variants_samples (
            variant_id INTEGER NOT NULL,
            sample_id INTEGER NOT NULL,
            FOREIGN KEY(variant_id) REFERENCES variants(id),
            FOREIGN KEY(sample_id) REFERENCES samples(id),
            PRIMARY KEY(variant_id, sample_id)
            )""")

    # Inserting entries into databases
    vf.process_variants(samples, all_vcf_entries, af_threshold, db_name)
    print(time.time()-global_now)
    print('Generating consensuses')
    out_list = []
    for chr, positions in pos_dict.items():

        cursor.execute(f"""
            SELECT id, seq FROM genome WHERE chr = '{chr}'
            """)
        sequence = cursor.fetchone()[1]
        len_seq = len(sequence)
        substrings = []
        for ip, position in enumerate(positions):
            if position+read_len-1 > len_seq:
                position = len_seq-read_len+1
                positions[ip] = position
        substrings = [sequence[position-1:position+read_len-1] for position in positions if position+read_len-1 ]

        now = time.time()
        for num_string, position in enumerate(positions):

            cursor.execute(f"""SELECT v.id, v.var_id, v.pos, v.alt_order, s.id, s.sample, s.gt  FROM 
                (SELECT * FROM variants WHERE chr = '{chr}' AND (pos BETWEEN {position} AND {position+read_len-1})) as v 
                JOIN variants_samples vs ON v.id = vs.variant_id
                JOIN samples s ON vs.sample_id = s.id
                """)
            data = cursor.fetchall()
            if len(data)==0:
                continue
            prev_id = data[0][0]
            prev_var_id = data[0][1]
            pos_list = []
            alts_list = []
            samples_gts = dd(list)
            samples_gts[data[0][5]].append(data[0][6])
            pos_list.append(data[0][2])
            alts_list.append(data[0][3])
            for entry in data:
                if entry[0]!=prev_id and entry[2] in pos_list:
                    continue
                elif entry[0] != prev_id:
                    pos_list.append(entry[2])
                    alts_list.append(entry[3])
                    prev_id = entry[0]
                    prev_var_id = entry[1]
                samples_gts[entry[5]].append(entry[6])


            # Generating consensus alleles for samples grouped by genotypes
            substring = substrings[num_string]
            sorted_indices = sorted(range(len(pos_list)), key=lambda k: pos_list[k])

            sorted_pos_list = sorted(pos_list)
            rel_pos_list = list(map(lambda x: x-position, sorted_pos_list))

            sorted_alts_list = [alts_list[i] for i in sorted_indices]
            sample_id = -1
            samples_of_genotypes_dict = dd(list)
            for samp, genotypes in samples_gts.items():
                sample_id+=1
                sorted_genotypes = [genotypes[i] for i in sorted_indices]
                gts_name = ','.join(sorted_genotypes)
                samples_of_genotypes_dict[gts_name].append(samp)
            for gts, samples in samples_of_genotypes_dict.items():
                gts = gts.split(',')
                read1 = []
                read2 = []
                prev_pos = 0
                reads_dict = dd(list)

                for i, rel_pos in enumerate(rel_pos_list):
                    gt = gts[i].replace('/', '|').split('|')
                    local_alts = sorted_alts_list[i].split(',')
                    reads_dict[f'read1'].append(substring[prev_pos:rel_pos])
                    reads_dict[f'read2'].append(substring[prev_pos:rel_pos])
                    ref_snp = local_alts[0]
                    prev_pos = rel_pos+len(ref_snp)
                    for allele in range(2):
                        snp_id = (gt[allele])
                        if snp_id != '.':
                            try:
                                snp = local_alts[int(snp_id)]
                            except:
                                snp = ref_snp
                        else:
                            snp = ref_snp

                        reads_dict[f'read{allele+1}'].append(snp)
                if prev_pos<len(substring):
                    reads_dict['read1'].append(substring[prev_pos:])
                    reads_dict['read2'].append(substring[prev_pos:])

                for allele in range(2):
                    header = f">chr={chr};pos={position-1}-{position+read_len-1};samples={','.join(samples)};allele={allele+1}"
                    new_seq = ''.join(reads_dict[f'read{allele+1}'])
                    out_list.append(header)
                    out_list.append(new_seq)

    fasta = open(f'{out_fasta}', 'w')
    fasta.write('\n'.join(out_list))
    fasta.close()
    print('Total time', time.time()-global_now)
    
    db.commit()
    db.close()

if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument('-n','--num_cons', help='Amount of consensus sequences to retrieve', type=int, required=True)
    parser.add_argument('-l','--length', help='Length of consensus', type=int, required=True) 
    parser.add_argument('-f','--frequency', metavar="[0-100]", help='Lower threshold for allele frequency (percents)',
                        type=int, choices=range(0,101), required=False, default=0)
    parser.add_argument('-o','--out', help='Name of output fasta', type=str, required=True)
    parser.add_argument('-g','--genome', help='Input genome fasta', type=str, required=True)
    parser.add_argument('--vcf', help='Input VCF file', type=str, required=True)
    parser.add_argument('--chr_pos', help="Provide desired chromosomes and positions as 'chr1:123,456;chr2:789,890' (without quotes) ", type=str, required=False, default=False)

    args = parser.parse_args()

    chr_pos = args.chr_pos
    if chr_pos:
        try:
            pos_dict = {}
            per_chr = chr_pos.split(';')
            for pos_list in per_chr:
                chr_name = pos_list.split(':')[0]
                positions = pos_list.split(':')[1].split(',')
                positions = list(map(int, positions))
                pos_dict[chr_name] = positions
        except:
            print('Wrong format of --chr_pos')
            raise Exception
    else:
        pos_dict = False



    main(args.genome, args.vcf, args.out, args.num_cons, args.length, args.frequency, pos_dict)
