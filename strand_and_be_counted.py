#!/usr/bin/env python
import argparse
import subprocess
from multiprocessing import Pool
import os
import pandas as pd

def run_command(cmd):
    print (" $",cmd)
    subprocess.run(cmd, shell=True, check=True)

def select_chrs(df, threshold):
    # Sort dataframe based on length in descending order
    sorted_df = df.sort_values(by='len', ascending=False)    

    total_len = 0
    selected_chrs = []
    # If the largest 'chr' exceeds the threshold, and no other 'chr' is below the threshold,
    # then select the smallest 'chr'
    if sorted_df.iloc[0]['len'] > threshold and all(x > threshold for x in sorted_df['len']):
        return [sorted_df.iloc[-1]['chr']]    
    # Iterate over the sorted dataframe
    for _, row in sorted_df.iterrows():
        if total_len + row['len'] <= threshold:
            total_len += row['len']
            selected_chrs.append(row['chr'])
        if total_len >= threshold:
            break
    return ' '.join(selected_chrs)

def get_concordant_count(fn):
    lines = open(fn).readlines()
    tot = 0
    conc = 0
    for line in lines:
        if 'aligned concordantly exactly 1 time' in line or 'aligned concordantly >1 times' in line:
            conc+=float(line.strip().split()[1].replace("(","").replace("%)",""))
        elif 'overall alignment rate' in line:
            tot+=float(line.strip().split()[0].replace("%",""))
    return (tot, conc, conc/tot) 

def parse_htseq_output(file_path):
    # Initialize the output dictionary with default values
    output = {
        '__no_feature': 0,
        '__ambiguous': 0,
        '__too_low_aQual': 0,
        '__not_aligned': 0,
        '__alignment_not_unique': 0,
        'assigned_reads': 0
    }
    
    with open(file_path, 'r') as f:
        for line in f:
            # Split each line into its two components: name and count
            name, count = line.strip().split()
            count = int(count)
            
            # Check if the line starts with '__'
            if name.startswith('__'):
                output[name] = count
            else:
                # Add to the total_mapped_reads if it's not a special line
                output['assigned_reads'] += count
    
    return output

message = "*** strand_and_be_counted ***\nA tool for determining the correct --fr/--rf/--ff options in HISAT2 and --strand options in htseq-count for RNA-seq projects.\nRequires the following installed, preferably by conda/mamba:\nmamba create -n strand_and_be_counted -c bioconda htseq hisat2 pandas samtools\n"

parser = argparse.ArgumentParser(description = message, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-r", "--ref", help="reference FASTA")
parser.add_argument("-g", "--gtf", help="reference GTF")
parser.add_argument("--r1", help="R1.fastq.gz")
parser.add_argument("--r2", help="R2.fastq.gz")
parser.add_argument("--num_reads", help="Number of reads to check. -1 = all [%(default)s]", default=-1, type=int)
parser.add_argument("--genome_prop", help="Proportion of genome to use [%(default)s]", default=1.0, type=float)
parser.add_argument("--idattr", help="htseq count gene name attribute [%(default)s]", default='gene_id')
parser.add_argument("--out", help="Output directory [%(default)s]", default='./out')
parser.add_argument("--threads", help="Threads [%(default)s]", default=64)
opts = parser.parse_args()

# 1 | Create output dir
if not os.path.exists(opts.out):
    os.mkdir(opts.out)

# 2 | Subsample reads
r1_fn = "{}/R1.fastq.gz".format(opts.out)
r2_fn = "{}/R2.fastq.gz".format(opts.out)
if opts.num_reads != -1:
    print ("* Subsampling reads")
    r1_fn = "{}/R1.fastq.gz".format(opts.out)
    r2_fn = "{}/R2.fastq.gz".format(opts.out)
    num_lines = opts.num_reads * 4
    subsample_commands = [
        'zcat {} | head -n {} | gzip > {}'.format(opts.r1, num_lines, r1_fn),
        'zcat {} | head -n {} | gzip > {}'.format(opts.r2, num_lines, r2_fn)
    ]
    with Pool(2) as p:
        p.map(run_command, subsample_commands)


else:
    os.symlink(os.path.abspath(opts.r1), os.path.abspath(r1_fn))
    os.symlink(os.path.abspath(opts.r2), os.path.abspath(r2_fn))


# 3 | Subsample genome if required
ref_fn="{}/ref.fasta".format(opts.out)
gtf_fn="{}/ref.gtf".format(opts.out)
if opts.genome_prop < 1.0:
    print ("* Subsampling genome")

    run_command('samtools faidx {} --fai-idx {}/ref.idx'.format(opts.ref, opts.out))
    # Get longest
    df = pd.read_csv('{}/ref.idx'.format(opts.out), sep='\t', names = ['chr','len','cumlen','d','d2'])
    gensize = int(df['cumlen'].max())
    threshold = gensize * opts.genome_prop
    print ("   - Selecting approx {} bp from {} bp genome".format(threshold, gensize))
    selected_chrs = select_chrs(df, threshold)
    print ("   - {}".format(selected_chrs))
    #max_chr = df.loc[df['len'].idxmax(),'chr']
    # Get ref fasta sequence
    run_command('samtools faidx {} {} > {}'.format(opts.ref, selected_chrs, ref_fn))
    run_command("cat {} | egrep -w '^{}' > {}".format(opts.gtf, selected_chrs.replace(" ","|"), gtf_fn))

else:

    os.symlink(os.path.abspath(opts.ref), os.path.abspath(ref_fn))
    os.symlink(os.path.abspath(opts.gtf), os.path.abspath(gtf_fn))
    pass
# 4 | index genome
print ("* Index genome")
run_command('hisat2-build {} {} -p {}'.format(ref_fn, ref_fn, opts.threads))

# 5 | Align 
print ("* Align reads")
rf_bam="{}/rf.bam".format(opts.out)
fr_bam="{}/fr.bam".format(opts.out)
ff_bam="{}/ff.bam".format(opts.out)
rf_stats="{}/rf.stats".format(opts.out)
fr_stats="{}/fr.stats".format(opts.out)
ff_stats="{}/ff.stats".format(opts.out)

aln_threads = int((opts.threads / 3) + 1)
aln_cmds = [
    'hisat2 -p {} -x {} -1 {} -2 {} --fr 2> {} | samtools view -bShF 4 | samtools sort > {}'.format(aln_threads, ref_fn, r1_fn, r2_fn, fr_stats, fr_bam),
    'hisat2 -p {} -x {} -1 {} -2 {} --rf 2> {} | samtools view -bShF 4 | samtools sort > {}'.format(aln_threads, ref_fn, r1_fn, r2_fn, rf_stats, rf_bam),
    'hisat2 -p {} -x {} -1 {} -2 {} --ff 2> {} | samtools view -bShF 4 | samtools sort > {}'.format(aln_threads, ref_fn, r1_fn, r2_fn, ff_stats, ff_bam)
]

with Pool(3) as p:
    p.map(run_command, aln_cmds)

# 6 | Get best hits
print ("* Get best configuration")
fr_metrics = get_concordant_count(fr_stats)
rf_metrics = get_concordant_count(rf_stats)
ff_metrics = get_concordant_count(ff_stats)
data = [
    fr_metrics,
    rf_metrics,
    ff_metrics
]
mdf = pd.DataFrame(data, columns=['total_map_rate', 'concordant_map_rate', 'proportion_concordant'], index=['fr','rf','ff'])
print (">>> *** HISAT2 results")
print ('\n'.join(['>>> ' + row for row in mdf.to_string().split('\n')]))
mdf['bam']=[fr_bam, rf_bam, ff_bam]
best_ori = mdf['proportion_concordant'].idxmax()

print (">>> Recommend using --{} for HISAT2 (Check table above to compare proportion_concordant values)".format(best_ori))
 
# 7 | Run htseq-count
print ("* Get counts")
bam = mdf.loc[best_ori, 'bam']
srt_bam = bam.replace(".bam", ".srt.bam")
run_command('samtools sort {} > {}'.format(bam, srt_bam))
run_command('samtools index {}'.format(srt_bam))
f_htseq = "{}/forward.htseq".format(opts.out)
r_htseq = "{}/reverse.htseq".format(opts.out)
u_htseq = "{}/unstrand.htseq".format(opts.out)
hts_cmds = [
    'htseq-count -t exon -r pos -s yes -i {} {} {} > {}'.format(opts.idattr, srt_bam, gtf_fn, f_htseq),
    'htseq-count -t exon -r pos -s reverse -i {} {} {} > {}'.format(opts.idattr, srt_bam, gtf_fn, r_htseq),
    'htseq-count -t exon -r pos -s no -i {} {} {} > {}'.format(opts.idattr, srt_bam, gtf_fn, u_htseq)
]
with Pool(3) as p:
    p.map(run_command, hts_cmds)

# 8 | Parse htseq-count output
print ("* Parse counts")
count_data = [
    parse_htseq_output(f_htseq),
    parse_htseq_output(r_htseq),
    parse_htseq_output(u_htseq)
]

cdf = pd.DataFrame(count_data)
cdf.index = ['yes','reverse','no']
print (">>>")
print (">>> *** htseq-count results using --{} bam data".format(best_ori))
print ('\n'.join(['>>> ' + row for row in cdf.to_string().split('\n')]))

# If non stranded, would expect the assignment rate to be twice either of the other two rates (+ or -)
# So we will reduce it be 33% (not 50%) as this is too strict) before determining which is best
mod_cdf = cdf.copy(deep=True)
mod_cdf.loc['no','assigned_reads'] = int(mod_cdf.loc['no', 'assigned_reads'] * 0.66)
best_strand_indices = mod_cdf[mod_cdf['assigned_reads'] == mod_cdf['assigned_reads'].max()].index.tolist()
if len(best_strand_indices) > 1 and 'no' in best_strand_indices:
    best_strand_indices.remove('no')

if len(best_strand_indices) >1:
    print (">>> Multiple best strands were identified {}. Stopping".format(best_strand_indices))
    sys.exit()
else:
    print (">>> Recommend using -s {} for htseq-count (Check table above to compare assigned_read values)".format(best_strand_indices[0]))

# 9 | Export results
print ("* Export results")
with open("{}/results.txt".format(opts.out), 'w') as o:
    o.write("*** HISAT2 results\n")
    o.write(mdf.drop(columns=['bam']).to_string())
    o.write("\n*** htseq-count results using --{} bam data\n".format(best_ori))
    o.write(cdf.to_string())
    o.write("\n\n")
    o.write("Recommend using --{} for HISAT2\n".format(best_ori))
    o.write("Recommend using -s {} for htseq-count\n".format(best_strand_indices[0]))
print ("* All done! For results and recommendations, check {}/results.txt".format(opts.out))
