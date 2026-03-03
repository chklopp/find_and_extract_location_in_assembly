import argparse
import os
import subprocess
import gzip
import shutil
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path

def move_files_by_pattern(source_dir, target_dir, pattern):
    """
    Moves files matching a pattern from source_dir to target_dir.

    :param source_dir: Path to the folder to search.
    :param target_dir: Path to the destination folder.
    :param pattern: Glob pattern (e.g., "*.pdf", "report_*.csv").
    """
    src = Path(source_dir)
    dest = Path(target_dir)

    # Ensure the target directory exists
    dest.mkdir(parents=True, exist_ok=True)

    # Counter for feedback
    moved_count = 0

    # Iterate through files matching the pattern
    for file_path in src.glob(pattern):
        if file_path.is_file():
            # Define the new path
            new_location = dest / file_path.name

            # Move the file
            shutil.move(str(file_path), str(new_location))
            moved_count += 1
            print(f"Moved: {file_path.name}")

    print(f"--- Completed. Moved {moved_count} files. ---")

def check_softwares(softwares):
    for software in softwares:
        #print(software)
        try:
            # Run the command with no arguments
            result = subprocess.run([software], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            # print(result.returncode)
            # Check the return code
            if result.returncode == 0 or result.returncode == 1 :
                print("The "+software+" command is accessible.")
            else:
                print("The "+software+" command is not accessible.")
                sys.exit(1)
                return False
        except FileNotFoundError:
            print("The "+software+" command is not found in the environment.")
            return False
        except Exception as e:
            print(f"An error occurred: {e}")
            return False

def get_fasta_reader(file_path):
    """Handles both compressed and uncompressed fasta files."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')

def run_vdj_insights(name, species, directory):
    """Runs vdj-insights on a given directory"""
    igtr = name[:2]
    cmd = ["vdj-insights", "annotation", "--receptor-type", "\""+igtr+"\"","--input", directory, "--species", "\""+species+"\"", "--output", directory+"/annotations"]
    with open(name+".log", "w") as f:
        #subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.DEVNULL)
        print(cmd)

def run_minimap2(ref_path, query_path, output_paf, threads, paf_dir):
    """Runs minimap2 and generates a PAF file."""
    cmd = ["minimap2", "-t", str(threads), "-x", "asm5", ref_path, query_path]
    with open(paf_dir+"/"+output_paf, "w") as f:
        subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.DEVNULL)

def process_paf(paf_file, max_dist, paf_dir):
    """
    Parses PAF and merges alignments on the same chromosome 
    if they are within max_dist. Returns the 'best' merged block.
    """
    alignments = defaultdict(list)
    
    with open(paf_dir+"/"+paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if not fields: continue
            # PAF columns: 0:query, 5:target, 7:t_start, 8:t_end, 9:matches, 10:blocklen, 4:strand
            q_id, t_id, strand = fields[0], fields[5], fields[4]
            t_start, t_end = int(fields[7]), int(fields[8])
            score = int(fields[9])
            alignments[(q_id, t_id, strand)].append([t_start, t_end, score])

    best_results = {}

    for (q_id, t_id, strand), locs in alignments.items():
        # Sort by start position
        locs.sort()
        
        merged = []
        if not locs: continue
        
        curr_start, curr_end, curr_score = locs[0]

        for i in range(1, len(locs)):
            next_start, next_end, next_score = locs[i]
            if next_start - curr_end <= max_dist:
                curr_end = max(curr_end, next_end)
                curr_score += next_score
            else:
                merged.append((curr_start, curr_end, curr_score))
                curr_start, curr_end, curr_score = next_start, next_end, next_score
        merged.append((curr_start, curr_end, curr_score))

        # Select the merged block with the highest match score
        best_block = max(merged, key=lambda x: x[2])
        
        if q_id not in best_results or best_block[2] > best_results[q_id]['score']:
            best_results[q_id] = {
                'target': t_id,
                'start': best_block[0],
                'end': best_block[1],
                'strand': strand,
                'score': best_block[2]
            }
            
    return best_results

def extract_sequence(genome_path, chrom, start, end, strand, out_path, asm_name, ref_name, flanquing):
    """Extracts sub-sequence and writes to file."""
    with get_fasta_reader(genome_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == chrom:
                if (start-flanquing) <= 0 or (end+flanquing) >= len(record.seq) :
                    flanquing = 0
                sub_seq = record.seq[start-flanquing:end+flanquing]
                if strand == '-':
                    sub_seq = sub_seq.reverse_complement()
                
                record.seq = sub_seq
                record.id = f"{asm_name}_{ref_name}_{record.id}_{start-flanquing}_{end+flanquing}"
                SeqIO.write(record, out_path, "fasta")
                break

def ensure_dir(directory_path):
    """
    Checks if a directory exists and creates it (including parents) if it doesn't.
    """
    Path(directory_path).mkdir(parents=True, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description="Align references to assemblies and extract best hits.")
    parser.add_argument("-a", "--assemblies", required=True, help="Directory containing genome assemblies")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for minimap2")
    parser.add_argument("-r", "--references", required=True, help="Directory containing reference sequences")
    parser.add_argument("-p", "--paf", required=False, default="PAF", help="Directory in which paf files will be written")
    parser.add_argument("-d", "--dist", type=int, default=10000, help="Max distance to merge alignments default 10Kb")
    parser.add_argument("-o", "--output", default="summary_results.tsv", help="Output summary file")
    parser.add_argument("-x", "--removepaf", default=True, help="Removes paf intermediate paf files")
    #parser.add_argument("-v", "--vdjinsights", default=True, help="Run vdj-insights on each reference")
    parser.add_argument("-f", "--flanquing", type=int, default=50000, help="Flanking region size to add")
    #parser.add_argument("-s", "--vdjspecies", required=False, help="Species used to run vdj-insights, use double quotes")
    
    args = parser.parse_args()

    assembly_files = [os.path.join(args.assemblies, f) for f in os.listdir(args.assemblies) 
                      if f.endswith(('.fasta', '.fa', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'))]
    ref_files = [os.path.join(args.references, f) for f in os.listdir(args.references) 
                 if f.endswith(('.fasta', '.fa', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'))]

    ensure_dir(args.paf)

    softwares = ["minimap2 -h"]
    check_softwares(softwares)
    #if args.vdjinsights == True :
    #    softwares = ["vdj-insights -h"]
    #    check_softwares(softwares)

    refs = []
    # running minimap2, filtering, fasta extraction
    with open(args.output, "w") as summary:
        summary.write("Assembly\tReference\tChr\tStart\tEnd\tStrand\n")

        for asm in assembly_files:
            asm_name = os.path.basename(asm).split('.')[0]+"."+os.path.basename(asm).split('.')[1]
            for ref in ref_files:
                ref_name = os.path.basename(ref).split('.')[0]
                refs.append(ref_name)
                paf_name = f"{asm_name}_{ref_name}.paf"
                
                print(f"Aligning {ref_name} against {asm_name}...")
                run_minimap2(asm, ref, paf_name, args.threads, args.paf)
                
                best_hits = process_paf(paf_name, args.dist, args.paf)
                
                for ref_id, data in best_hits.items():
                    summary.write(f"{asm_name}\t{ref_id}\t{data['target']}\t{data['start']}\t{data['end']}\t{data['strand']}\n")
                    
                    # Extract FASTA
                    out_fasta = f"{asm_name}_{ref_name}.fasta"
                    extract_sequence(asm, data['target'], data['start'], data['end'], data['strand'], out_fasta, asm_name, ref_name, args.flanquing)
                
                # Cleanup temporary PAF
                if os.path.exists(args.paf+"/"+paf_name) and args.removepaf == True: 
                    os.remove(args.paf+"/"+paf_name)

    # running vdj-insigth and producing graphics
    #if args.vdjinsights == True :
    vdjdir = "vdj-insights"
    ensure_dir(vdjdir)
    for r in refs :
        ensure_dir(vdjdir+"/"+r)
        move_files_by_pattern(".",vdjdir+"/"+r,"*_"+r+".fasta")
        #run_vdj_insights(r, args.vdjspecies, vdjdir+"/"+r)

if __name__ == "__main__":
    main()
