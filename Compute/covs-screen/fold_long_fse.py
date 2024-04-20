import os
import sys

seq_dir = "long_fse_seqs"
ct_dir = "long_fse_structs"
file_list = "long_fse_seqs_files.txt"
batch_script = "fold_long_fse_seqs.sh"


def write_folding_script():
    if not os.path.isdir(ct_dir):
        os.mkdir(ct_dir)
    files = [f for f in os.listdir(seq_dir) if f.endswith(".fasta")]
    n_files = len(files)
    with open(file_list, "w") as f:
        f.write("\n".join(files))
    with open(batch_script, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-00:15
#SBATCH -p short
#SBATCH --array=1-{n_files}
#SBATCH -o job_%A_%a.out
#SBATCH -e job_%A_%a.err

python fold_long_fse.py $SLURM_ARRAY_TASK_ID""")
    cmd = f"sbatch {batch_script}"
    os.system(cmd)


def fold_rna(fasta_file, max_structs=100):
    ct_file = os.path.join(ct_dir, fasta_file.replace(".fasta", ".ct"))
    fasta_file = os.path.join(seq_dir, fasta_file)
    fold_cmd = f"Fold -M {max_structs} {fasta_file} {ct_file}"
    os.system(fold_cmd)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        write_folding_script()
    elif len(sys.argv) == 2:
        file_number = int(sys.argv[1])
        with open(file_list) as f:
            files = [l.strip() for l in f]
        fasta_file = files[file_number - 1]
        fold_rna(fasta_file)
    else:
        raise ValueError("invalid number of arguments")

