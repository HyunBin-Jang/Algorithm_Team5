import os
import random

def generate_reference_to_file(filename, length, chunk_size=10_000_000):
    bases = ['A', 'C', 'G', 'T']
    with open(filename, 'w') as f:
        written = 0
        while written < length:
            to_write = min(chunk_size, length - written)
            seq_chunk = ''.join(random.choices(bases, k=to_write))
            f.write(seq_chunk)
            written += to_write

def simulate_ancient_reads_to_file(
    reference_file,
    reads_filename,
    truth_filename,
    read_len,
    num_reads,
    mutation_rate=0.02
):
    with open(reference_file, 'r') as f:
        reference_seq = f.read().strip()
    L = len(reference_seq)

    with open(reads_filename, 'w') as rf, open(truth_filename, 'w') as tf:
        for i in range(num_reads):
            start = random.randint(0, L - read_len)
            original = list(reference_seq[start : start + read_len])
            mutated = []
            for j, base in enumerate(original):
                prob = mutation_rate
                if j < 5 or j >= read_len - 5:
                    prob *= 5
                if random.random() < prob:
                    if base == 'C':
                        mutated.append('T')
                    elif base == 'G':
                        mutated.append('A')
                    else:
                        mutated.append(random.choice(['A', 'C', 'G', 'T']))
                else:
                    mutated.append(base)
            read_str = ''.join(mutated)
            rf.write(read_str + '\n')
            tf.write(f"{start}\n")

def run_simulation():
    settings = [
        (1_000_000, 10_000,"1M",   "10K"),
        (10_000_000, 100_000,"10M",  "100K"),
        (100_000_000, 1_000_000,"100M", "1M"),
        (1_000_000_000,10_000_000,"1G",  "10M"),
    ]

    # Reference 생성
    for ref_len, _, ref_tag, _ in settings:
        ref_filename = f"reference_{ref_tag}.txt"
        if not os.path.exists(ref_filename):
            generate_reference_to_file(ref_filename, ref_len)

    for ref_len, num_reads, ref_tag, read_tag in settings:
        ref_filename   = f"reference_{ref_tag}.txt"
        reads_filename = f"mammoth_reads_{read_tag}.txt"
        truth_filename = f"ground_truth_{read_tag}.txt"
        read_len = 100

        if not (os.path.exists(reads_filename) and os.path.exists(truth_filename)):
            simulate_ancient_reads_to_file(
                reference_file=ref_filename,
                reads_filename=reads_filename,
                truth_filename=truth_filename,
                read_len=read_len,
                num_reads=num_reads,
                mutation_rate=0.02
            )

if __name__ == "__main__":
    run_simulation()
