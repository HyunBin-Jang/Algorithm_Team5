import os


def load_reference_to_string(filename):
    with open(filename, 'r') as f:
        return f.read().strip()

def load_reads_from_file(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f]

def load_truth_positions(filename):
    return [int(line.strip()) for line in open(filename, 'r')]


def build_minimizer_index(reference, k=20, w=8, max_occ=500):
    index = {}
    num_kmers = len(reference) - k + 1
    for i in range(num_kmers - w + 1):
        window_kmers = [reference[j:j + k] for j in range(i, i + w)]
        mn = min(window_kmers)
        mn_pos = i + window_kmers.index(mn)
        index.setdefault(mn, []).append(mn_pos)

    filtered = {m: poses for m, poses in index.items() if len(poses) <= max_occ}
    return filtered


def minimizer_match(
    reference,index,read,
    k=20,w=8,
    max_mismatch=2,seed_min=2
):
    read_kmers = [read[i:i + k] for i in range(len(read) - k + 1)]
    read_min_pairs = []
    for i in range(len(read_kmers) - w + 1):
        window = read_kmers[i:i + w]
        mn = min(window)
        mn_idx = window.index(mn)
        read_pos = i + mn_idx
        read_min_pairs.append((mn, read_pos))

    delta_counts = {}
    for mn, rpos in read_min_pairs:
        if mn not in index:
            continue
        for refpos in index[mn]:
            delta = refpos - rpos
            delta_counts[delta] = delta_counts.get(delta, 0) + 1

    candidates = [d for d, cnt in delta_counts.items() if cnt >= seed_min]
    best_pos, best_mm = -1, max_mismatch + 1

    for delta in candidates:
        if delta < 0 or delta + len(read) > len(reference):
            continue
        window_seq = reference[delta : delta + len(read)]
        mismatches = sum(1 for a, b in zip(read, window_seq) if a != b)
        if mismatches <= max_mismatch and mismatches < best_mm:
            best_mm = mismatches
            best_pos = delta

    return best_pos, best_mm


def reconstruct_genome_with_reads(
    reference,
    reads,
    truth_positions,
    index,
    k=20, w=8, max_mismatch=2, seed_min=2
):
    reconstructed = list(reference)
    matched_reads = 0
    total_reads = len(reads)

    for i, read in enumerate(reads):
        true_pos = truth_positions[i]
        pred_pos, mm = minimizer_match(
            reference, index, read,
            k=k, w=w, max_mismatch=max_mismatch, seed_min=seed_min
        )
        # mismatch 허용 범위에서 복원
        if pred_pos == true_pos and mm <= max_mismatch:
            for j, base in enumerate(read):
                reconstructed[pred_pos + j] = base
            matched_reads += 1

    return ''.join(reconstructed), matched_reads, total_reads


def evaluate_reconstruction(reference, reconstructed):
    assert len(reference) == len(reconstructed)
    L = len(reference)
    matches = sum(1 for a, b in zip(reference, reconstructed) if a == b)
    return matches / L


def run_mapping_and_evaluation():

    K = 20
    W = 8
    MAX_OCC = 500
    MAX_MM = 2
    SEED_MIN = 2

    pairs = [
        ("reference_1M.txt",  "mammoth_reads_100K.txt",  "ground_truth_100K.txt"),
        ("reference_10M.txt", "mammoth_reads_1M.txt",   "ground_truth_1M.txt"),
        ("reference_100M.txt","mammoth_reads_10M.txt",  "ground_truth_10M.txt"),
        #("reference_1T.txt",  "mammoth_reads_100M.txt", "ground_truth_100M.txt"),
    ]

    for ref_file, read_file, truth_file in pairs:
        if not (os.path.exists(ref_file) and os.path.exists(read_file) and os.path.exists(truth_file)):
            print(f"> Skipping {ref_file} / {read_file} / {truth_file}: 파일이 존재하지 않음.")
            continue

        print(f"\n=== Processing {ref_file} & {read_file} ===")
        reference = load_reference_to_string(ref_file)
        reads = load_reads_from_file(read_file)
        truth_positions = load_truth_positions(truth_file)

        print("> Building minimizer index ...")
        index = build_minimizer_index(reference, k=K, w=W, max_occ=MAX_OCC)
        print(f"  - Unique minimizers after filtering: {len(index)}")

        print("> Performing mapping & reconstruction ...")
        reconstructed, matched_reads, total_reads = reconstruct_genome_with_reads(
            reference, reads, truth_positions, index,
            k=K, w=W, max_mismatch=MAX_MM, seed_min=SEED_MIN
        )

        read_level_acc = matched_reads / total_reads * 100
        base_level_acc = evaluate_reconstruction(reference, reconstructed) * 100

        print(f"  * Read-level mapping accuracy (exact 위치 매칭 & ≤{MAX_MM} mismatch): "
              f"{matched_reads}/{total_reads} = {read_level_acc:.2f}%")
        print(f"  * Base-level reconstruction accuracy: {base_level_acc:.2f}%")

if __name__ == "__main__":
    run_mapping_and_evaluation()
