import os
import time

# 파일 입출력
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
    num_kmers = len(reference) - k + 1 #k-mer분할
    for i in range(num_kmers - w + 1):
        window_kmers = [reference[j:j + k] for j in range(i, i + w)]
        mn = min(window_kmers)
        mn_pos = i + window_kmers.index(mn)
        index.setdefault(mn, []).append(mn_pos)

    filtered = {m: poses for m, poses in index.items() if len(poses) <= max_occ}
    return filtered

def minimizer_match(
    reference, index, read,
    k=20, w=8, max_mismatch=2, seed_min=2
):
    read_kmers = [read[i:i + k] for i in range(len(read) - k + 1)] # read k-mer 추출
    read_min_pairs = []
    for i in range(len(read_kmers) - w + 1):
        window = read_kmers[i:i + w]
        mn = min(window)
        mn_idx = window.index(mn)
        read_pos = i + mn_idx
        read_min_pairs.append((mn, read_pos))

    delta_counts = {} # delta 카운팅
    for mn, rpos in read_min_pairs:
        if mn not in index:
            continue
        for refpos in index[mn]:
            delta = refpos - rpos
            delta_counts[delta] = delta_counts.get(delta, 0) + 1

    candidates = [d for d, cnt in delta_counts.items() if cnt >= seed_min]
    best_pos, best_mm = -1, max_mismatch + 1

    for delta in candidates: #후보 delta 중 최적 위치 탐색
        if delta < 0 or delta + len(read) > len(reference):
            continue
        window_seq = reference[delta : delta + len(read)]
        mismatches = sum(1 for a, b in zip(read, window_seq) if a != b)
        if mismatches <= max_mismatch and mismatches < best_mm:
            best_mm = mismatches
            best_pos = delta

    return best_pos, best_mm

def reconstruct_genome_with_reads(
    reference, reads, truth_positions, index,
    k=20, w=8, max_mismatch=2, seed_min=2
):
    reconstructed = list(reference)
    matched_reads = 0
    total_reads = len(reads)

    for i, read in enumerate(reads): # 각 read 순회하며 매핑 및 재구성
        true_pos = truth_positions[i]
        pred_pos, mm = minimizer_match(
            reference, index, read,
            k=k, w=w, max_mismatch=max_mismatch, seed_min=seed_min
        )
        if pred_pos == true_pos and mm <= max_mismatch: #매핑 성공 조건
            for j, base in enumerate(read):
                reconstructed[pred_pos + j] = base
            matched_reads += 1

    return ''.join(reconstructed), matched_reads, total_reads

#원본 reference와 비교
def evaluate_reconstruction(reference, reconstructed): 
    assert len(reference) == len(reconstructed)
    matches = sum(1 for a, b in zip(reference, reconstructed) if a == b)
    return matches / len(reference)

def run_mapping_and_evaluation():
    K = 20
    W = 8
    MAX_OCC = 500
    MAX_MM = 2
    SEED_MIN = 2

    pairs = [
        ("1_reference_10M.txt",   "1_1_mammoth_reads_10K.txt",  "1_1_ground_truth_10K.txt"),
        ("1_reference_10M.txt",  "1_2_mammoth_reads_100K.txt",   "1_2_ground_truth_100K.txt"),
        ("1_reference_10M.txt", "1_3_mammoth_reads_1M.txt",  "1_3_ground_truth_1M.txt"),
    ]

    for ref_file, read_file, truth_file in pairs:
        if not (os.path.exists(ref_file) and os.path.exists(read_file) and os.path.exists(truth_file)):
            print(f"> Skipping {ref_file} / {read_file} / {truth_file}: 파일이 존재하지 않음.")
            continue

        print(f"\n=== Processing {ref_file} & {read_file} ===")
        start_time = time.time()

        # 레퍼런스 로딩
        reference = load_reference_to_string(ref_file)
        # Reads & Ground truth 로딩
        reads = load_reads_from_file(read_file)
        truth_positions = load_truth_positions(truth_file)

        # Minimizer 인덱스 생성
        print("> Building minimizer index ...")
        idx_start = time.time()
        index = build_minimizer_index(reference, k=K, w=W, max_occ=MAX_OCC)
        idx_elapsed = time.time() - idx_start
        print(f"  - Unique minimizers after filtering: {len(index)}")
        print(f"  - (Index build time: {idx_elapsed:.2f} sec)")

        # 매핑 및 재구성
        print("> Performing mapping & reconstruction ...")
        recon_start = time.time()
        reconstructed, matched_reads, total_reads = reconstruct_genome_with_reads(
            reference, reads, truth_positions, index,
            k=K, w=W, max_mismatch=MAX_MM, seed_min=SEED_MIN
        )
        recon_elapsed = time.time() - recon_start

        # 정확도 계산
        read_level_acc = matched_reads / total_reads * 100
        base_level_acc = evaluate_reconstruction(reference, reconstructed) * 100

        # 결과 출력
        print(f"  * Read-level mapping accuracy (≦{MAX_MM} mismatch): " # read 예측 위치 == 실제 위치 정확도
              f"{matched_reads}/{total_reads} = {read_level_acc:.2f}%")
        print(f"  * Base-level reconstruction accuracy: {base_level_acc:.2f}%") # 원래 reference와 일치하는 염기의 비율 정확도
        print(f"  * (Mapping & reconstruction time: {recon_elapsed:.2f} sec)")

        total_elapsed = time.time() - start_time
        print(f"  => Total elapsed for this pair: {total_elapsed:.2f} sec "
              f"({total_elapsed/60:.2f} min)\n")

if __name__ == "__main__":
    run_mapping_and_evaluation()
