import time
import pydivsufsort

def build_suffix_array(text: str) -> list:
    return list(pydivsufsort.divsufsort(text.encode()))

def _compare_prefix(text: str, pos: int, pattern: str) -> int:
    slice_ref = text[pos : pos + len(pattern)]
    if slice_ref == pattern:
        return 0
    return -1 if slice_ref < pattern else 1

def find_exact_matches(sa: list, text: str, pattern: str) -> list:
    n = len(sa)
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        if _compare_prefix(text, sa[mid], pattern) < 0:
            left = mid + 1
        else:
            right = mid
    start = left

    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        if _compare_prefix(text, sa[mid], pattern) <= 0:
            left = mid + 1
        else:
            right = mid
    end = left

    return [sa[i] for i in range(start, end)]


def count_mismatches(s1: str, s2: str) -> int:
    return sum(1 for a, b in zip(s1, s2) if a != b)


def align_single_read(sa: list, reference: str, read: str, k: int) -> int:
    L = len(read)
    t = k + 1
    chunk_size = L // t
    candidates = set()

    for i in range(t):
        start_i = i * chunk_size
        end_i = L if (i == t - 1) else (i + 1) * chunk_size
        chunk = read[start_i:end_i]
        if not chunk:
            continue

        matches = find_exact_matches(sa, reference, chunk)
        for match_pos in matches:
            cand_start = match_pos - start_i
            if cand_start < 0 or cand_start + L > len(reference):
                continue
            candidates.add(cand_start)

    best_pos = -1
    best_mm = k + 1
    for pos in candidates:
        seg = reference[pos : pos + L]
        mm = count_mismatches(seg, read)
        if mm < best_mm:
            best_mm = mm
            best_pos = pos
            if mm == 0:
                break

    return best_pos if best_mm <= k else -1


def align_all_reads(reference: str, reads: list, k: int) -> dict:
    print("🔄 Suffix array 생성 중… (C 라이브러리 divsufsort 사용)")
    t0 = time.time()
    sa = build_suffix_array(reference)
    t1 = time.time()
    print(f"✅ SA 생성 완료 (소요: {t1-t0:.2f}초, length={len(sa)})\n")

    total_reads = len(reads)
    alignments = {}
    PRINT_INTERVAL = max(1, total_reads // 100)

    for idx, read in enumerate(reads):
        pos = align_single_read(sa, reference, read, k)
        alignments[idx] = pos
        if (idx + 1) % PRINT_INTERVAL == 0 or (idx + 1) == total_reads:
            pct = (idx + 1) / total_reads * 100
            print(f"🔄 처리 중: {idx + 1}/{total_reads} reads ({pct:.1f}%) 완료")

    print()
    return alignments


def reconstruct_mammoth(reference: str, reads: list, alignments: dict) -> str:
    ref_list = list(reference)
    if not reads:
        return reference

    L = len(reads[0])
    for idx, read in enumerate(reads):
        pos = alignments.get(idx, -1)
        if pos < 0:
            continue
        for i, base in enumerate(read):
            ref_list[pos + i] = base

    return "".join(ref_list)


def main():
    k = 2  # 허용 mismatch 개수

    start_time = time.time()

    with open("../genome_generation/4_reference_1B.txt", "r") as f_ref:
        reference = f_ref.read().strip()

    with open("../genome_generation/4_mammoth_reads_10M.txt", "r") as f_reads:
        reads = [line.strip() for line in f_reads if line.strip()]

    with open("../genome_generation/4_ground_truth_10M.txt", "r") as f_gt:
        ground_truth = [int(line.strip()) for line in f_gt if line.strip()]

    alignments = align_all_reads(reference, reads, k)

    total_reads = len(reads)
    correct_count = 0
    for idx, true_pos in enumerate(ground_truth):
        pred_pos = alignments.get(idx, -1)
        if pred_pos == true_pos:
            correct_count += 1

    accuracy = correct_count / total_reads * 100
    print(f"\n🔍 Total reads: {total_reads}")
    print(f"🔍 Correctly aligned reads: {correct_count}")
    print(f"🔍 Alignment accuracy: {accuracy:.2f}%\n")

    reconstructed = reconstruct_mammoth(reference, reads, alignments)
    with open("reconstructed_mammoth_dna.txt", "w") as f_out:
        f_out.write(reconstructed)
    print("✅ 복원된 mammoth DNA 시퀀스가 'reconstructed_mammoth_dna.txt'에 저장되었습니다.\n")

    end_time = time.time()
    print(f"⏱ 전체 실행 시간: {end_time - start_time:.2f}초")


if __name__ == "__main__":
    main()