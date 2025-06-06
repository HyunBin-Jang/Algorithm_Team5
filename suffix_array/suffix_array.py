import time
import functools


def build_suffix_array(text: str) -> list:
    """
    슬라이스 복사 없이 비교하는 메모리 절약형 suffix array 생성
    """
    def compare(i, j):
        while i < len(text) and j < len(text):
            if text[i] != text[j]:
                return -1 if text[i] < text[j] else 1
            i += 1
            j += 1
        return -1 if i == len(text) and j != len(text) else (1 if j == len(text) and i != len(text) else 0)

    return sorted(range(len(text)), key=functools.cmp_to_key(compare))


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


def align_all_reads_blockwise(reference: str, reads: list, k: int, block_size: int = 500_000, overlap: int = 100) -> dict:
    """
    reference를 block_size 단위로 나누고, 블록 간 overlap을 추가하여 suffix array 생성 및 매칭
    """
    total_reads = len(reads)
    alignments = {}
    matched = [False] * total_reads

    num_blocks = (len(reference) + block_size - 1) // block_size

    for block_idx in range(num_blocks):
        block_start = block_idx * block_size
        block_end = min((block_idx + 1) * block_size + overlap, len(reference))  # ← overlap 추가
        block_ref = reference[block_start:block_end]

        print(f"🔄 Block {block_idx + 1}/{num_blocks} (위치 {block_start}~{block_end}) 처리 중...")
        sa = build_suffix_array(block_ref)

        for idx, read in enumerate(reads):
            if matched[idx]:
                continue

            pos_in_block = align_single_read(sa, block_ref, read, k)
            if pos_in_block >= 0:
                global_pos = block_start + pos_in_block
                if global_pos + len(read) <= len(reference):  # 전체 범위 초과 방지
                    alignments[idx] = global_pos
                    matched[idx] = True

    # 매칭 실패한 read는 -1로 저장
    for idx in range(total_reads):
        if not matched[idx]:
            alignments[idx] = -1

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
            if pos + i < len(ref_list):
                ref_list[pos + i] = base

    return "".join(ref_list)


def main():
    k = 2  # 허용 mismatch 개수
    BLOCK_SIZE = 500_000  # reference 블록 크기
    OVERLAP = 100  # 블록 간 오버랩

    start_time = time.time()

    with open("../genome_generation/3_3_reference_100M.txt", "r") as f_ref:
        reference = f_ref.read().strip()

    with open("../genome_generation/3_3_mammoth_reads_100K.txt", "r") as f_reads:
        reads = [line.strip() for line in f_reads if line.strip()]

    with open("../genome_generation/3_3_ground_truth_100K.txt", "r") as f_gt:
        ground_truth = [int(line.strip()) for line in f_gt if line.strip()]

    alignments = align_all_reads_blockwise(reference, reads, k, block_size=BLOCK_SIZE, overlap=OVERLAP)

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