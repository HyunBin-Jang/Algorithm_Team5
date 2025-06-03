# -*- coding: utf-8 -*-
"""
suffix array 기반으로 elephant(reference) DNA를 활용하여
mammoth_reads_100을 정위치 탐색하고,
ground truth와 비교하여 정확도를 계산한 뒤, 복원된 전체 맘모스 시퀀스를 출력하는 예시

실행에 걸린 시간과 함께, 중간중간 몇 개의 read를 처리했는지 진행 상황도 함께 출력합니다.
"""

import time  # 시간 측정을 위해 불러옴


def build_suffix_array(text: str) -> list:
    """
    reference 문자열(text)의 suffix array를 구성하여 반환
    """
    suffixes = [(text[i:], i) for i in range(len(text))]
    suffixes.sort(key=lambda x: x[0])
    sa = [pos for (_suffix, pos) in suffixes]
    return sa


def _compare_prefix(text: str, pos: int, pattern: str) -> int:
    """
    reference[text]의 pos 위치에서 pattern을 사전식으로 비교
    """
    slice_ref = text[pos : pos + len(pattern)]
    if slice_ref == pattern:
        return 0
    if slice_ref < pattern:
        return -1
    else:
        return 1


def find_exact_matches(sa: list, text: str, pattern: str) -> list:
    """
    suffix array(sa)와 reference(text)를 활용하여 pattern이 정확히 매칭되는
    모든 시작 인덱스를 반환.
    """
    n = len(sa)

    # lower bound: prefix >= pattern
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        pos_mid = sa[mid]
        cmp = _compare_prefix(text, pos_mid, pattern)
        if cmp < 0:
            left = mid + 1
        else:
            right = mid
    start = left

    # upper bound: prefix > pattern
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        pos_mid = sa[mid]
        cmp = _compare_prefix(text, pos_mid, pattern)
        if cmp <= 0:
            left = mid + 1
        else:
            right = mid
    end = left

    return [sa[i] for i in range(start, end)]


def count_mismatches(s1: str, s2: str) -> int:
    """
    두 문자열 s1, s2 길이는 같다고 가정. 서로 다른 문자 개수(mismatch)를 반환
    """
    return sum(1 for a, b in zip(s1, s2) if a != b)


def align_single_read(sa: list, reference: str, read: str, k: int) -> int:
    """
    하나의 read를 suffix array 기반으로 reference 위에 approximate matching
    허용 mismatch 개수 = k
    """
    L = len(read)
    t = k + 1
    chunk_size = L // t
    candidates = set()

    # read를 (k+1)개 chunk로 분할하여 정확히 일치하는 후보 위치 수집
    for i in range(t):
        start_i = i * chunk_size
        end_i = L if i == t - 1 else (i + 1) * chunk_size
        chunk = read[start_i:end_i]
        if len(chunk) == 0:
            continue

        matches = find_exact_matches(sa, reference, chunk)
        for match_pos in matches:
            cand_start = match_pos - start_i
            if cand_start < 0 or cand_start + L > len(reference):
                continue
            candidates.add(cand_start)

    # 후보 위치마다 전체 read와 비교하여 mismatch 개수 계산
    best_pos = -1
    best_mm = k + 1
    for pos in candidates:
        segment = reference[pos : pos + L]
        mm = count_mismatches(segment, read)
        if mm < best_mm:
            best_mm = mm
            best_pos = pos
            if best_mm == 0:
                break

    return best_pos if best_mm <= k else -1


def align_all_reads(reference: str, reads: list, k: int) -> dict:
    """
    전체 read들을 차례로 정위치 탐색하면서, 중간중간 진행 상황을 출력합니다.
    반환값: {read_index: aligned_start_pos} (매칭 실패 시 value = -1)
    """
    sa = build_suffix_array(reference)
    alignments = {}
    total_reads = len(reads)

    # 진행 상황을 출력할 간격 (예: 매 10개씩 출력)
    PRINT_INTERVAL = 10

    for idx, read in enumerate(reads):
        pos = align_single_read(sa, reference, read, k)
        alignments[idx] = pos

        # 진행 상황 출력
        if (idx + 1) % PRINT_INTERVAL == 0 or (idx + 1) == total_reads:
            percent = (idx + 1) / total_reads * 100
            print(f"🔄 처리 중: {idx + 1}/{total_reads} reads ({percent:.1f}%) 완료")

    return alignments


def reconstruct_mammoth(reference: str, reads: list, alignments: dict) -> str:
    """
    alignments 결과를 바탕으로 reference를 수정하여 reconstructed 시퀀스를 생성
    """
    ref_list = list(reference)
    L = len(reads[0]) if reads else 0

    for idx, read in enumerate(reads):
        pos = alignments.get(idx, -1)
        if pos is None or pos < 0:
            continue
        for i, base in enumerate(read):
            ref_list[pos + i] = base

    return "".join(ref_list)


def main():
    k = 2  # 허용 mismatch 개수

    # --- 타이머 시작 ---
    start_time = time.time()

    # 1) 파일 로드
    with open("../genome_generation/reference_100K.txt", "r") as f_ref:
        reference = f_ref.read().strip()

    with open("../genome_generation/mammoth_reads_1K.txt", "r") as f_reads:
        reads = [line.strip() for line in f_reads if line.strip()]

    # 2) ground truth 로드 (정확도 측정용)
    with open("../genome_generation/ground_truth_1K.txt", "r") as f_gt:
        ground_truth = [int(line.strip()) for line in f_gt if line.strip()]

    # 3) 모든 read 정위치 탐색 (진행 상황 출력 포함)
    alignments = align_all_reads(reference, reads, k)

    # 4) 정확도 계산 (ground truth와 예측 위치가 같은지만 판단)
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

    # 5) 복원된 시퀀스 생성
    reconstructed = reconstruct_mammoth(reference, reads, alignments)

    # 6) 파일로 저장
    with open("reconstructed_mammoth_dna.txt", "w") as f_out:
        f_out.write(reconstructed)

    print("✅ 복원된 mammoth DNA 시퀀스가 'reconstructed_mammoth_dna.txt'에 저장되었습니다.")

    # --- 타이머 종료 및 경과 시간 출력 ---
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"⏱ 전체 실행 시간: {elapsed:.2f}초")


if __name__ == "__main__":
    main()