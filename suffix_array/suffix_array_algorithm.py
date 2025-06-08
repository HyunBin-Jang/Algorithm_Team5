import time
import functools # 비교 함수를 정렬에 사용할 수 있게 변환하는 데 사용 (cmp_to_key)

# 주어진 문자열의 접미사를 사전 순으로 정렬한 suffix array 생성
def build_suffix_array(text: str) -> list:
    # i번째 접미사와 j번째 접미사를 비교
    def compare(i, j):
        while i < len(text) and j < len(text):
            # 문자 다르면 사전순 비교
            if text[i] != text[j]:
                return -1 if text[i] < text[j] else 1
            i += 1
            j += 1
        # 한쪽이 먼저 끝난 경우 비교 결과 처리
        if i == len(text) and j != len(text):
            return -1
        elif j == len(text) and i != len(text):
            return 1
        else:
            return 0

    # 접미사 시작 위치(index)들을 사전 순으로 정렬
    return sorted(range(len(text)), key=functools.cmp_to_key(compare))

# pattern과 reference의 특정 위치 substring을 사전 순으로 비교
def compare_pattern_with_mammoth(text: str, pos: int, pattern: str) -> int:
    # reference의 해당 위치에서 pattern 길이만큼 잘라서 비교
    ref_chunk = text[pos : pos + len(pattern)]
    if ref_chunk == pattern:
        return 0
    return -1 if ref_chunk < pattern else 1

# suffix array를 이용해 정확히 일치하는 pattern의 위치 리스트 반환
def search_exact_matches_in_mammoth(sa: list, text: str, pattern: str) -> list:
    n = len(sa)

    # lower bound (pattern 이상인 첫 위치) 찾기
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        if compare_pattern_with_mammoth(text, sa[mid], pattern) < 0:
            left = mid + 1
        else:
            right = mid
    lower = left

    # upper bound (pattern 초과인 첫 위치) 찾기
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        if compare_pattern_with_mammoth(text, sa[mid], pattern) <= 0:
            left = mid + 1
        else:
            right = mid
    upper = left

    # 정렬된 suffix array에서 pattern과 정확히 일치하는 위치들 반환
    return [sa[i] for i in range(lower, upper)]

# 두 문자열 간의 mismatch 개수 반환
def count_mismatches_in_bases(s1: str, s2: str) -> int:
    # 같은 위치의 문자가 다를 때 count
    return sum(1 for a, b in zip(s1, s2) if a != b)

# 하나의 read를 reference에 정렬, 최대 k mismatch 허용
def align_read_to_mammoth_reference(sa: list, reference: str, read: str, k: int) -> int:
    read_len = len(read)
    num_chunks = k + 1  # pigeonhole 원리: 하나는 무조건 맞아야 하므로 k+1 개로 쪼갬
    chunk_size = read_len // num_chunks
    candidates = set()  # 후보 시작 위치 저장

    # read를 여러 chunk로 나누고, 각 chunk에 대해 exact match 검색
    for i in range(num_chunks):
        chunk_start = i * chunk_size
        chunk_end = read_len if i == num_chunks - 1 else (i + 1) * chunk_size
        chunk = read[chunk_start:chunk_end]
        if not chunk:
            continue

        # suffix array에서 exact match 위치 찾기
        matches = search_exact_matches_in_mammoth(sa, reference, chunk)

        for match_pos in matches:
            # read 전체가 시작될 수 있는 위치 환산
            candidate_pos = match_pos - chunk_start
            if candidate_pos < 0 or candidate_pos + read_len > len(reference):
                continue
            candidates.add(candidate_pos)

    # 후보 위치들 중에서 실제 mismatch 가장 적은 위치 선택
    best_pos = -1
    min_mismatches = k + 1  # 초기값은 허용 범위 밖
    for pos in candidates:
        segment = reference[pos : pos + read_len]
        mismatch_count = count_mismatches_in_bases(segment, read)
        if mismatch_count < min_mismatches:
            min_mismatches = mismatch_count
            best_pos = pos
            # 완전히 일치하면 바로 종료
            if mismatch_count == 0:
                break

    # best_pos 반환 (k 이내 mismatch일 경우만)
    return best_pos if min_mismatches <= k else -1

# reference를 블록 단위로 나누어 read 전체 정렬 수행
def align_all_reads_to_mammoth_blocks(reference: str, reads: list, k: int, block_size: int = 500_000, overlap: int = 100) -> dict:
    total_reads = len(reads)
    alignments = {}  # 결과 저장: read index -> 위치
    matched_flags = [False] * total_reads  # 이미 매칭된 read는 생략

    # reference를 block_size 크기로 나누되, 블록 경계에는 overlap 추가
    total_blocks = (len(reference) + block_size - 1) // block_size

    for block_index in range(total_blocks):
        block_start = block_index * block_size
        block_end = min((block_index + 1) * block_size + overlap, len(reference))
        ref_block = reference[block_start:block_end]

        print(f"Block {block_index + 1}/{total_blocks} (위치 {block_start}~{block_end}) 처리 중")

        # 현재 블록의 suffix array 생성
        sa = build_suffix_array(ref_block)

        for i, read in enumerate(reads):
            if matched_flags[i]:
                continue

            # 해당 read를 현재 블록에 정렬 시도
            pos_in_block = align_read_to_mammoth_reference(sa, ref_block, read, k)
            if pos_in_block >= 0:
                absolute_pos = block_start + pos_in_block
                if absolute_pos + len(read) <= len(reference):
                    alignments[i] = absolute_pos
                    matched_flags[i] = True

    # 정렬 실패한 read는 -1로 처리
    for i in range(total_reads):
        if not matched_flags[i]:
            alignments[i] = -1

    return alignments

# 정렬된 read들을 reference에 덮어써서 최종 서열 복원
def rebuild_mammoth_with_aligned_reads(reference: str, reads: list, alignments: dict) -> str:
    reference_list = list(reference)
    if not reads:
        return reference

    read_len = len(reads[0])

    for i, read in enumerate(reads):
        pos = alignments.get(i, -1)
        if pos < 0:
            continue
        # reference의 해당 위치를 read로 덮어씀
        for j, base in enumerate(read):
            if pos + j < len(reference_list):
                reference_list[pos + j] = base

    return "".join(reference_list)

# 전체 정렬 및 복원 작업 흐름
def main():
    max_mismatches = 2           # 허용 mismatch 개수
    BLOCK_SIZE = 500_000         # reference 블록 크기
    OVERLAP = 100                # 블록 경계 겹침

    start_time = time.time()

    # reference, read, 정답 위치 로드
    with open("../genome_generation/3_1_reference_1M.txt", "r") as f_ref:
        reference = f_ref.read().strip()
    with open("../genome_generation/3_1_mammoth_reads_100K.txt", "r") as f_reads:
        reads = [line.strip() for line in f_reads if line.strip()]
    with open("../genome_generation/3_1_ground_truth_100K.txt", "r") as f_gt:
        ground_truth = [int(line.strip()) for line in f_gt if line.strip()]

    # 전체 read 정렬
    alignments = align_all_reads_to_mammoth_blocks(reference, reads, max_mismatches, block_size=BLOCK_SIZE, overlap=OVERLAP)

    # 정확도 평가
    total_reads = len(reads)
    correct_matches = 0
    for i, true_pos in enumerate(ground_truth):
        predicted_pos = alignments.get(i, -1)
        if predicted_pos == true_pos:
            correct_matches += 1
    accuracy = correct_matches / total_reads * 100

    print(f"Alignment accuracy: {accuracy:.2f}%\n")

    # read로 reference 복원
    reconstructed = rebuild_mammoth_with_aligned_reads(reference, reads, alignments)

    # 복원 결과 저장
    with open("reconstructed_mammoth_dna.txt", "w") as f_out:
        f_out.write(reconstructed)
    print("복원된 mammoth DNA 시퀀스 저장 완료\n")

    # 각 read의 정렬 위치 출력
    with open("read_alignment_positions.txt", "w") as f_pos:
        for i in range(total_reads):
            f_pos.write(f"{alignments[i]}\n")

    end_time = time.time()
    print(f"전체 실행 시간: {end_time - start_time:.2f}초")

# 메인 실행 지점
if __name__ == "__main__":
    main()
