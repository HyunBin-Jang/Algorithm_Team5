import time

def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    m, n = len(seq1), len(seq2)

    # 점수 테이블을 0으로 초기화
    score = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0
    max_pos = (0, 0)

    # DP 테이블 채우기
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match = match_score
            else:
                match = mismatch_penalty

            score[i][j] = max(
                0,
                score[i - 1][j - 1] + match,
                score[i - 1][j] + gap_penalty,
                score[i][j - 1] + gap_penalty
            )

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    return max_score, max_pos

def seed_and_extend(reference, read, kmer_index, seed_len=10):
    L = len(read)
    ref_len = len(reference)
    start_idx = L // 2 - seed_len // 2
    seed = read[start_idx:start_idx + seed_len]
    positions = kmer_index.get(seed, [])

    best_score = -1
    best_pos = -1

    for pos in positions:
        win_start = max(0, pos - start_idx)
        win_end = min(ref_len, win_start + len(read))
        window = reference[win_start : win_end]
        score, _ = smith_waterman(read, window)
        if score > best_score:
            best_score = score
            best_pos = pos - start_idx

    return best_pos, best_score

def build_kmer_index(reference, k):
    """
    reference 문자열에서 k-mer 인덱스를 딕셔너리 형태로 생성
    key: k-mer 문자열
    value: [등장 위치 인덱스 리스트]
    """
    index = {}
    for i in range(len(reference) - k + 1):
        kmer = reference[i:i + k]
        if kmer in index:
            index[kmer].append(i)
        else:
            index[kmer] = [i]
    return index

def evaluate_accuracy(true_positions, predicted_positions):
    correct = 0
    total = min(len(true_positions), len(predicted_positions))
    for i in range(total):
        if true_positions[i] == predicted_positions[i]:
            correct += 1
    accuracy = (correct / total) * 100
    return accuracy, correct, total

def load_reference(filename="reference_10M.txt"):
    with open(filename, "r") as f:
        return f.read().strip()

def load_reads(filename="mammoth_reads_100K.txt"):
    with open(filename, "r") as f:
        return [line.strip() for line in f.readlines()]

def load_ground_truth(filename="ground_truth_100K.txt"):
    with open(filename, "r") as f:
        return [int(line.strip()) for line in f.readlines()]

# 1. 파일 로딩
reference = load_reference("reference_100M.txt")
reads = load_reads("mammoth_reads_1M.txt")
true_positions = load_ground_truth("ground_truth_1M.txt")

# 2. 매칭
start_time = time.time()   # 매칭 시작 시간 기록
k = 20
kmer_index = build_kmer_index(reference, k)
predicted_positions = []

for i in range(10000):
    pred_pos, _ = seed_and_extend(reference, reads[i], kmer_index, seed_len=k)
    predicted_positions.append(pred_pos)
end_time = time.time()  #매칭 종료 시간 기록
elapsed_time = end_time - start_time
print(f"Total Matching Time : {elapsed_time:.2f} seconds")

# 3. 정확도 평가
accuracy, correct, total = evaluate_accuracy(true_positions, predicted_positions)
print(f"\n Accuracy: {accuracy:.2f}% ({correct}/{total} matched)")