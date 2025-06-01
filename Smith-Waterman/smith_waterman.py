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
                score[i - 1][j - 1] + match,  # 대각선 (match/mismatch)
                score[i - 1][j] + gap_penalty,  # 위 (deletion)
                score[i][j - 1] + gap_penalty   # 왼쪽 (insertion)
            )

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    return max_score, max_pos


def seed_and_extend(reference, read, seed_len=10):
    seed = read[:seed_len]
    positions = []

    # reference에서 seed와 일치하는 위치 찾기
    for i in range(len(reference) - seed_len + 1):
        if reference[i:i + seed_len] == seed:
            positions.append(i)

    best_score = -1
    best_pos = -1

    for pos in positions:
        # read보다 조금 더 긴 구간을 비교 대상으로 설정
        window = reference[pos : pos + len(read) + 20]
        score, _ = smith_waterman(read, window)

        if score > best_score:
            best_score = score
            best_pos = pos

    return best_pos, best_score

def evaluate_accuracy(true_positions, predicted_positions, tolerance=5):
    correct = 0
    total = min(len(true_positions), len(predicted_positions))
    for i in range(total):
        if abs(true_positions[i] - predicted_positions[i]) <= tolerance:
            correct += 1
    accuracy = (correct / total) * 100
    return accuracy, correct, total

def load_reference(filename="reference_1M.txt"):
    with open(filename, "r") as f:
        return f.read().strip()

def load_reads(filename="mammoth_reads_10K.txt"):
    with open(filename, "r") as f:
        return [line.strip() for line in f.readlines()]

def load_ground_truth(filename="ground_truth_10K.txt"):
    with open(filename, "r") as f:
        return [int(line.strip()) for line in f.readlines()]

# 1. 파일 로딩
reference = load_reference("reference_1M.txt")
reads = load_reads("mammoth_reads_10K.txt")
true_positions = load_ground_truth("ground_truth_10K.txt")

# 2. 매칭
predicted_positions = []
for i in range(1000):
    pred_pos, _ = seed_and_extend(reference, reads[i])
    predicted_positions.append(pred_pos)

# 3. 정확도 평가
accuracy, correct, total = evaluate_accuracy(true_positions[:100], predicted_positions)
print(f"\n✅ Accuracy: {accuracy:.2f}% ({correct}/{total} matched within ±5bp)")