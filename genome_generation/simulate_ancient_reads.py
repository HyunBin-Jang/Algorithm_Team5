import random

# reference txt 파일 불러오기
def load_reference(filename="reference_1M"):
    with open(filename, "r") as f:
        return f.read().strip()  # 줄바꿈 제거

# 고대 DNA read 생성 함수
def simulate_ancient_reads(reference_seq, read_len=100, num_reads=10000, mutation_rate=0.01):
    """
    reference_seq에서 M개의 read를 랜덤하게 자르고,
    고대 DNA 특성에 따라 read 말단에 C→T, G→A 변이를 삽입합니다.
    """
    reads = []
    truth = []
    print(len(reference_seq))
    for i in range(num_reads):
        start = random.randint(0, len(reference_seq) - read_len)
        original = list(reference_seq[start:start + read_len])
        mutated = []
        mutation_count = 0

        for j in range(read_len):
            prob = mutation_rate
            if j < 5 or j >= read_len - 5:
                prob *= 5  # 말단에서 변이 확률 증가

            base = original[j]
            if random.random() < prob:
                if base == 'C':
                    mutated.append('T')
                elif base == 'G':
                    mutated.append('A')
                else:
                    mutated.append(random.choice(['A', 'C', 'G', 'T']))
                mutation_count += 1
            else:
                mutated.append(base)

        read_str = ''.join(mutated)
        reads.append(read_str)
        truth.append((f'read_{i}', start, mutation_count))
    return reads, truth

# reference 불러오기
reference = load_reference("reference_1M.txt")

# read 생성
reads, ground_truth = simulate_ancient_reads(reference)

# read 저장
with open("mammoth_reads_10K.txt", "w") as f:
    for read in reads:
        f.write(read + "\n")

# ground_truth (index만) 저장
with open("ground_truth_10K.txt", "w") as f:
    for _, pos, _ in ground_truth:
        f.write(f"{pos}\n")