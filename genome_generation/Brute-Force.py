def load_reference(filename="reference_10M.txt"):
    with open(filename, "r") as f:
        return f.read().strip()

def load_reads(filename="mammoth_reads_10K.txt"):
    with open(filename, "r") as f:
        return [line.strip() for line in f]

def load_ground_truth(filename="ground_truth_10K.txt"):
    with open(filename, "r") as f:
        return [int(line.strip()) for line in f]

def brute_force_match(reference, read, max_mismatch=2):
    best_index = -1
    best_mismatch = len(read) + 1
    for i in range(len(reference) - len(read) + 1):
        window = reference[i:i+len(read)]
        mismatches = sum(1 for a, b in zip(read, window) if a != b)
        if mismatches <= max_mismatch and mismatches < best_mismatch:
            best_index = i
            best_mismatch = mismatches

    return best_index, best_mismatch

# 파일 불러오기
reference = load_reference("reference_1M.txt")
reads = load_reads("mammoth_reads_10K.txt")
true_positions = load_ground_truth("ground_truth_10K.txt")

# 테스트 수 제한 (예: 상위 100개만)
test_n = 10000
correct = 0

print("read_id\ttrue_pos\tpredicted_pos\tmismatch")

for i in range(test_n):
    read = reads[i]
    true_pos = true_positions[i]
    pred_pos, mismatch = brute_force_match(reference, read, max_mismatch=2)

    if pred_pos == true_pos:
        correct += 1
    print(f"read_{i}\t{true_pos}\t{pred_pos}\t{mismatch}")

accuracy = correct / test_n * 100
print(f"\n✅ Accuracy (±5bp tolerance): {accuracy:.2f}%")
