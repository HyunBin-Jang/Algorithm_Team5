def count_mismatches(s1, s2):
    """두 문자열 간 mismatch 수를 셈"""
    return sum(1 for a, b in zip(s1, s2) if a != b)

# 🧬 1. 파일 불러오기
with open("genome_generation/mammoth_reads_10K.txt") as f:
    reads = [line.strip() for line in f]

with open("genome_generation/ground_truth_10K.txt") as f:
    true_positions = [int(line.strip()) for line in f]

with open("myers_edit_distance/reconstructed_mammoth_dna.txt") as f:
    reconstructed = f.read().strip()

# 🎯 2. 정확도 측정 (mismatch ≤ 2 허용)
correct = 0
for read, pos in zip(reads, true_positions):
    segment = reconstructed[pos:pos + len(read)]
    mismatches = count_mismatches(read, segment)
    if mismatches <= 2:
        correct += 1

# ✅ 3. 정확도 출력
total = len(reads)
accuracy = correct / total * 100

print(f"✅ 허용 mismatch ≤ 2 기준 복원 성공 수: {correct} / {total}")
print(f"🎯 복원 정확도 (mismatch ≤ 2): {accuracy:.2f}%")