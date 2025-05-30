import os

DATA_DIR = "../genome_generation"

def load_reference(filename="reference_1M.txt"):
    with open(os.path.join(DATA_DIR, filename), "r") as f:
        return f.read().strip()

def load_reads(filename="mammoth_reads_10K.txt"):
    with open(os.path.join(DATA_DIR, filename), "r") as f:
        return [line.strip() for line in f]

def load_ground_truth(filename="ground_truth_10K.txt"):
    with open(os.path.join(DATA_DIR, filename), "r") as f:
        return [int(line.strip()) for line in f]

# 🧬 1. 파일 불러오기
reference = load_reference()
reads = load_reads()
true_positions = load_ground_truth()

# 🧬 2. 코끼리 reference를 복사해서 메머드 복원용 배열 생성
reconstructed = list(reference)

# 🧬 3. 각 read를 해당 위치에 복원
for read, pos in zip(reads, true_positions):
    for i in range(len(read)):
        if pos + i < len(reconstructed):
            reconstructed[pos + i] = read[i]  # 덮어쓰기 방식

# 🧬 4. 최종 복원된 메머드 DNA 생성
mammoth_dna = "".join(reconstructed)

# 🧬 5. 결과 저장
with open("reconstructed_mammoth_dna.txt", "w") as f:
    f.write(mammoth_dna)

print("✅ 메머드 DNA 복원 완료 → reconstructed_mammoth_dna.txt")