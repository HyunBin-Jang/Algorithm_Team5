# Algorithm_Team5 - 유전체 복원 알고리즘 프로젝트

알고리즘 팀 프로젝트 (5조)  
메머드의 read 데이터와 코끼리의 reference genome을 기반으로,  
다양한 알고리즘을 활용해 메머드 DNA를 복원합니다.

---

```bash
## 📁 프로젝트 구조
Algorithm_Team5/
├── BWT-FMIndex/
├── genome_generation/
├── .gitignore
├── README.md
```

---

## 📌 사용한 알고리즘 (브랜치별)

| 알고리즘 | 설명 | 브랜치 이름 |
|----------|------|--------------|
| Myers' Bit-parallel Edit Distance | 비트 연산 기반의 빠른 근사 문자열 매칭 | `feat/myers-algorithm` |
| FM Index | 접미사 기반 BWT 생성 방식을 이용한 FMIndex 방식 | `feat/bwt` |
| Brute-force Matching | 모든 위치를 순차 비교하는 기본 매칭 방식 | `main` or `brute-force` |
| (추가 알고리즘) | ... | ... |

---

## ⚙️ 실행 방법

### 가상환경 설정

```bash
# 가상환경 생성
python3.11 -m venv venv

# 가상환경 활성화 (macOS/Linux)
source venv/bin/activate

# 패키지 설치
pip install -r requirements.txt

---

## 🧪 복원 정확도 평가 스크립트 사용법

`reconstructed_mammoth_dna.txt` 파일이 생성된 이후,  
복원이 얼마나 정확했는지 평가하려면 아래 스크립트를 실행하세요.

### ✅ 실행 방법

```bash
python evaluate_reconstruction.py
```

### ✅ 기능 설명

- `mammoth_reads_10K.txt`와 `ground_truth_10K.txt`를 불러와,
- 복원된 DNA(`reconstructed_mammoth_dna.txt`)에서 각 read가 **제 위치에 얼마나 정확히 복원되었는지** 평가합니다.
- **최대 mismatch 2개까지 허용**하여 read가 성공적으로 복원되었는지를 판단합니다.

### ✅ 출력 예시

```
✅ 허용 mismatch ≤ 2 기준 복원 성공 수: 9938 / 10000
🎯 복원 정확도 (mismatch ≤ 2): 99.38%
```
