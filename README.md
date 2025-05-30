# Algorithm_Team5 - 유전체 복원 알고리즘 프로젝트

알고리즘 팀 프로젝트 (5조)  
메머드의 read 데이터와 코끼리의 reference genome을 기반으로,  
다양한 알고리즘을 활용해 메머드 DNA를 복원합니다.

---

## 📁 프로젝트 구조
Algorithm_Team5/
├── genome_generation/           # 데이터 생성 및 전처리 코드
├── myers_edit_distance/         # Myers 알고리즘 구현 및 실행
├── venv/                        # Python 가상환경 (추적 제외)
├── requirements.txt             # 필요 패키지 목록
├── .gitignore
└── README.md

---

## 📌 사용한 알고리즘 (브랜치별)

| 알고리즘 | 설명 | 브랜치 이름 |
|----------|------|--------------|
| Myers' Bit-parallel Edit Distance | 비트 연산 기반의 빠른 근사 문자열 매칭 | `feat/myers-algorithm` |
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