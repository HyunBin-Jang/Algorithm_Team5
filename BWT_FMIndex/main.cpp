#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <chrono>

using namespace std;


// 접미사 배열 생성
vector<int> buildSuffixArray(const string& s) {
	int n = s.length();
	vector<int> sa(n); // 접미사의 시작 인덱스를 저장하는 배열
	vector<int> rank(n);
	/* rank[i]: 접미사 s[i:]의 앞 k글자에 해당하는 정렬 순위
	 * rank[i + k]: 그 다음 k글자 (s[i+k : i+2k])의 정렬 순위
	 * (rank[i], rank[i + k])는 접미사 s[i:]의 앞 2k글자"를 두 부분으로 쪼개서 비교한 것과 동일 */
	vector<int> tmp(n);

	for (int i = 0; i < n; i++) {
		sa[i] = i; // [0,1,2,3,4,5,6]
		rank[i] = s[i]; // rank = [66, 65, 78, 65, 78, 65, 36] BANANA$
	}

	for (int k = 1; k < n; k *= 2) { // 접미사의 왼쪽 끝부터 2k글자씩 비교
		auto cmp = [&](int i, int j) {
			// 앞 글자가 다르면 그 순위로 비교
			if (rank[i] != rank[j]) return rank[i] < rank[j];
			// 앞 글자가 같으면 다음 k글자 비교
			int ri = i + k < n ? rank[i + k] : -1;
			int rj = j + k < n ? rank[j + k] : -1;
			return ri < rj;
			};
		sort(sa.begin(), sa.end(), cmp);

		tmp[sa[0]] = 0;
		for (int i = 1; i < n; i++) {
			tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
			/* 이전 접미사와 같으면 같은 rank, 다르면 + 1 증가
			 * 이전 접미사가 (1, 2)이고 현재 접미사가 (1, 2)라면 동일한 순위를 가짐
			 * tmp[sa[i]] : sa[i]번째부터 시작하는 접미사의 순위. ex) temp[sa[2]] = 3의 뜻은 2번째부터 시작하는 접미사의 순위가 3위라는 뜻 */
		}
		rank = tmp;
		// k = 4일 때 앞의 4글자씩을 비교해야하는데 이미 k = 2일 때 앞 두 글자의 상대 순위를 구해놨으므로
	}
	return sa;
}

// BWT 생성
string buildBWT(const string& s, const vector<int>& sa) {
	string bwt;
	int n = s.length();
	for (int i = 0; i < n; i++) {
		bwt += (sa[i] == 0) ? s[n - 1] : s[sa[i] - 1];
	}
	return bwt;
}

// FM-Index 클래스
class FMIndex {
private:
	string bwt;
	vector<int> sa;
	map<char, int> first; // ACGT가 원래 텍스트의 접미사들을 정렬한 것의 첫 열에서 처음 나타나는 위치. LF 매핑할 때 시작점으로 사용됨
	vector<map<char, int>> occ; // BWT의 각 위치에서 ACGT의 누적 개수. occ[i][c]는 BWT[0:i]에서 c가 나타난 횟수

public:
	FMIndex(const string& ref) { // ref : reference sequence
		sa = buildSuffixArray(ref + '$');
		//for (int i = 0; i < sa.size(); i++)
		//	cout << "Suffix Array[" << i << "]: " << sa[i] << endl; // 디버깅용 출력
		bwt = buildBWT(ref + '$', sa);
		//cout << "BWT: " << bwt << endl; // 디버깅용 출력

		map<char, int> char_count;
		for (char c : bwt) char_count[c]++;
		int sum = 0;
		for (auto& p : char_count) {
			first[p.first] = sum;
			sum += p.second;
		}
		//for (auto& p : first)
		//	cout << "First[" << p.first << "]: " << p.second << endl; // 디버깅용 출력

		occ.resize(bwt.length());
		map<char, int> running_count;
		for (size_t i = 0; i < bwt.length(); i++) {
			running_count[bwt[i]]++;
			occ[i] = running_count;
		}
		//for (size_t i = 0; i < occ.size(); i++) {
		//	cout << "occ[" << i << "]: ";
		//	for (const auto& p : occ[i]) {
		//		cout << p.first << ":" << p.second << " ";
		//	}
		//	cout << endl; // 디버깅용 출력
		//}
	}

	// 정확한 매칭을 위한 범위 계산. 패턴의 한 문자를 처리하여 BWT에서 해당 문자의 범위를 반환
	bool getRange(char c, int& top, int& bottom) { // 현재 BWT에서 검색 범위 : top ~ bottom
		if (first.find(c) == first.end()) return false; // c가 BWT에 없으면 false 반환
		top = first[c] + (top > 0 ? occ[top - 1][c] : 0); // c의 첫 번째 위치에 top 이전 까지의 c 등장 횟수 더하기
		bottom = first[c] + occ[bottom][c] - 1;
		return top <= bottom;
	}

	// 패턴을 뒤에서부터 처리. 최대 k개의 mismatch를 허용하여 매칭 위치를 찾음.
	void approxSearch(const string& pattern, int k, int pos, int top, int bottom, vector<int>& results) {
		/* pattern: 검색할 패턴.
		 * k: 허용할 최대 mismatch 수.
		 * pos: 현재 처리 중인 패턴의 인덱스(뒤에서부터 처리).
		 * top, bottom: 현재 BWT에서의 검색 범위.
		 * results: 매칭된 위치를 저장할 벡터. */

		if (pos < 0) { // 패턴을 모두 처리한 경우.
			for (int i = top; i <= bottom; i++) {
				if (sa[i] == pattern.length()) continue; // '$' 제외
				results.push_back(sa[i]);
			}
			return;
		}

		// 정확한 매칭 시도
		int new_top = top, new_bottom = bottom;
		if (getRange(pattern[pos], new_top, new_bottom)) {
			approxSearch(pattern, k, pos - 1, new_top, new_bottom, results);
		}

		// k가 남아있으면 mismatch 시도
		if (k > 0) {
			for (char c : {'A', 'C', 'G', 'T'}) {
				if (c == pattern[pos]) continue;
				new_top = top;
				new_bottom = bottom;
				if (getRange(c, new_top, new_bottom)) {
					approxSearch(pattern, k - 1, pos - 1, new_top, new_bottom, results);
				}
			}
		}
	}

	vector<int> searchWithMismatch(const string& pattern, int k) {
		vector<int> results;
		approxSearch(pattern, k, pattern.length() - 1, 0, bwt.length() - 1, results);
		sort(results.begin(), results.end());
		return results;
	}
};

// 파일 읽기 함수
string readReference(const string& filename) {
	ifstream file(filename);
	string ref;
	getline(file, ref);
	file.close();
	return ref;
}

vector<string> readPatterns(const string& filename) {
	ifstream file(filename);
	vector<string> patterns;
	string line;
	while (getline(file, line)) {
		patterns.push_back(line);
	}
	file.close();
	return patterns;
}

vector<int> readGroundTruth(const string& filename) {
	ifstream file(filename);
	vector<int> truth;
	int pos;
	while (file >> pos) {
		truth.push_back(pos);
	}
	file.close();
	return truth;
}

string formatDuration(chrono::milliseconds ms) {
	using namespace std::chrono;

	auto total_seconds = duration_cast<seconds>(ms).count();
	int hours = total_seconds / 3600;
	int minutes = (total_seconds % 3600) / 60;
	int seconds = total_seconds % 60;

	std::string result;
	if (hours > 0) {
		result += std::to_string(hours) + "시간 ";
	}
	result += std::to_string(minutes) + "분 ";
	result += std::to_string(seconds) + "초";

	return result;
}

void FMIndexBWT(const string referenceSize, const string patternSize) {
	cout << "===============================" << endl;
	cout << "Reference Size: " << referenceSize << ", Number of patterns : " << patternSize << endl;
	auto start = chrono::high_resolution_clock::now();

	// 파일 읽기
	string ref = readReference("reference_" + referenceSize + ".txt");
	vector<string> patterns = readPatterns("short_reads/mammoth_reads_" + referenceSize + "_" + patternSize + ".txt");
	vector<int> ground_truth = readGroundTruth("ground_truth/ground_truth_" + referenceSize + "_" + patternSize + ".txt");

	// FM-Index 생성
	FMIndex fm(ref);

	// 최대 mismatch 허용 수
	int k = 2; // 필요시 조정 가능

	// 결과 출력 파일
	/*vector<string> outputs;
	ofstream out("results_" + referenceSize + "_" + patternSize + ".txt");*/

	// 패턴 검색 및 검증
	int correct = 0;
	for (size_t i = 0; i < patterns.size(); i++) {
		vector<int> positions = fm.searchWithMismatch(patterns[i], k);

		/*for (int i = 0; i < positions.size(); i++)
			cout << positions[i] << " ";
		cout << endl;*/

		bool found = false;
		for (int pos : positions) {
			if (pos == ground_truth[i]) {
				found = true;
				break;
			}
		}
		if (found) correct++;
		//outputs.push_back(("Pattern " + to_string(i + 1) + ": " + (found ? "Match at " + to_string(ground_truth[i]) : "No match")));
	}
	auto end = chrono::high_resolution_clock::now();

	/*for (const string& output : outputs) {
		out << output << endl;
	}*/
	cout << "Accuracy: " << (double)correct / patterns.size() * 100 << "%" << endl;

	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
	string executionTime = formatDuration(duration);
	cout << "Execution time: " << executionTime << endl;

	//out.close();
}

int main()
{
	string references[] = { "1M", "10M", "100M" };
	string patterns[] = { "10K", "100K", "1M", "10M" };
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			FMIndexBWT(references[i], patterns[j]);
		}
	}
	return 0;
}