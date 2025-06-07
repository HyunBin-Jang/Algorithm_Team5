#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <chrono>

using namespace std;

vector<int> buildSuffixArray(const string& s) {
	int n = s.length();
	vector<int> sa(n); // ���̻��� ���� �ε����� �����ϴ� �迭
	vector<int> rank(n);
	/* rank[i]: ���̻� s[i:]�� �� k���ڿ� �ش��ϴ� ���� ����
	 * rank[i + k]: �� ���� k���� (s[i+k : i+2k])�� ���� ����
	 * (rank[i], rank[i + k])�� ���̻� s[i:]�� �� 2k����"�� �� �κ����� �ɰ��� ���� �Ͱ� ���� */
	vector<int> tmp(n);

	for (int i = 0; i < n; i++) {
		sa[i] = i; // [0,1,2,3,4,5,6]
		rank[i] = s[i]; // rank = [66, 65, 78, 65, 78, 65, 36] BANANA$
	}

	for (int k = 1; k < n; k *= 2) { // ���̻��� ���� ������ 2k���ھ� ��
		auto cmp = [&](int i, int j) {
			// �� ���ڰ� �ٸ��� �� ������ ��
			if (rank[i] != rank[j]) return rank[i] < rank[j];
			// �� ���ڰ� ������ ���� k���� ��
			int ri = i + k < n ? rank[i + k] : -1;
			int rj = j + k < n ? rank[j + k] : -1;
			return ri < rj;
			};
		sort(sa.begin(), sa.end(), cmp);

		tmp[sa[0]] = 0;
		for (int i = 1; i < n; i++) {
			tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
			/* ���� ���̻�� ������ ���� rank, �ٸ��� + 1 ����
			 * ���� ���̻簡 (1, 2)�̰� ���� ���̻簡 (1, 2)��� ������ ������ ����
			 * tmp[sa[i]] : sa[i]��°���� �����ϴ� ���̻��� ����. ex) temp[sa[2]] = 3�� ���� 2��°���� �����ϴ� ���̻��� ������ 3����� �� */
		}
		rank = tmp;
		// k = 4�� �� ���� 4���ھ��� ���ؾ��ϴµ� �̹� k = 2�� �� �� �� ������ ��� ������ ���س����Ƿ�
	}
	return sa;
}

// BWT ����
string buildBWT(const string& s, const vector<int>& sa) {
	string bwt;
	int n = s.length();
	for (int i = 0; i < n; i++) {
		bwt += (sa[i] == 0) ? s[n - 1] : s[sa[i] - 1];
	}
	return bwt;
}

// FM-Index Ŭ����
class FMIndex {
private:
	string bwt;
	vector<int> sa;
	map<char, int> first; // ACGT�� ���� �ؽ�Ʈ�� ���̻���� ������ ���� ù ������ ó�� ��Ÿ���� ��ġ. LF ������ �� ���������� ����
	vector<map<char, int>> occ; // BWT�� �� ��ġ���� ACGT�� ���� ����. occ[i][c]�� BWT[0:i]���� c�� ��Ÿ�� Ƚ��

public:
	FMIndex(const string& ref) { // ref : reference sequence
		sa = buildSuffixArray(ref + '$');
		bwt = buildBWT(ref + '$', sa);

		map<char, int> char_count;
		for (char c : bwt) char_count[c]++;
		int sum = 0;
		for (auto& p : char_count) {
			first[p.first] = sum;
			sum += p.second;
		}

		occ.resize(bwt.length());
		map<char, int> running_count;
		for (size_t i = 0; i < bwt.length(); i++) {
			running_count[bwt[i]]++;
			occ[i] = running_count;
		}
	}

	// ��Ȯ�� ��Ī�� ���� ���� ���. ������ �� ���ڸ� ó���Ͽ� BWT���� �ش� ������ ������ ��ȯ
	bool getRange(char c, int& top, int& bottom) { // ���� BWT���� �˻� ���� : top ~ bottom
		if (first.find(c) == first.end()) return false; // c�� BWT�� ������ false ��ȯ
		top = first[c] + (top > 0 ? occ[top - 1][c] : 0); // c�� ù ��° ��ġ�� top ���� ������ c ���� Ƚ�� ���ϱ�
		bottom = first[c] + occ[bottom][c] - 1;
		return top <= bottom;
	}

	// ������ �ڿ������� ó��. �ִ� k���� mismatch�� ����Ͽ� ��Ī ��ġ�� ã��.
	void approxSearch(const string& pattern, int k, int pos, int top, int bottom, vector<int>& results) {
		/* pattern: �˻��� ����.
		 * k: ����� �ִ� mismatch ��.
		 * pos: ���� ó�� ���� ������ �ε���(�ڿ������� ó��).
		 * top, bottom: ���� BWT������ �˻� ����.
		 * results: ��Ī�� ��ġ�� ������ ����. */

		if (pos < 0) { // ������ ��� ó���� ���.
			for (int i = top; i <= bottom; i++) {
				if (sa[i] == pattern.length()) continue; // '$' ����
				results.push_back(sa[i]);
			}
			return;
		}

		// ��Ȯ�� ��Ī �õ�
		int new_top = top, new_bottom = bottom;
		if (getRange(pattern[pos], new_top, new_bottom)) {
			approxSearch(pattern, k, pos - 1, new_top, new_bottom, results);
		}

		// k�� ���������� mismatch �õ�
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

// ���� �б� �Լ�
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
		result += std::to_string(hours) + "�ð� ";
	}
	result += std::to_string(minutes) + "�� ";
	result += std::to_string(seconds) + "��";

	return result;
}

void FMIndexBWT(const string ref, const vector<string> patterns, const vector<int> ground_truth, const int readLength) {
	cout << "============================================" << endl;
	cout << "reference length: " << ref.length() << ", patterns length: " << readLength << ", number of patterns: " << patterns.size() << endl;
	auto start = chrono::high_resolution_clock::now();

	// ���� �б�

	// FM-Index ����
	FMIndex fm(ref);

	// �ִ� mismatch ��� ��
	int k = 2; // �ʿ�� ���� ����

	// ���� �˻� �� ����
	int correct = 0;
	size_t total = patterns.size();
	size_t nextProgress = total / 10;
	size_t progressCheck = nextProgress;

	for (size_t i = 0; i < total; i++) {
		vector<int> positions = fm.searchWithMismatch(patterns[i], k);

		bool found = false;
		for (int pos : positions) {
			if (pos == ground_truth[i]) {
				found = true;
				break;
			}
		}
		if (found) correct++;
		if (i + 1 >= progressCheck) {
			cout << (progressCheck * 100 / total) << "% completed..." << endl;
			progressCheck += nextProgress;
		}
	}
	auto end = chrono::high_resolution_clock::now();

	cout << "Accuracy: " << (double)correct / patterns.size() * 100 << "%" << endl;

	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
	string executionTime = formatDuration(duration);
	cout << "Execution time: " << executionTime << endl;
}

void first() {
	string sizes[] = { "10K", "100K", "1M" };
	string ref = readReference("1/1_reference_10M.txt");
	for (int i = 0; i < 3; i++) {
		vector<string> patterns = readPatterns("1/1_" + to_string(i + 1) + "_mammoth_reads_" + sizes[i] + ".txt");
		vector<int> ground_truth = readGroundTruth("1/1_" + to_string(i + 1) + "_ground_truth_" + sizes[i] + ".txt");
		FMIndexBWT(ref, patterns, ground_truth, 100);
	}
}

void second() {
	int sizes[] = { 40, 70, 100 };
	string ref = readReference("2/2_reference_10M.txt");
	for (int i = 0; i < 3; i++) {
		vector<string> patterns = readPatterns("2/2_" + to_string(i + 1) + "_mammoth_reads_100K.txt");
		vector<int> ground_truth = readGroundTruth("2/2_" + to_string(i + 1) + "_ground_truth_100K.txt");
		FMIndexBWT(ref, patterns, ground_truth, sizes[i]);
	}
}

void third() {
	string sizes[] = { "1M", "10M", "100M" };
	for (int i = 0; i < 3; i++) {
		string ref = readReference("3/3_" + to_string(i + 1) + "_reference_" + sizes[i] + ".txt");
		vector<string> patterns = readPatterns("3/3_" + to_string(i + 1) + "_mammoth_reads_100K.txt");
		vector<int> ground_truth = readGroundTruth("3/3_" + to_string(i + 1) + "_ground_truth_100K.txt");
		FMIndexBWT(ref, patterns, ground_truth, 100);
	}
}

int main()
{
	first();
	second();
	third();
	return 0;
}