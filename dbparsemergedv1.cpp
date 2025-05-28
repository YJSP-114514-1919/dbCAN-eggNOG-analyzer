#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <regex>
#include <map>
#include <string>
#include <set>
#include <tuple>
#include <algorithm>

namespace fs = std::filesystem;
using namespace std;

// ===================== 原有 CazymeInfo 类及相关函数（不变） =====================
class CazymeInfo {
public:
    string strain_name;            // 菌株名
    map<string, int> cazyme_count; // CAZyme 出现次数

    CazymeInfo(string strain) : strain_name(move(strain)) {}

    void addCazyme(const string& cazyme) {
        cazyme_count[cazyme]++;
    }
};

int notfound = 0;
vector<string> notfoundstrains;

// 提取菌株名函数（不变）
string extractStrainName(const string& fna_file) {
    ifstream fna(fna_file);
    if (!fna.is_open()) {
        cerr << "Could not open .fna file: " << fna_file << endl;
        return "";
    }
    set<string> known_strains = {//已知的带空格菌株名
        "Lactobacillus crispatus strain ATCC 33820",
        "Lactobacillus crispatus strain DSM 108970",
        "Lactobacillus crispatus DSM 20584 = JCM 1185 = ATCC 33820 strain DSM 20584",
        "Lactobacillus crispatus DSM 20584 = JCM 1185 = ATCC 33820",
        "Lactobacillus crispatus strain LMG 11415",
        "Lactobacillus crispatus strain LMG 12005",
        "Lactobacillus crispatus strain LMG 18189",
        "Lactobacillus crispatus strain LMG 18199",
        "Lactobacillus crispatus strain LMG 11440",
        "Lactobacillus crispatus strain M247",
        "Lactobacillus crispatus strain LMG 18200",
        "Lactobacillus crispatus strain DSM 29598",
        "Lactobacillus crispatus strain CIRM-BIA 2233",
        "Lactobacillus crispatus strain CIP 104459",
        "Lactobacillus crispatus strain LMG 12003",
        "Lactobacillus crispatus strain CIRM-BIA 2111",
        "Lactobacillus crispatus strain CIRM-BIA 2112",
        "Lactobacillus crispatus strain CIRM-BIA 523",
        "Lactobacillus crispatus isolate metabat2_6755_HET_CTRL_metaflye.524 MAG",
        "Lactobacillus iners DSM 13335",
        "Lactobacillus iners ATCC 55195",
        "Lactobacillus iners LEAF 2062A-h1",
        "Lactobacillus iners UPII 143-D",
        "Lactobacillus iners UPII 60-B",
        "Lactobacillus iners SPIN 1401G",
        "Lactobacillus iners LEAF 3008A-a",
        "Lactobacillus iners LEAF 2053A-b",
        "Lactobacillus iners LEAF 2052A-d",
        "Lactobacillus iners SPIN 2503V10-D",
        "Lactobacillus iners LactinV 11V1-d",
        "Lactobacillus iners LactinV 09V1-c",
        "Lactobacillus iners LactinV 03V1-b",
        "Lactobacillus iners LactinV 01V1-a",
        "Lactobacillus sp. 7_1_47FAA"
    };

    string line, strain_name;
    regex strain_pattern1(R"(Lactobacillus\s+[A-Za-z]+\s+\S+)");
    regex strain_pattern2(R"(Lactobacillus\s+[A-Za-z]+\s+strain\s+\S+)");
    regex strain_pattern3(R"(Lactobacillus\s+[A-Za-z]+\s+isolate\s+\S+)");

    if (getline(fna, line)) {
        smatch match;
        // 先检查已知菌株名
        for (const auto& strain : known_strains) {
            if (line.find(strain) != string::npos) {
                std::cout << "Known strain matched: " << strain << endl;
                strain_name = strain;
                break;
            }
        }
        // 如果没找到已知菌株名，再尝试正则
        if (strain_name.empty()) {
            if (regex_search(line, match, strain_pattern1)) {
                strain_name = match[0];
                std::cout << "strain_pattern1 matched" << endl;
            }
            if ((strain_name == "Lactobacillus crispatus strain" || strain_name == "Lactobacillus iners strain") &&
                regex_search(line, match, strain_pattern2)) {
                strain_name = match[0];
                std::cout << "strain_pattern2 matched" << endl;
            }
            if ((strain_name == "Lactobacillus crispatus isolate" || strain_name == "Lactobacillus iners isolate") &&
                regex_search(line, match, strain_pattern3)) {
                strain_name = match[0];
                std::cout << "strain_pattern3 matched" << endl;
            }
        }
    }
    if (!strain_name.empty()) {
        std::cout << "Strain_name: " << strain_name << " extracted" << endl;
    }
    else {
        std::cout << "Strain_name in fna file: " << fna_file << " not found" << endl;
        notfound++;
        notfoundstrains.push_back(fna_file);
    }
    fna.close();
    return strain_name;
}

// 处理 overview 文件（不变）
void processOverview(const string& overview_file, CazymeInfo& cazyme_info) {
    ifstream overview(overview_file);
    if (!overview.is_open()) {
        cerr << "Could not open overview.txt file: " << overview_file << endl;
        return;
    }
    std::cout << "Beginning to process the file: " << overview_file << endl;

    string line;
    const regex pattern("([A-Za-z]+\\d+)");
    const regex pattern_sub("([A-Za-z]+\\d+_[A-Za-z]*\\d+)");

    getline(overview, line); // 跳过第一行
    while (getline(overview, line)) {
        stringstream ss(line);
        string gene_id, ec, hmmer, dbcan_sub, diamond;
        int tools;
        const string EC = " EC:";
        getline(ss, gene_id, '\t');
        getline(ss, ec, '\t');
        getline(ss, hmmer, '\t');
        getline(ss, dbcan_sub, '\t');
        getline(ss, diamond, '\t');
        ss >> tools;

        smatch dbcansubmatch;
        smatch diamondmatch;
        // 以下为你的原处理逻辑
        if (dbcan_sub.find("+") != string::npos && diamond.find("+") != string::npos) {
            cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
        }
        else if (dbcan_sub.find("+") != string::npos && diamond.find("+") == string::npos) {
            if (regex_search(diamond, diamondmatch, pattern_sub)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (regex_search(diamond, diamondmatch, pattern)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (diamond == "-") {
                cazyme_info.addCazyme(dbcan_sub + EC + ec);
            }
        }
        else if (dbcan_sub.find("+") == string::npos && diamond.find("+") != string::npos) {
            if (regex_search(dbcan_sub, dbcansubmatch, pattern_sub)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (regex_search(dbcan_sub, dbcansubmatch, pattern)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (dbcan_sub == "-") {
                cazyme_info.addCazyme(diamond + EC + ec);
            }
        }
        else {
            if (regex_search(dbcan_sub, dbcansubmatch, pattern_sub) && regex_search(diamond, diamondmatch, pattern_sub)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (regex_search(dbcan_sub, dbcansubmatch, pattern_sub) && !regex_search(diamond, diamondmatch, pattern_sub)) {
                if (regex_search(diamond, diamondmatch, pattern)) {
                    cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
                }
                else if (diamond == "-") {
                    cazyme_info.addCazyme(dbcan_sub + EC + ec);
                }
            }
            else if (!regex_search(dbcan_sub, dbcansubmatch, pattern_sub) && regex_search(diamond, diamondmatch, pattern_sub)) {
                if (regex_search(dbcan_sub, dbcansubmatch, pattern)) {
                    cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
                }
                else if (dbcan_sub == "-") {
                    cazyme_info.addCazyme(diamond + EC + ec);
                }
            }
            else if (regex_search(dbcan_sub, dbcansubmatch, pattern) && regex_search(diamond, diamondmatch, pattern)) {
                cazyme_info.addCazyme(dbcan_sub + " / " + diamond + EC + ec);
            }
            else if (regex_search(dbcan_sub, dbcansubmatch, pattern) && !regex_search(diamond, diamondmatch, pattern)) {
                cazyme_info.addCazyme(dbcan_sub + EC + ec);
            }
            else if (!regex_search(dbcan_sub, dbcansubmatch, pattern) && regex_search(diamond, diamondmatch, pattern)) {
                cazyme_info.addCazyme(diamond + EC + ec);
            }
            else if (!hmmer.empty() && hmmer.find("-") != string::npos) {
                smatch hmmermatch;
                if (regex_search(hmmer, hmmermatch, pattern)) {
                    cazyme_info.addCazyme("only hmmer: " + hmmermatch[0].str() + EC + ec);
                }
            }
        }
    }
    overview.close();
    std::cout << "File " << overview_file << " processed and closed" << endl;
}
// 用于第一层的行信息
struct lineforgrouping {
    vector<string> cazymain;
    vector<string> cazysub;
    vector<string> ec;
    string line; // 原行（含计数列）
};

// 第一层聚合
struct group_layer1 {
    vector<lineforgrouping> lines; // 该分组内的所有行
    string cazymain;               // 该分组共用的主类（或空字符串）
    string cazysub;                // 该分组共用的亚类（或空字符串）
};

// 第二层聚合
struct group_layer2 {
    vector<group_layer1> groups; // 包含多个 group_layer1
    string cazymain;             // 该大组的主类
};
// ============== 新增的函数 ==============
void doGrouping() {
    using namespace std;

    // 1) 读取 countDatamerged2.csv
    ifstream fin("countDatamerged2.csv");
    if (!fin.is_open()) {
        cerr << "Cannot open countDatamerged2.csv for grouping.\n";
        return;
    }

    // 读取表头
    string header;
    if (!std::getline(fin, header)) {
        cerr << "countDatamerged2.csv is empty?\n";
        return;
    }

    // 存放所有行(不含表头)
    vector<lineforgrouping> allLines;

    // 正则
    static const regex cazyMain_grouping_regex(R"([A-Z]+\d+)");
    static const regex cazySub_grouping_regex(R"([A-Z]+\d+_\d+)");

    // 逐行读取
    string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        lineforgrouping obj;
        obj.line = line;  // 保留整个行，以便后面输出

        // 拆出第一列：找 "\t" 第一次出现位置
        // 或者用更安全的方法：先 split
        vector<string> tokens;
        {
            stringstream ss(line);
            string tk;
            while (std::getline(ss, tk, '\t')) {
                tokens.push_back(tk);
            }
        }
        if (tokens.empty()) {
            // 该行可能有问题
            continue;
        }
        // tokens[0] => "CAZyme" + " EC:..."
        string cazyPart = tokens[0];

        // 解析 EC: ...
        // 先找到 "EC:"
        string ecPart;
        {
            size_t pos = cazyPart.find("EC:");
            if (pos != string::npos) {
                ecPart = cazyPart.substr(pos + 3); // 跳过 "EC:"
                cazyPart = cazyPart.substr(0, pos);
            }
        }
        // 解析 ecPart => 用 '|' 拆
        if (!ecPart.empty()) {
            // 去除前后空格
            while (!ecPart.empty() && (ecPart[0] == ' ' || ecPart[0] == ':')) {
                ecPart.erase(ecPart.begin());
            }
            stringstream ssec(ecPart);
            string e;
            while (std::getline(ssec, e, '|')) {
                // 去掉 : 后数字
                size_t cpos = e.find(':');
                if (cpos != string::npos) {
                    e = e.substr(0, cpos);
                }
                if (!e.empty()) {
                    obj.ec.push_back(e);
                }
            }
        }

        // 用正则匹配 cazyMain
        {
            auto begin = sregex_iterator(cazyPart.begin(), cazyPart.end(), cazyMain_grouping_regex);
            auto end = sregex_iterator();
            for (auto it = begin; it != end; ++it) {
                obj.cazymain.push_back(it->str());
            }
        }
        // 用正则匹配 cazySub
        {
            auto begin = sregex_iterator(cazyPart.begin(), cazyPart.end(), cazySub_grouping_regex);
            auto end = sregex_iterator();
            for (auto it = begin; it != end; ++it) {
                obj.cazysub.push_back(it->str());
            }
        }

        // 去重一下(如果有需要)
        // 例如 GH13_14+GH13_14 如果重复可能 push 两次
        sort(obj.cazymain.begin(), obj.cazymain.end());
        obj.cazymain.erase(unique(obj.cazymain.begin(), obj.cazymain.end()), obj.cazymain.end());

        sort(obj.cazysub.begin(), obj.cazysub.end());
        obj.cazysub.erase(unique(obj.cazysub.begin(), obj.cazysub.end()), obj.cazysub.end());

        sort(obj.ec.begin(), obj.ec.end());
        obj.ec.erase(unique(obj.ec.begin(), obj.ec.end()), obj.ec.end());

        allLines.push_back(obj);
    }
    fin.close();

    // 记录总行数(不含表头)
    int totalLines = (int)allLines.size();

    // 2) 第一层聚合 => vector<group_layer1>
    //    规则：如果行有 cazySub => 按 cazySub 分组；若无 => 按 cazyMain 分组
    //    但要注意，如果有多个 sub，可能你要拆分或者只取第一个。此处示例“假设每行只有 0 或 1 个 sub”。
    //    如确实有多 sub 的复杂场景，你需自行设计拆分或更多逻辑。
    map<string, group_layer1> grouping1_map;
    // key = cazySub (若存在)，否则 "#NO_SUB# + cazyMain" 或类似区分

    // 遍历 allLines
    for (auto& lf : allLines) {
        // 先看 lf.cazysub 是否为空
        if (!lf.cazysub.empty()) {
            // 假设只取第一个 cazySub[0] （如果你可能有多个 sub，要么按第一个，要么拆分多行等）
            string subName = lf.cazysub[0];
            // 同时 cazyMain 如果也有多个，这里也只能取第一个
            string mainName = lf.cazymain.empty() ? "" : lf.cazymain[0];

            // map key: subName
            auto it = grouping1_map.find(subName);
            if (it == grouping1_map.end()) {
                group_layer1 gl1;
                gl1.cazymain = mainName; // 可能是 GH13
                gl1.cazysub = subName;  // GH13_14
                gl1.lines.push_back(lf);
                grouping1_map[subName] = gl1;
            }
            else {
                // 已存在
                it->second.lines.push_back(lf);
            }
        }
        else {
            // 无 cazySub => 以 cazyMain 分组
            // 可能出现多 mainName；这里示例只取第一个
            // 若一个行里有 GH25,GH50 两个 main，你要自行决定如何聚合
            string mainName = lf.cazymain.empty() ? "#UNKNOWN#" : lf.cazymain[0];
            // key 
            string keyName = string("#NO_SUB#_") + mainName;
            auto it = grouping1_map.find(keyName);
            if (it == grouping1_map.end()) {
                group_layer1 gl1;
                gl1.cazymain = mainName;
                gl1.cazysub = "";  // no sub
                gl1.lines.push_back(lf);
                grouping1_map[keyName] = gl1;
            }
            else {
                it->second.lines.push_back(lf);
            }
        }
    }

    // 收集所有 group_layer1
    vector<group_layer1> all_group_layer1;
    all_group_layer1.reserve(grouping1_map.size());
    for (auto& pr : grouping1_map) {
        all_group_layer1.push_back(pr.second);
    }

    // 3) 第二层聚合 => map< string(cazymain), group_layer2 >
    //    group_layer2.cazymain = GH13 (例如)
    //    里面存放若干 group_layer1(可能sub=GH13_14, GH13_29等)
    map<string, group_layer2> grouping2_map;
    for (auto& g1 : all_group_layer1) {
        string m = g1.cazymain;
        if (m.empty()) m = "#EMPTY_MAIN#"; // 如果确实没有 main
        auto it = grouping2_map.find(m);
        if (it == grouping2_map.end()) {
            group_layer2 gl2;
            gl2.cazymain = m;
            gl2.groups.push_back(g1);
            grouping2_map[m] = gl2;
        }
        else {
            it->second.groups.push_back(g1);
        }
    }

    // 现在我们就拿到了多个 group_layer2，每个 group_layer2 包含若干 group_layer1

    // 4) 输出到文件夹
    //    - 创建 "groups" 文件夹
    {
        fs::path groupsDir("groups");
        if (!fs::exists(groupsDir)) {
            fs::create_directory(groupsDir);
        }
    }

    // 对每个 group_layer2 => 以它的 cazymain 命名文件夹
    // 对里面每个 group_layer1 => 以它的 cazysub 或 NO_SUB 命名 CSV
    // CSV 的第一行写 header，然后写 lines
    for (auto& kv : grouping2_map) {
        auto& gl2 = kv.second;
        // 创建文件夹 groups/<cazymain>
        // 注意: 如果 cazymain 有特殊字符，需替换或避免
        fs::path subDir = fs::path("groups") / fs::path(gl2.cazymain);
        // 若为空，就用 #EMPTY_MAIN# 作为文件夹名
        if (gl2.cazymain == "#EMPTY_MAIN#") {
            subDir = fs::path("groups") / fs::path("EMPTY_MAIN");
        }
        if (!fs::exists(subDir)) {
            fs::create_directory(subDir);
        }

        // 遍历 group_layer1
        for (auto& g1 : gl2.groups) {
            // 如果 g1.cazysub 非空 => GH13_29.csv
            // 否则 => GH13.csv 或 NO_SUB.csv
            string filename;
            if (!g1.cazysub.empty()) {
                filename = g1.cazysub + ".csv";
            }
            else {
                // 也可以用 e.g. "NO_SUB.csv" 或干脆用 gl2.cazymain+".csv"
                filename = gl2.cazymain + ".csv";
                if (gl2.cazymain == "#EMPTY_MAIN#") {
                    filename = "NO_SUB.csv";
                }
            }
            fs::path outPath = subDir / fs::path(filename);

            // 写 CSV
            ofstream ofs(outPath.string());
            if (!ofs.is_open()) {
                cerr << "Cannot write " << outPath << "\n";
                continue;
            }
            // 写表头
            ofs << header << "\n";
            // 写 lines
            for (auto& ln : g1.lines) {
                ofs << ln.line << "\n";
            }
            ofs.close();
        }
    }

    // 5) 统计信息 => group_summary.txt
    {
        ofstream gsum("group_summary.txt");
        if (!gsum.is_open()) {
            cerr << "Cannot write group_summary.txt\n";
            return;
        }
        // 统计“第一层分组”数量
        gsum << "Total lines (countDatamerged2.csv, excluding header): " << totalLines << "\n";
        gsum << "Number of group_layer1: " << all_group_layer1.size() << "\n";
        gsum << "Number of group_layer2: " << grouping2_map.size() << "\n\n";

        // 打印每个 group_layer2 -> 其下 group_layer1 数量
        int sum_g1_lines = 0;
        for (auto& kv : grouping2_map) {
            auto& gl2 = kv.second;
            gsum << "group_layer2( main=" << gl2.cazymain << " ) has " << gl2.groups.size() << " group_layer1\n";
            // 小结一下 group_layer1
            for (auto& g1 : gl2.groups) {
                gsum << "  -> group_layer1( sub=" << (g1.cazysub.empty() ? "NO_SUB" : g1.cazysub)
                    << " ) #lines=" << g1.lines.size() << "\n";
                sum_g1_lines += (int)g1.lines.size();
            }
            gsum << "\n";
        }
        // 检查所有行都被放入 group_layer1
        gsum << "Sum of all group_layer1.lines.size() = " << sum_g1_lines << "\n";
        if (sum_g1_lines == totalLines) {
            gsum << "All lines accounted for. OK.\n";
        }
        else {
            gsum << "Warning: mismatch?\n";
        }
        gsum.close();
    }
    // =============== 6) 额外：输出一个 grouped.csv ===============
   // 将所有行合并到一个 CSV 中，并且：
   // 1) 写表头
   // 2) 每个 group_layer2（按照 map 的遍历顺序） => 其下 group_layer1 => 其下多行
   // 3) group_layer1 的行彼此相邻
   // 4) group_layer2 与 group_layer2 之间的行也是相邻，但可以考虑是否插空行分隔

    {
        ofstream groupedOut("grouped.csv");
        if (!groupedOut.is_open()) {
            cerr << "Cannot write grouped.csv\n";
            return;
        }
        // 写表头
        groupedOut << header << "\n";

        // 遍历 grouping2_map
        for (auto& kv : grouping2_map) {
            auto& gl2 = kv.second;
            // 这里若想在每个 group_layer2 之间插空行，可这么做：
            // groupedOut << "\n"; // 视需求而定

            // 依次输出 group_layer1
            for (auto& g1 : gl2.groups) {
                // 同理，可在不同 group_layer1 之间插空行
                // groupedOut << "\n";

                // 输出该 group_layer1 包含的所有行
                for (auto& ln : g1.lines) {
                    groupedOut << ln.line << "\n";
                }
            }
        }

        groupedOut.close();
    }
    // =============== 7) 生成 annotation_row.txt ===============
   // 思路：
   //  (1) 要与 grouped.csv 的行顺序一致
   //  (2) 行名 (row.names) = 该行第一列字符串
   //  (3) 若 group_layer1.cazysub 不为空 => 标注用 g1.cazysub，否则用 g1.cazymain
   // 先收集 rowNamesVec, labelVec

    vector<string> rowNamesVec;
    vector<string> labelVec;

    // 小函数：从整行文本中提取“第一列”字符串
    auto getFirstColumn = [&](const string& fullLine) {
        // 找到第一个 '\t'
        // 或者把 fullLine split
        size_t tabPos = fullLine.find('\t');
        if (tabPos != string::npos) {
            return fullLine.substr(0, tabPos);
        }
        else {
            return fullLine; // 如果没有 '\t'，就返回整行
        }
        };

    // 同样循环一遍 grouping2_map => group_layer1 => lines
    // 与写 grouped.csv 一致
    for (auto& kv : grouping2_map) {
        auto& gl2 = kv.second;
        for (auto& g1 : gl2.groups) {
            // 优先用 g1.cazysub，如果为空则用 g1.cazymain
            string label = g1.cazysub.empty() ? g1.cazymain : g1.cazysub;

            for (auto& ln : g1.lines) {
                string rowName = getFirstColumn(ln.line);
                rowNamesVec.push_back(rowName);
                labelVec.push_back(label);
            }
        }
    }

    // 现在 rowNamesVec.size() = labelVec.size() = totalLines
    // 输出到 annotation_row.txt，形如：
    // annotation_row = data.frame(
    //   Label = c("GH13_29","GH13_29","GH1","GH1", ...)
    // )
    // row.names(annotation_row) = c("GH13_29 EC:-","GH13_29 EC:3.2.1.93", "GH1 EC:-", ...)
    // annotation_row

    {
        ofstream annoOut("annotation_row.txt");
        if (!annoOut.is_open()) {
            cerr << "Cannot write annotation_row.txt\n";
            return;
        }
        annoOut << "annotation_row = data.frame(\n";
        annoOut << "  Label = c(";

        // 写 labelVec
        for (size_t i = 0; i < labelVec.size(); i++) {
            // 需要加引号，逗号分隔
            annoOut << "\"" << labelVec[i] << "\"";
            if (i + 1 < labelVec.size()) annoOut << ",";
        }
        annoOut << ")\n)\n\n";

        // row.names(annotation_row) = c(...)
        annoOut << "row.names(annotation_row) = c(";
        for (size_t i = 0; i < rowNamesVec.size(); i++) {
            annoOut << "\"" << rowNamesVec[i] << "\"";
            if (i + 1 < rowNamesVec.size()) annoOut << ",";
        }
        annoOut << ")\n\n";

        annoOut << "annotation_row\n";
        annoOut.close();
    }
    cout << "Done grouping -> see ./groups folder, group_summary.txt, and grouped.csv"<<endl;
}
// ===================== 主函数 =====================
int main() {
    // ------- 1) 从两个目录收集数据 -> 输出原始 countDataEC.csv -------
    string main_folder1 = "F:\\large\\all_cris_out778";
    string main_folder2 = "F:\\large\\all_iners_out405";

    vector<CazymeInfo> all_cazyme_info;
    set<string> all_cazymes;

    // 遍历第一个目录
    {
        int count = 0;
        std::cout << "Beginning to iterate the folder all_cris_out" << endl;
        for (const auto& dir_entry : fs::directory_iterator(main_folder1)) {
            if (fs::is_directory(dir_entry)) {
                string subfolder = dir_entry.path().string();
                string overview_file, fna_file;
                for (const auto& file_entry : fs::directory_iterator(subfolder)) {
                    string file_name = file_entry.path().string();
                    if (file_name.ends_with("overview.txt")) {
                        overview_file = file_name;
                    }
                    else if (file_name.ends_with(".fna")) {
                        fna_file = file_name;
                    }
                }
                if (!overview_file.empty() && !fna_file.empty()) {
                    string strain_name = extractStrainName(fna_file);
                    if (!strain_name.empty()) {
                        CazymeInfo cazyme_info(strain_name);
                        processOverview(overview_file, cazyme_info);
                        all_cazyme_info.push_back(cazyme_info);
                        for (auto& kv : cazyme_info.cazyme_count) {
                            all_cazymes.insert(kv.first);
                        }
                    }
                }
            }
            count++;
            std::cout << "crispatus process: " << count << endl;
        }
    }

    // 遍历第二个目录
    {
        int count = 0;
        std::cout << "Beginning to iterate the folder all_iners_out" << endl;
        for (const auto& dir_entry : fs::directory_iterator(main_folder2)) {
            if (fs::is_directory(dir_entry)) {
                string subfolder = dir_entry.path().string();
                string overview_file, fna_file;
                for (const auto& file_entry : fs::directory_iterator(subfolder)) {
                    string file_name = file_entry.path().string();
                    if (file_name.ends_with("overview.txt")) {
                        overview_file = file_name;
                    }
                    else if (file_name.ends_with(".fna")) {
                        fna_file = file_name;
                    }
                }
                if (!overview_file.empty() && !fna_file.empty()) {
                    string strain_name = extractStrainName(fna_file);
                    if (!strain_name.empty()) {
                        CazymeInfo cazyme_info(strain_name);
                        processOverview(overview_file, cazyme_info);
                        all_cazyme_info.push_back(cazyme_info);
                        for (auto& kv : cazyme_info.cazyme_count) {
                            all_cazymes.insert(kv.first);
                        }
                    }
                }
            }
            count++;
            std::cout << "iners process: " << count << endl;
        }
    }

    // ------- 输出原始矩阵 countDataEC.csv -------
    {
        ofstream countDataFile("countDataEC.csv");
        if (!countDataFile.is_open()) {
            cerr << "Cannot write countDataEC.csv" << endl;
            return 1;
        }
        // 表头
        countDataFile << "CAZyme";
        for (auto& cinfo : all_cazyme_info) {
            countDataFile << "\t" << cinfo.strain_name;
        }
        countDataFile << "\n";

        // 行
        for (auto& cazy : all_cazymes) {
            // 排除以 e 开头
            if (regex_match(cazy, regex("e\\d+"))) {
                continue;
            }
            countDataFile << cazy;
            for (auto& cinfo : all_cazyme_info) {
                auto it = cinfo.cazyme_count.find(cazy);
                if (it != cinfo.cazyme_count.end()) {
                    countDataFile << "\t" << ((it->second>0)? 1:0);
                }
                else {
                    countDataFile << "\t0";
                }
            }
            countDataFile << "\n";
        }
        countDataFile.close();
    }

    // ------- 2) 第一次合并 -> 输出 countDatamerged.csv + 更新 summary.txt -------
    // 先读取 countDataEC.csv 进行第一次合并
    ifstream fin("countDataEC.csv");
    if (!fin.is_open()) {
        cerr << "Failed to open countDataEC.csv" << endl;
        return 1;
    }

    int originalRowCount = 0;
    int originalColCount = 0;

    // 读表头
    string headerLine;
    if (!std::getline(fin, headerLine)) {
        cerr << "countDataEC.csv is empty?\n";
        return 1;
    }
    originalRowCount++;
    vector<string> headerCols;
    {
        stringstream ss(headerLine);
        string col;
        while (std::getline(ss, col, '\t')) {
            headerCols.push_back(col);
        }
    }
    originalColCount = (int)headerCols.size();

    // 定义 ParsedLine (第一次合并需要的解析)
    struct ParsedLine {
        vector<int> counts;    // 对应各菌株计数(列)
        set<string> mainSet;
        set<string> subSet;
        set<string> eSet;
        set<string> ecSet;
    };

    // 一些正则
    const regex eCAMI_regex(R"([A-Z]+\d+_e\d+)");
    const regex cazySub_regex(R"([A-Z]+\d+_\d+)");
    const regex cazyMain_regex(R"([A-Z]+\d+)");

    auto stripEC = [&](const string& ec_in) {//可以把EC:1.10.3.2:77中的1.10.3.2提取出来
        size_t pos = ec_in.find(':');
        if (pos != string::npos) {
            return ec_in.substr(0, pos);
        }
        return ec_in;
        };

    // 用于 key 的序列化
    auto setToString = [&](const set<string>& s) {
        ostringstream oss;
        bool first = true;
        for (auto& elem : s) {
            if (!first) oss << ",";
            oss << elem;
            first = false;
        }
        return oss.str();
        };

    struct MergedGroup {
        set<string> mainSet;
        set<string> subSet;
        set<string> eSet;
        set<string> mergedEC;
        vector<int> mergedCounts;
    };

    map< tuple<string, string, string>, MergedGroup > mergedMap;//这一个map最后会保存合并后的结果，tuple中的三个string就是eCAMI, CAZysub, CAZymain

    // 开始读取数据行 (countDataEC.csv)
    {
        string line;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            originalRowCount++;

            vector<string> tokens;
            {
                stringstream ss(line);
                string t;
                while (std::getline(ss, t, '\t')) {//以'\t'为分隔把一行的每个元素分开，存到一个vector里
                    tokens.push_back(t);
                }
            }
            if ((int)tokens.size() != originalColCount) {
                cerr << "Line has inconsistent col count. skip\n";
                continue;
            }

            // 解析
            ParsedLine pl;
            pl.counts.resize(originalColCount - 1);

            // 从第2列开始是计数
            for (int i = 1; i < originalColCount; i++) {
                int val = 0;
                try { val = stoi(tokens[i]); }
                catch (...) {}
                pl.counts[i - 1] = val;//这里是为什么pl.counts.resize()时传入了originalColCount-1，不用储存第一列
            }

            // 拆分第一列 -> cazyPart & ecPart
            string cazyPart = tokens[0];
            string ecPart;
            {
                size_t pos = cazyPart.find("EC:");
                if (pos != string::npos) {
                    ecPart = cazyPart.substr(pos + 3);//EC: 之后
                    cazyPart = cazyPart.substr(0, pos);//EC: 之前
                }
            }

            // 解析 EC
            if (!ecPart.empty()) {
                // 去除空格/冒号
                while (!ecPart.empty() && (ecPart[0] == ' ' || ecPart[0] == ':')) {
                    ecPart.erase(ecPart.begin());
                }
                // 用 '|'
                stringstream ssec(ecPart);
                string tok;
                while (std::getline(ssec, tok, '|')) {
                    tok = stripEC(tok);//以'|'分隔的EC号是3.2.1.41:25|3.2.1.1:2|3.2.1.68:2 的形式，要去掉:后的数字
                    if (!tok.empty()) {
                        pl.ecSet.insert(tok);
                    }
                }
            }

            // 找 eCAMI
            {
                auto begin = sregex_iterator(cazyPart.begin(), cazyPart.end(), eCAMI_regex);
                auto end = sregex_iterator();
                for (auto it = begin; it != end; ++it) {
                    pl.eSet.insert(it->str());//从cazyPart中找出所有能匹配eCAMI_regex的部分并全部insert到pl.sSet中
                }
            }
            // 找 sub
            {
                auto begin = sregex_iterator(cazyPart.begin(), cazyPart.end(), cazySub_regex);
                auto end = sregex_iterator();
                for (auto it = begin; it != end; ++it) {
                    pl.subSet.insert(it->str());
                }
            }
            // 找 main
            {
                auto begin = sregex_iterator(cazyPart.begin(), cazyPart.end(), cazyMain_regex);
                auto end = sregex_iterator();
                for (auto it = begin; it != end; ++it) {
                    string m = it->str();
                    // 如果 m 已出现在 eSet/subSet，则跳过
                    if (pl.eSet.count(m)) continue;
                    if (pl.subSet.count(m)) continue;//这里不要改，涉及到后面输出
                    pl.mainSet.insert(m);
                }
            }

            // 放进 mergedMap
            string mk = setToString(pl.mainSet);
            string sk = setToString(pl.subSet);
            string ek = setToString(pl.eSet);//以这三个set组成tuple，然后作为mergedMap的key
            auto key = make_tuple(mk, sk, ek);
            //以下涉及到合并操作
            auto itM = mergedMap.find(key);//mergedMap本来是空的，从mergedMap
            if (itM == mergedMap.end()) {//如果没有在mergeMap中没找到，建一个新的MergedGroup
                MergedGroup mg;
                mg.mainSet = pl.mainSet;
                mg.subSet = pl.subSet;
                mg.eSet = pl.eSet;
                mg.mergedEC = pl.ecSet;
                mg.mergedCounts = pl.counts;
                mergedMap[key] = mg;
            }
            else {//如果找到了，则合并
                MergedGroup& mg = itM->second;
                //对所有counts进行or操作 counts: or
                for (size_t i = 0; i < mg.mergedCounts.size(); i++) {
                    if (mg.mergedCounts[i] > 0 || pl.counts[i] > 0) {
                        mg.mergedCounts[i] = 1;
                    }
                    else {
                        mg.mergedCounts[i] = 0;
                    }
                }
                // ec
                for (auto& e : pl.ecSet) {
                    mg.mergedEC.insert(e);//这里的insert操作实际会检查是否重合，重合的不会加进去
                }
            }
        }//while (std::getline(fin, line))结束，每行都进行了如上操作
    }
    fin.close();

    // 输出第一次合并结果到 countDatamerged.csv
    ofstream fout("countDatamerged.csv");
    if (!fout.is_open()) {
        cerr << "Cannot open countDatamerged.csv\n";
        return 1;
    }

    // 写表头
    fout << headerLine << "\n";

    // mergedMap.size() = 第一次合并后去重的行数（除表头）
    int firstPassMergedRowCount = (int)mergedMap.size();
    int firstPassMergedColCount = originalColCount; // 与原来相同

    // 写数据，注意这里第一列的输出实际上不太好，但是在第二次合并被改善
    for (auto& kv : mergedMap) {
        const auto& mg = kv.second;
        // 拼第一列
        // 这里仍然保留 eCAMI 输出(括号)，如你当前的逻辑
        // 只是演示： main+sub + ( eSet ) + " EC:..."
        // -----------
        // 先准备 mainVec/subVec/eVec
        vector<string> mainVec(mg.mainSet.begin(), mg.mainSet.end());//把set转为vector
        vector<string> subVec(mg.subSet.begin(), mg.subSet.end());
        vector<string> eVec(mg.eSet.begin(), mg.eSet.end());

        // 合并 mainVec+subVec => combined
        vector<string> combined;
        combined.insert(combined.end(), mainVec.begin(), mainVec.end());
        combined.insert(combined.end(), subVec.begin(), subVec.end());

        // 用 '+' 拼起来
        ostringstream firstCol;
        bool firstItem = true;
        for (auto& c : combined) {
            if (!firstItem) firstCol << "+";
            firstCol << c;
            firstItem = false;
        }
        // 如果有 eVec，则用括号逗号拼
        if (!eVec.empty()) {
            firstCol << "(";
            for (size_t i = 0; i < eVec.size(); i++) {
                if (i > 0) firstCol << ",";
                firstCol << eVec[i];
            }
            firstCol << ")";
        }
        // 如果有 EC
        if (!mg.mergedEC.empty()) {
            firstCol << " EC:";
            bool firstEC = true;
            for (auto& ec : mg.mergedEC) {
                if (!firstEC) firstCol << "|";
                firstCol << ec;
                firstEC = false;
            }
        }

        // 输出第一列，这里第一列的输出逻辑实际上不好，在第二次合并被改善
        fout << firstCol.str();

        // 输出计数列
        for (auto val : mg.mergedCounts) {
            fout << "\t" << val;
        }
        fout << "\n";
    }
    fout.close();

    // ------ 输出第一次合并后的 summary （原始表 + 一次合并后） ------
    ofstream summaryFile("summary.txt");
    if (!summaryFile.is_open()) {
        cerr << "Cannot open summary.txt\n";
        return 1;
    }
    summaryFile << "Original row count: " << originalRowCount << "\n";
    summaryFile << "Original col count: " << originalColCount << "\n";
    summaryFile << "Merged row count (1st pass): " << (1 + firstPassMergedRowCount) << "\n";
    summaryFile << "Merged col count (1st pass): " << firstPassMergedColCount << "\n";
    // （此时尚未进行二次合并，所以先把 2nd pass 留在后面）
    summaryFile.close();

    // ===================== 3) 第二次合并 =====================
    // 需求：若两行或多行第一列的 CAZymain 相同，CAZysub 相同，EC 相同，则可以合并，忽略 eCAMI；计数做位或。
    // 并且如果 sub 存在 GH13_14，就不再单独输出 main 中的 GH13；省得重复。

    ifstream fin2("countDatamerged.csv");
    if (!fin2.is_open()) {
        cerr << "Cannot open countDatamerged.csv for reading second pass.\n";
        return 1;
    }

    // 读取表头
    string headerLine2;
    if (!std::getline(fin2, headerLine2)) {
        cerr << "countDatamerged.csv empty?\n";
        return 1;
    }

    // 第二次合并的列数
    vector<string> headerCols2;
    {
        stringstream ss(headerLine2);
        string col;
        while (std::getline(ss, col, '\t')) {
            headerCols2.push_back(col);
        }
    }
    int secondPassColCount = (int)headerCols2.size();

    // 为第二次合并定义一个简单的结构
    struct ParsedLine2 {
        vector<int> counts;         // 从第2列起
        set<string> mainSet;        // 不再包含 eCAMI
        set<string> subSet;
        set<string> ecSet;          // 解析到的 EC
    };

    // 我们的规则：
    //  - 忽略 eCAMI（即 [A-Z]+\\d+_e\\d+）
    //  - 提取 sub: [A-Z]+\\d+_\\d+
    //  - 提取 main: [A-Z]+\\d+
    //  - 解析 EC
    //  - 如果 sub 中出现 GH13_14，就把 main 里的 GH13 删掉；(例如 sub="GH13_14" => main 去除 "GH13")
    //  - 行合并条件： mainSet 相同 & subSet 相同 & ecSet 相同 => 计数做 OR

    const regex sub_regex_2(R"([A-Z]+\d+_\d+)");
    const regex main_regex_2(R"([A-Z]+\d+)");

    // 用来做 key
    auto setToString2 = [&](const set<string>& s) {
        // 还是和之前类似
        ostringstream oss;
        bool first = true;
        for (auto& elem : s) {
            if (!first) oss << ",";
            oss << elem;
            first = false;
        }
        return oss.str();
        };

    struct MergedGroup2 {
        set<string> mainSet;
        set<string> subSet;
        set<string> ecSet;
        vector<int> counts;
    };

    map< tuple<string, string, string>, MergedGroup2 > mergedMap2;

    int secondPassOriginalRowCount = 0; // 二次合并前(即 countDatamerged.csv 的总行数)
    // 读数据
    {
        string line;
        while (std::getline(fin2, line)) {
            if (line.empty()) continue;
            secondPassOriginalRowCount++;

            // 拆分
            vector<string> tokens;
            {
                stringstream ss(line);
                string t;
                while (std::getline(ss, t, '\t')) {
                    tokens.push_back(t);
                }
            }
            if ((int)tokens.size() != secondPassColCount) {
                cerr << "In second pass, line col != header col. skip.\n";
                continue;
            }

            // 第 1 个 token => "CAZyme" + " EC:..."
            // 剩下的是计数
            ParsedLine2 pl2;
            pl2.counts.resize(secondPassColCount - 1, 0);

            for (int i = 1; i < secondPassColCount; i++) {
                int val = 0;
                try { val = stoi(tokens[i]); }
                catch (...) {}
                pl2.counts[i - 1] = val;
            }

            // 拆出 cazyPart & ecPart
            string cazyPart2 = tokens[0];
            string ecPart2;
            {
                size_t pos = cazyPart2.find("EC:");
                if (pos != string::npos) {
                    ecPart2 = cazyPart2.substr(pos + 3);
                    cazyPart2 = cazyPart2.substr(0, pos);
                }
            }
            // 解析 ecPart2
            if (!ecPart2.empty()) {
                while (!ecPart2.empty() && (ecPart2[0] == ' ' || ecPart2[0] == ':')) {
                    ecPart2.erase(ecPart2.begin());
                }
                stringstream ssec(ecPart2);
                string ecTok;
                while (std::getline(ssec, ecTok, '|')) {
                    // 去除 : 后数字，这里可能已经不需要了
                    size_t cpos = ecTok.find(':');
                    if (cpos != string::npos) {
                        ecTok = ecTok.substr(0, cpos);
                    }
                    if (!ecTok.empty()) {
                        pl2.ecSet.insert(ecTok);
                    }
                }
            }
            //以下仍然从每行第一列找出所有的可匹配项
            // 1) 忽略 eCAMI => 不做任何匹配 [A-Z]+\\d+_e\\d+
            // 2) 匹配 sub [A-Z]+\\d+_\\d+
            {
                auto begin = sregex_iterator(cazyPart2.begin(), cazyPart2.end(), sub_regex_2);
                auto end = sregex_iterator();
                for (auto it = begin; it != end; ++it) {
                    pl2.subSet.insert(it->str());
                }
            }
            // 3) 匹配 main [A-Z]+\\d+ 
            {
                auto begin = sregex_iterator(cazyPart2.begin(), cazyPart2.end(), main_regex_2);
                auto end = sregex_iterator();
                for (auto it = begin; it != end; ++it) {
                    string m = it->str();
                    // 可能是 GH13, 也可能已经被 sub 包含(如 GH13_14)
                    // 先加到 mainSet, 待会再判断移除
                    pl2.mainSet.insert(m);
                }
            }

            // 如果 subSet 中含有 GH13_14，就去除 mainSet 中的 GH13
            // 例如 sub="GH13_14" => sub.substr(0, sub.find('_')) => "GH13", 则在 mainSet.erase("GH13")
            // 做法：对每个 sub, 如果有 '_', 取下划线前面的字符串 => cand
            // 然后 if(mainSet.count(cand)) => mainSet.erase(cand)
            for (auto& subItem : pl2.subSet) {
                // 找到下划线
                size_t pos = subItem.find('_');
                if (pos != string::npos) {
                    string cand = subItem.substr(0, pos); // GH13
                    // 如果 mainSet 有这项，则删除
                    if (pl2.mainSet.count(cand) > 0) {
                        pl2.mainSet.erase(cand);
                    }
                }
            }

            // 计算 key: (mainSet, subSet, ecSet)
            string mk2 = setToString2(pl2.mainSet);
            string sk2 = setToString2(pl2.subSet);
            string ek2 = setToString2(pl2.ecSet);
            auto key2 = make_tuple(mk2, sk2, ek2);

            // 合并
            auto itM2 = mergedMap2.find(key2);
            if (itM2 == mergedMap2.end()) {
                // 新建
                MergedGroup2 mg2;
                mg2.mainSet = pl2.mainSet;
                mg2.subSet = pl2.subSet;
                mg2.ecSet = pl2.ecSet;
                mg2.counts = pl2.counts;
                mergedMap2[key2] = mg2;
            }
            else {
                // OR
                auto& mg2 = itM2->second;
                for (size_t i = 0; i < mg2.counts.size(); i++) {
                    if (mg2.counts[i] > 0 || pl2.counts[i] > 0) {
                        mg2.counts[i] = 1;
                    }
                    else {
                        mg2.counts[i] = 0;
                    }
                }
            }
        }
    }
    fin2.close();

    // 现在 mergedMap2.size() 就是二次合并后的行数(不含表头)
    int secondPassMergedRowCount = (int)mergedMap2.size();

    // 写二次合并结果到 countDatamerged2.csv
    ofstream fout2("countDatamerged2.csv");
    if (!fout2.is_open()) {
        cerr << "Cannot open countDatamerged2.csv for writing.\n";
        return 1;
    }
    // 写表头(与 countDatamerged.csv 一样)
    fout2 << headerLine2 << "\n";

    // 输出
    for (auto& kv : mergedMap2) {
        const auto& mg2 = kv.second;
        // 拼第一列 => main+sub
        // 不输出 eCAMI (因为我们已经忽略了它)
        // 用 '+' 连接 mainVec + subVec
        // 若有 EC => " EC:..."
        vector<string> mainVec(mg2.mainSet.begin(), mg2.mainSet.end());
        vector<string> subVec(mg2.subSet.begin(), mg2.subSet.end());

        // 拼接 mainVec+subVec
        vector<string> combined;
        combined.insert(combined.end(), mainVec.begin(), mainVec.end());
        combined.insert(combined.end(), subVec.begin(), subVec.end());

        ostringstream firstCol;
        bool firstItem = true;
        for (auto& c : combined) {
            if (!firstItem) firstCol << "+";
            firstCol << c;
            firstItem = false;
        }

        if (!mg2.ecSet.empty()) {
            firstCol << " EC:";
            bool firstEC = true;
            for (auto& ec : mg2.ecSet) {
                if (!firstEC) firstCol << "|";
                firstCol << ec;
                firstEC = false;
            }
        }

        // 写第一列
        fout2 << firstCol.str();
        // 写 counts
        for (auto val : mg2.counts) {
            fout2 << "\t" << val;
        }
        fout2 << "\n";
    }
    fout2.close();

    // 最后，更新 summary.txt，追加二次合并的行列数
    // 这里为了简单，直接以追加模式打开
    ofstream summaryFile2("summary.txt", ios::app);
    if (!summaryFile2.is_open()) {
        cerr << "Cannot append summary.txt\n";
        return 1;
    }
    // secondPassOriginalRowCount 包含表头, secondPassMergedRowCount 不含表头
    summaryFile2 << "\n--- 2nd pass merging info ---\n";
    summaryFile2 << "Merged row count (2nd pass, final): " << (1 + secondPassMergedRowCount) << "\n";
    summaryFile2 << "Merged col count (2nd pass, final): " << secondPassColCount << "\n";
    summaryFile2.close();
    
    // secondPassMergedRowCount 写入 summary.txt 等操作都结束后：
    // ------------------------------------------
    // 在这里新增调用:
    doGrouping();
    // 打印提示
    std::cout << "strain names not found: " << notfound << endl;
    for (auto& nf : notfoundstrains) {
        cout << nf << endl;
    }
    cout << "First pass => countDatamerged.csv generated.\n";
    cout << "Second pass => countDatamerged2.csv generated.\n";
    cout << "summary.txt updated.\n";
    cout << "Program finished.\n";
    return 0;
}
