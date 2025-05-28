#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <set>
#include <sstream>
#include <filesystem>
#include <map>
#include <algorithm>
#include <stack>
#include <unordered_set>
#include <array>
#include <opencv4/opencv2/opencv.hpp>

namespace fs = std::filesystem;
using std::set, std::vector, std::string, std::map, std::ifstream,std::endl,std::string;

const string KEGG_Pathway_HTML_PNG = "D:\\kegg_pathway";

std::ofstream g_logFile("Logs.txt");
std::ofstream g_errorFile("Errors.txt");

//KEGG通路图中有两种数据类型，圆形代表物质，矩形代表基因（酶）
enum class Shape {DEFAULT, RECT, CIRCLE };

//这个基类包含了圆形或矩形都应有的属性，id和形状分类
class ElementBase {
public:
    std::string area_id;
    Shape shape;
    ElementBase():area_id(""),shape(Shape::DEFAULT){}
    virtual void parseCoords(const std::string& coords) = 0;
    virtual void parseTitle(const std::string& title) = 0;
    virtual ~ElementBase() = default;
};
//矩形有两个坐标，分别确定左上点和右下点，还有ko与基因名的映射表，EC号，反应号
class RectElement : public ElementBase {
public:
    std::pair<int, int> point1;  // 左上角
    std::pair<int, int> point2;  // 右下角
    std::map<std::string, std::string> ko_to_geneName;
    std::string EC;
    std::string reaction;

    //parseCoords函数假设坐标只有四个数字
    void parseCoords(const std::string& coords) override {
        std::stringstream ss(coords);
        std::string token;
        getline(ss, token, ','); point1.first = std::stoi(token)*2-3;
        getline(ss, token, ','); point1.second = std::stoi(token)*2-3;
        getline(ss, token, ','); point2.first = std::stoi(token)*2-3;
        getline(ss, token); point2.second = std::stoi(token)*2-3;
    }
    //这里parseTitle函数也假设title有一个特定的格式，进行此函数前应检查title
    void parseTitle(const std::string& title) override {
        std::stringstream ss(title);
        std::string item;
        while (getline(ss, item, ',')) {//html源码中这一部分由","分隔
            item.erase(0, item.find_first_not_of(" \t"));
            item.erase(item.find_last_not_of(" \t") + 1);//清空两边的空格

            if (item.empty()) continue;

            size_t parenOpen = item.find('(');
            size_t parenClose = item.find(')');
            if (parenOpen != std::string::npos && parenClose != std::string::npos) {
                std::string ko = item.substr(0, parenOpen - 1);//ko号："("之前的部分
                std::string gene = item.substr(parenOpen + 1, parenClose - parenOpen - 1);//括号里的部分，不一定是基因名
                ko_to_geneName[ko] = gene;
            }
            else if (item.find('.') != std::string::npos) {
                EC=item;//如果有"."，说明是EC号
            }
            else if (item[0] == 'R') {
                reaction=item;//有“R”说明是反应号
            }
        }
    }
};
//一个圆形有一个圆心坐标和半径，以及物质的id与名称
class CircleElement : public ElementBase {
public:
    std::pair<int, int> center;
    int radius;
    std::string compound_id;
    std::string compound_name;

    void parseCoords(const std::string& coords) override {
        std::stringstream ss(coords);
        std::string token;
        getline(ss, token, ','); center.first = std::stoi(token);
        getline(ss, token, ','); center.second = std::stoi(token);
        getline(ss, token); radius = std::stoi(token);
    }

    void parseTitle(const std::string& title) override {
        size_t parenOpen = title.find('(');
        size_t parenClose = title.find(')');
        if (parenOpen != std::string::npos && parenClose != std::string::npos) {
            compound_id = title.substr(0, parenOpen - 1);
            compound_name = title.substr(parenOpen + 1, parenClose - parenOpen - 1);
        }
    }
};

// 这个提取菌株名的函数从.fna文件的首行提取菌株名，但
// 进行多菌种比对时，下面这个extractStrainName就需要更改了，以GCA编号从json文件中提取菌株名
std::string extractStrainName(const std::string& fna_file) {
    static int notfound = 0;
    static std::vector<std::string> notfoundstrains;

    std::ifstream fna(fna_file);
    if (!fna.is_open()) {
        std::cerr << "Could not open .fna file: " << fna_file << std::endl;
        g_errorFile << "Could not open .fna file: " << fna_file << std::endl;
        notfound++;
        notfoundstrains.push_back(fna_file);
        return "";
    }

    static std::vector<std::string> known_strains = {
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

    std::string line, strain_name;
    std::regex strain_pattern1(R"(Lactobacillus\s+[A-Za-z]+\s+\S+)");
    std::regex strain_pattern2(R"(Lactobacillus\s+[A-Za-z]+\s+strain\s+\S+)");
    std::regex strain_pattern3(R"(Lactobacillus\s+[A-Za-z]+\s+isolate\s+\S+)");

    if (std::getline(fna, line)) {
        std::smatch match;

        // 先检查是否在已知菌株名中
        for (const auto& strain : known_strains) {
            if (line.find(strain) != std::string::npos) {
                g_logFile << "Known strain matched: " << strain << std::endl;
                strain_name = strain;
                break;
            }
        }
        // 如果没找到已知菌株名，则用正则
        if (strain_name.empty()) {
            if (std::regex_search(line, match, strain_pattern1)) {
                strain_name = match[0];
                g_logFile << "strain_pattern1 matched: " << strain_name << std::endl;
            }
            if ((strain_name.find("Lactobacillus crispatus strain") != std::string::npos ||
                strain_name.find("Lactobacillus iners strain") != std::string::npos) &&
                std::regex_search(line, match, strain_pattern2))
            {
                strain_name = match[0];
                g_logFile << "strain_pattern2 matched: " << strain_name << std::endl;
            }
            if ((strain_name.find("Lactobacillus crispatus isolate") != std::string::npos ||
                strain_name.find("Lactobacillus iners isolate") != std::string::npos) &&
                std::regex_search(line, match, strain_pattern3))
            {
                strain_name = match[0];
                g_logFile << "strain_pattern3 matched: " << strain_name << std::endl;
            }
        }
    }

    if (strain_name.empty()) {
        g_logFile << "Strain_name in fna file: " << fna_file << " not found" << std::endl;
        g_errorFile << "Strain_name in fna file: " << fna_file << " not found" << std::endl;
        notfound++;
        notfoundstrains.push_back(fna_file);
    }

    fna.close();
    return strain_name;
}
//Gene类定义不需要改变，表示一个基因，eggNOG注释结果中一个基因有如下信息
class Gene {
public:
    std::string gene_id;
    std::string seed_ortholog;
    double evalue;
    double score;
    std::string eggnog_ogs;
    std::string max_annot_lvl;
    std::string cog_category;
    std::string description;
    std::string preferred_name;

    std::vector<std::string> GOs;
    std::vector<std::string> EC;
    std::vector<std::string> kegg_ko;
    std::vector<std::string> kegg_pathway;
    std::vector<std::string> kegg_module;
    std::vector<std::string> kegg_reaction;
    std::vector<std::string> kegg_rclass;
    std::vector<std::string> brite;
    std::vector<std::string> kegg_tc;
    std::vector<std::string> cazy;
    std::vector<std::string> bigg_reaction;
    std::vector<std::string> pfams;

    std::string origin_strain_name; // 来源菌株名称

    // 构造函数
    Gene() {}

    // shrink_to_fit 和其他函数不变
    void shrinkAllVectors() {
        GOs.shrink_to_fit();
        EC.shrink_to_fit();
        kegg_ko.shrink_to_fit();
        kegg_pathway.shrink_to_fit();
        kegg_module.shrink_to_fit();
        kegg_reaction.shrink_to_fit();
        kegg_rclass.shrink_to_fit();
        brite.shrink_to_fit();
        kegg_tc.shrink_to_fit();
        cazy.shrink_to_fit();
        bigg_reaction.shrink_to_fit();
        pfams.shrink_to_fit();
    }
};

// 从根文件夹中提取菌种名
std::vector<std::string> getSpeciesFromFolder(const std::string& rootFolderPath) {
    std::vector<std::string> species;

    // 遍历总文件夹中的所有子文件夹
    for (const auto& entry : fs::directory_iterator(rootFolderPath)) {
        if (fs::is_directory(entry)) {
            species.push_back(entry.path().filename().string());  // 将子文件夹的名字作为菌种名称加入 species 向量
        }
    }
    g_logFile << "Species: ";
    std::cout << "Species: ";
    for (const auto& s : species) {
        g_logFile << s << "  ";
        std::cout << s << "  ";
    }
    std::cout << std::endl;
    g_logFile << std::endl;
    return species;
}

//一个菌株的基因集
class Strain {
public:
    string species;
    std::vector<Gene> geneset;
    std::string strain_name;  // 从 .fna 文件中解析到的菌株名
};

//一个菌种的所有菌株
class Species_Strain_set {
public:
    string species;               // 这个 Strain_set 所属菌种
    std::vector<Strain> all_strains; // 一个菌种的所有菌株

    // 如果使用GeneMap（std::map<std::string, std::vector<Gene>>）比较占用内存
    // 考虑改为map<string,vector<Gene*>>，以指针指向all_strains中Strain的geneset包含的Gene对象
    // 下面这个gname_vec代表这一个菌种内每个基因名到Gene指针的映射
    // 相当于从每个基因名向所有菌株的基因的映射
    map<string, vector<Gene*>> gname_vec;

    //下面这个clusters是在gname_vec中的vector<Gene*>进行基因间对比后，把各类信息完全相同的Gene*给聚在一起了
    map<string, std::vector<std::vector<Gene*>>> clusters;

    // KEGG编号-->(基因名-->基因cluster)，比较复杂
    map<string, map<string, vector<vector<Gene*>>>> KEGGtoGeneclusters;

    // 上面这些map有基因名到所有基因，基因名到聚类后的基因，还有KEGG到GeneMap

    // 现在需要EC号到菌株的映射，KO号到菌株，反应号到菌株
    // 这些结构是为了之后在KEGG图中标记时，确定某个EC号，KO号，反应号在哪些菌种中出现
    map<string, set<string>> ECtoStrains;
    map<string, set<string>> KOtoStrains;
    map<string, set<string>> ReactiontoStrains;

    std::set<string> geneNames; // 保存这些菌株出现过的所有基因名
};

void ModifyECKOReac(vector<Species_Strain_set>& SSsvec) {
    for (Species_Strain_set& SSs : SSsvec) {
        // 从clusters中得到信息
        for(const auto& [gname,gptrvec]:SSs.gname_vec) {
            // 遍历这个基因名对应的基因指针vector
            for (const auto& gptr : gptrvec) {
                //每个基因的EC也是个vector
                for (const auto& EC : gptr->EC) {
                    SSs.ECtoStrains[EC].insert(gptr->origin_strain_name);
                }
                for (const auto& KO : gptr->kegg_ko) {
                    SSs.KOtoStrains[KO].insert(gptr->origin_strain_name);
                }
                for (const auto& Reac : gptr->kegg_reaction) {
                    SSs.ReactiontoStrains[Reac].insert(gptr->origin_strain_name);
                }
            }
        }
    }
}

// 高度注意这个函数，很有可能问题出在这里
void modifyKEGGtoGeneclusters(Species_Strain_set& strainset) {
    // 假设此时strainset的clusters已经被填充
    // 先把strainset中所有的KEGGPathway编号都找出来
    set<string> all_KEGG_Pathway_code; // set有自动去重功能
    for (const auto& [_, clusters] : strainset.clusters) {
        for (const auto& onecluster : clusters) {
            for (const auto& keggPathway : onecluster[0]->kegg_pathway) {
                // 检查 keggPathway 中是否包含 "map"
                if (keggPathway.find("map") != std::string::npos) {
                    // 如果包含 "map"，则插入到 all_KEGG_Pathway_code
                    all_KEGG_Pathway_code.insert(keggPathway);
                }
            }
        }
    }
    
    g_logFile << "All KEGGPathway codes: " << all_KEGG_Pathway_code.size() << std::endl;
    for (const auto& onecode : all_KEGG_Pathway_code) {
        g_logFile << onecode << "  ";
    }
    g_logFile << std::endl;

    for (const auto& onecode : all_KEGG_Pathway_code) {
        for ( auto& [gname, clusters] : strainset.clusters) {
            for ( auto& c : clusters) {
                //如果在一个cluster的KEGG_Pathway编号中发现了onecode，当前的KEGG_Pathway编号
                if (std::find(c[0]->kegg_pathway.begin(), c[0]->kegg_pathway.end(),onecode) != c[0]->kegg_pathway.end()) {
                    strainset.KEGGtoGeneclusters[onecode][gname].push_back(c); //
                }
            }
            /*
            g_logFile << "strainset: " << strainset.species << "KEGGtoGeneclusters[" << onecode << "][" << gname << "]: " << strainset.KEGGtoGeneclusters[onecode][gname].size() <<"  sizes of one cluster: ";
            for (const auto& Onecluster : strainset.KEGGtoGeneclusters[onecode][gname]) {
                g_logFile << Onecluster.size() << "  ";
            }
            g_logFile << std::endl;
            */
        }
    }
}

//这个类用于KeggPathwayGeneSet的单个基因表示
class GeneforKPGS {
public:
    string genename;

    // 菌种名-->菌株数量
    map<string, int> species_strain_num;

    set<string> KOs;
    set<string> ECs;
    set<string> reactions;
};

// 这个类是为了生成PathwayGeneSet.txt文件服务的，一个pathwaygeneset要有一个KEGG通路编号，一些基因
// 以及这些基因在各菌种的菌株数
class KeggPathwayGeneSet {
public:
    std::string kegg_pathway;  // KEGG 路径编号

    vector<GeneforKPGS> genes;

    // 现在可以写一个函数，遍历每个菌种的clusters，根据cluster所在的KEGG_Pathway为它们归类

    // 构造函数
    KeggPathwayGeneSet() : kegg_pathway("") {}

    KeggPathwayGeneSet(const std::string& pathway) : kegg_pathway(pathway) {}
};

void GenerateKeggPathwaySet(vector<Species_Strain_set>& species_strain_sets,
    set<string>& all_KEGG_Pathway_code,
    vector<KeggPathwayGeneSet>& KPGSvec) {

    //onespecies代表一个Species_Strain_set
    for (const auto& onespecies : species_strain_sets) {
        //每个菌种的gname_vec和clusters是需要的，遍历clusters
        for (const auto& [gname, gclusters] : onespecies.clusters) {
            //下面的indexes表示每个cluster
            for (const auto& onecluster : gclusters) {
                //每个cluster中各个Gene都是相同的，因此随便取一个即可，取第一个作为索引
                for (const auto& keggPathway : onecluster[0]->kegg_pathway) {
                    // 检查 keggPathway 中是否包含 "map"
                    if (keggPathway.find("map") != std::string::npos) {
                        // 如果包含 "map"，则插入到 all_KEGG_Pathway_code
                        all_KEGG_Pathway_code.insert(keggPathway);
                    }
                }
            }
        }
    }

    // 现在得到了 all_KEGG_Pathway_code，可以输出 pathwaygeneset.txt
    std::ofstream outputFile("KEGG_Pathway_GeneSet.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file for writing." << std::endl;
        return;
    }

    // 需要遍历all_KEGG_Pathway_code
    for (const string& one_pathway : all_KEGG_Pathway_code) {
        outputFile << "Pathway: " << one_pathway << std::endl;
        KeggPathwayGeneSet kpgs(one_pathway);
        //然后填充Species_to_GeneMap，首先要得知菌种信息
        for (auto& onespecies : species_strain_sets) {
            // Species_Strain_set 类有成员 KEGGtoGeneclusters，可以充分利用这个成员，因为当前在一个pathway中，不需要遍历Species_Strain_set的整个gname_vec
            // 首先以one_pathway索引到KEGGtoGeneclusters的相应map<string,vector<vector<Gene*>>>，这个map保存了从基因名到Gene*的cluster
            // 然后可以遍历这个map<string,vector<vector<Gene*>>>，以每个基因名（每个键值对的键）建立新的GeneforKPGS，或者在kpgs的genes中查找到这个基因名对应的GeneforKPGS
            // 遍历map<string,vector<vector<Gene*>>>时，还要遍历vector<vector<Gene*>>里的每个Gene*
            // 把每个Gene*的origin_strain_name都insert到set中，然后以这个set的size()作为相应GeneforKPGS的species_strain_num的菌株数量，注意赋值给species_strain_num[onespecies.species]
            // 还需要给GeneforKPGS的KOs, ECs, Reactions赋值，可以遍历当前vector<vector<Gene*>>的每个vector<Gene*>的第一个Gene*，用三个set<string>分别insert所有vector<Gene*>的第一个Gene*的kegg_ko, EC, kegg_reaction
            // 然后把这三个set分别赋给GeneforKPGS的KOs, ECs, Reactions
            // 在species_strain_sets遍历完成时，kpgs也填充完毕了，然后就可以在outputFile中输出这个one_pathway的各个基因了

            // 遍历 KEGGtoGeneclusters，索引 KEGG Pathway 对应的基因信息
            for (const auto& [gname, geneClusters] : onespecies.KEGGtoGeneclusters[one_pathway]) {
                
                std::set<string> uniqueStrains;  // 用来存储不重复的菌株名
                set<string> ECs_onespecies, KOs_onespecies, Reacs_onespecies;
                for (const auto& geneCluster : geneClusters) {
                    // 从这个基因名对应的geneClusters的所有gene指针的origin_strain_name插入到uniqueStrains，这表示了这个基因的菌株数
                    for (const auto& geneptr : geneCluster) {
                        uniqueStrains.insert(geneptr->origin_strain_name);
                    }

                    // 下面从每个基因cluster中收集EC, KO, reaction，然后放到对应的set中，这些set代表了一个菌种的这个基因的所有EC, KO, Reaction
                    ECs_onespecies.insert(geneCluster[0]->EC.begin(), geneCluster[0]->EC.end());
                    KOs_onespecies.insert(geneCluster[0]->kegg_ko.begin(), geneCluster[0]->kegg_ko.end());
                    Reacs_onespecies.insert(geneCluster[0]->kegg_reaction.begin(), geneCluster[0]->kegg_reaction.end());
                }
                // 现在uniqueStrains被填充完毕

                bool geneExists = false;
                //找到kpgs.genes中的相应基因
                for (auto& geneForKPGS : kpgs.genes) {
                    if (geneForKPGS.genename == gname) {
                        // 将菌株数量赋值给 species_strain_num[onespecies.species]
                        geneForKPGS.species_strain_num[onespecies.species] = uniqueStrains.size();
                        // 向相关容器添加元素
                        geneForKPGS.ECs.insert(ECs_onespecies.begin(), ECs_onespecies.end());
                        geneForKPGS.KOs.insert(KOs_onespecies.begin(), KOs_onespecies.end());
                        geneForKPGS.reactions.insert(Reacs_onespecies.begin(), Reacs_onespecies.end());

                        // 标记为已经找到
                        geneExists = true;
                        break;
                    }
                }
                // 如果上面没有在kpgs.genes中找到这恶鬼基因名对应的基因，则新建一个基因
                if (!geneExists) {
                    // 如果基因不存在，则新建一个基因条目
                    GeneforKPGS newGene;
                    newGene.genename = gname;

                    // 将菌株数量赋值给 species_strain_num[onespecies.species]
                    newGene.species_strain_num[onespecies.species] = uniqueStrains.size();

                    newGene.ECs.insert(ECs_onespecies.begin(), ECs_onespecies.end());
                    newGene.KOs.insert(KOs_onespecies.begin(), KOs_onespecies.end());
                    newGene.reactions.insert(Reacs_onespecies.begin(), Reacs_onespecies.end());

                    // 把新建的基因推入kpgs.genes
                    kpgs.genes.push_back(newGene);
                }
            }
        }

        // 向装载KeggPathwayGeneSet的vector中推入这个KeggPathwayGeneSet，用于后续进行KEGG通路图标记
        KPGSvec.push_back(kpgs);


        //现在可以输出kpgs的各基因了
        // 输出每个基因的信息
        for (const auto& geneForKPGS : kpgs.genes) {
            outputFile << "Gene: " << geneForKPGS.genename << "\t";

            // 输出每个菌种的菌株数量
            bool first = true;
            for (const auto& [species, strainCount] : geneForKPGS.species_strain_num) {
                if (!first) {
                    outputFile << " ";
                }
                first = false;
                outputFile << "[" << species << ": " << strainCount << "]";
            }

            // 输出 KO、EC 和 Reaction
            outputFile << "\tKO: ";
            for (auto it = geneForKPGS.KOs.begin(); it != geneForKPGS.KOs.end(); ++it) {
                outputFile << *it;
                if (std::next(it) != geneForKPGS.KOs.end()) {
                    outputFile << ",";
                }
            }

            outputFile << "\tEC: ";
            for (auto it = geneForKPGS.ECs.begin(); it != geneForKPGS.ECs.end(); ++it) {
                outputFile << *it;
                if (std::next(it) != geneForKPGS.ECs.end()) {
                    outputFile << ",";
                }
            }

            outputFile << "\tReaction: ";
            for (auto it = geneForKPGS.reactions.begin(); it != geneForKPGS.reactions.end(); ++it) {
                outputFile << *it;
                if (std::next(it) != geneForKPGS.reactions.end()) {
                    outputFile << ",";
                }
            }

            outputFile << std::endl;
        }
    }

    std::cout << "KEGG_Pathway_GeneSet.txt finished" << std::endl;
    outputFile.close();
}

//该函数以delimiter为分隔将string分成小份并存入vector
std::vector<std::string> splitString(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    size_t start = 0, end = 0;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        tokens.push_back(s.substr(start, end - start));
        start = end + 1;
    }
    // 最后一个
    if (start < s.size()) {
        tokens.push_back(s.substr(start));
    }
    return tokens;
}

// 通用的函数，将符合条件的字符串添加到目标容器中
template <typename T>
void addValidItems(const std::string& field, std::vector<T>& container) {
    auto items = splitString(field, ',');
    for (auto& item : items) {
        if (!item.empty() && item != "-") {
            container.push_back(item);
        }
    }
}

//解析eggnog annotations文件的单行
bool parseAnnotationLine(const std::string& line, Gene& gene) {
    // 按 '\t' 分割
    // 如果line第一个字符是"#"，即注释行，返回false
    if (line.empty() || line[0] == '#') {
        return false;
    }
    auto fields = splitString(line, '\t');
    //一般都是21个元素，少于21个很可能是错误行，目前还没见过这条语句输出
    if (fields.size() < 21) {
        std::cerr << "[WARN] parseAnnotationLine: fields.size() < 21. line=" << line << "\n";
        return false;
    }
    //下面将信息读入
    gene.gene_id = fields[0];
    gene.seed_ortholog = fields[1];

    try {
        gene.evalue = std::stod(fields[2]);
    }
    catch (...) {
        std::cerr << "evalue assignment error" << std::endl;
        gene.evalue = 0.0;
    }
    try {
        gene.score = std::stod(fields[3]);
    }
    catch (...) {
        std::cerr << "score assignment error" << std::endl;
        gene.score = 0.0;
    }

    gene.eggnog_ogs = fields[4];
    gene.max_annot_lvl = fields[5];
    gene.cog_category = fields[6];
    gene.description = fields[7];
    gene.preferred_name = fields[8];

    // 使用通用函数处理字段，减少重复代码
    addValidItems(fields[9], gene.GOs);
    addValidItems(fields[10], gene.EC);
    addValidItems(fields[11], gene.kegg_ko);
    // 现在遍历 gene.kegg_ko，去掉所有 string 头部的 ko:
    for (std::string& ko : gene.kegg_ko) {
        if (ko.find("ko:") == 0) {  // 检查是否以 "ko:" 开头
            ko.erase(0, 3);  // 移除前3个字符 ("ko:")
        }
    }
    addValidItems(fields[12], gene.kegg_pathway);
    addValidItems(fields[13], gene.kegg_module);
    addValidItems(fields[14], gene.kegg_reaction);
    addValidItems(fields[15], gene.kegg_rclass);
    addValidItems(fields[16], gene.brite);
    addValidItems(fields[17], gene.kegg_tc);
    addValidItems(fields[18], gene.cazy);
    addValidItems(fields[19], gene.bigg_reaction);
    addValidItems(fields[20], gene.pfams);

    gene.shrinkAllVectors();  // 减少内存占用
    return true;
}

// 读取单个 .annotations 文件，生成 Strain
Strain parseAnnotationsFile(const std::string& annotationsFilePath,
    const std::string& fnaFilePath,
    const string& species)
{
    Strain strain;
    strain.species = species;//确定这个菌株的菌种
    strain.strain_name = extractStrainName(fnaFilePath);//从fna文件中提取菌株名
    g_logFile << "parseAnnotationsFile(): " << annotationsFilePath << " Strain name: " << strain.strain_name << std::endl;
    std::cout<<"parseAnnotationsFile(): " << annotationsFilePath << " Strain name: " << strain.strain_name << std::endl;

    std::ifstream fin(annotationsFilePath);
    if (!fin.is_open()) {
        g_errorFile << "Failed to open annotations file: " << annotationsFilePath << std::endl;
        std::cerr << "Failed to open annotations file: " << annotationsFilePath << std::endl;
        return strain; // 返回一个空的 Strain
    }

    // 跳过前5行，这一部分我认为可以不要了，前面的函数已经引入了判断首字符'#'，前5行首字符都是'#'
    std::string line;
    for (int i = 0; i < 5; i++) {
        if (!std::getline(fin, line)) {
            g_errorFile << "Maybe empty annotations file: " << annotationsFilePath << std::endl;
            std::cerr << "Maybe empty annotations file: " << annotationsFilePath << std::endl;
            return strain;// 如果std::getline失败了，可能这个文件没有5行，很有可能是个空的
        }
    }

    // 从第6行起逐行读取
    while (true) {
        if (!std::getline(fin, line)) break;//如果读取不到一行，比如到达文件末尾，break
        if (line.empty()) continue;

        Gene g;
        if (parseAnnotationLine(line, g)) {
            // 如果 preferred_name 不是"-"，则说明有特定基因名
            if (g.preferred_name != "-") {
                g.origin_strain_name = strain.strain_name; // 赋值菌株名
                strain.geneset.push_back(g);
            }
            if (g.preferred_name == "-") {//如果没有基因名，则使用pfam部分作为临时基因名，我想不到很好的办法
                if (!g.kegg_pathway.empty() && !g.pfams.empty()) {
                    g.origin_strain_name = strain.strain_name;
                    std::string temp_preferred_name;
                    for (auto it = g.pfams.begin(); it != g.pfams.end(); ++it) {
                        temp_preferred_name += *it;
                        // 如果不是最后一个元素，添加一个逗号
                        if (std::next(it) != g.pfams.end()) {
                            temp_preferred_name += ",";
                        }
                    }
                    g.preferred_name = temp_preferred_name;
                    strain.geneset.push_back(g);
                }
            }
        }
    }

    fin.close();
    return strain;
}

// key: preferred_name（基因名）
// value: 这个基因名的同一菌种内所有菌株的 "Gene" 的集合
using GeneMap = std::map<std::string, std::vector<Gene>>;
using GeneMapptr = map<string, vector<Gene*>>;

std::ofstream groupGenesByNameLog("groupGenesByName.txt");

// 把一个菌种中所有菌株的同一个基因名对应的所有基因都聚集起来，形成从基因名映射到基因vector的map
void groupGenesByName(Species_Strain_set& sset)
{
    std::cout << "groupGenesByName(): " << sset.species << std::endl;
    g_logFile << "groupGenesByName(): " << sset.species << std::endl;
    groupGenesByNameLog << "groupGenesByName(): " << sset.species << std::endl;

    for (auto& st : sset.all_strains) { //遍历sset中所有菌株
        for (auto& g : st.geneset) { //某个菌株的geneset
            sset.gname_vec[g.preferred_name].push_back(&g); // 这个gname_vec里面保存的是基因的地址
        }
    }
    // 再次遍历，进行shrink_to_fit节省内存
    for (auto& [_, gvec] : sset.gname_vec) {
        gvec.shrink_to_fit();
    }
}

void OutputAftergroupGenesByName(Species_Strain_set& sset) {
    for (const auto& [geneName, vec] : sset.gname_vec) {
        groupGenesByNameLog << geneName<<"  "<<vec.size()<<std::endl;
    }
}
//比对两个基因是否相同
bool compareGeneFields(const Gene& g1, const Gene& g2) {
    // 除 gene_id, evalue, score 外的字段都要比对
    //直接将g1和g2的各string相加，如果不相等则返回false，这是较方便的方法
    if (g1.seed_ortholog + g1.eggnog_ogs + g1.max_annot_lvl + g1.cog_category + g1.description != g2.seed_ortholog + g2.eggnog_ogs + g2.max_annot_lvl + g2.cog_category + g2.description) return false;

    // vector类型的比较可以使用lambda函数
    // 使用 unordered_set 进行无序比较
    auto vecEq = [](const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
        return std::unordered_set<std::string>(v1.begin(), v1.end()) ==
            std::unordered_set<std::string>(v2.begin(), v2.end());
        };

    if (!vecEq(g1.GOs, g2.GOs)) return false;
    if (!vecEq(g1.EC, g2.EC)) return false;
    if (!vecEq(g1.kegg_ko, g2.kegg_ko)) return false;
    if (!vecEq(g1.kegg_pathway, g2.kegg_pathway)) return false;
    if (!vecEq(g1.kegg_module, g2.kegg_module)) return false;
    if (!vecEq(g1.kegg_reaction, g2.kegg_reaction)) return false;
    if (!vecEq(g1.kegg_rclass, g2.kegg_rclass)) return false;
    if (!vecEq(g1.brite, g2.brite)) return false;
    if (!vecEq(g1.kegg_tc, g2.kegg_tc)) return false;
    if (g1.cazy != g2.cazy) return false;
    if (g1.bigg_reaction != g2.bigg_reaction) return false;
    if (!vecEq(g1.pfams, g2.pfams)) return false;

    return true;
}

// 验证同一菌种的所有菌株的某个基因名的所有Gene是否都是一致的
// 我不知道这个函数能不能有帮助，如果不能提升效率可以不使用这个函数
// 这个函数没用到
bool checkAllSame(const std::vector<Gene>& genes) {
    if (genes.size() < 2) return true; // 只有一个或空，就算全一样
    for (size_t i = 1; i < genes.size(); i++) {
        if (!compareGeneFields(genes[0], genes[i])) {//一旦vector<Gene> genes中有一个 Gene是与众不同的，就输出False
            return false;
        }
    }
    return true;
}

// 返回所有子集，每个子集包含一组完全相同的Gene的索引
std::vector<std::vector<Gene*>> groupGenesByFields(const std::vector<Gene*>& geneVec) {

    std::vector<std::vector<Gene*>> groups;//保存所有的group，每个group是一个Gene*集合
    std::vector<bool> visited(geneVec.size(), false);//保存geneVec中每个Gene是否被访问过的情况

    for (size_t i = 0; i < geneVec.size(); i++) {
        if (visited[i]) continue;//如果被访问过则略过

        // 启动一个新子集
        std::vector<Gene*> group;
        group.push_back(geneVec[i]);
        visited[i] = true;

        // 跟后面的基因比对
        for (size_t j = i + 1; j < geneVec.size(); j++) {
            if (!visited[j]) {
                if (compareGeneFields(*geneVec[i], *geneVec[j])) {
                    group.push_back(geneVec[j]);
                    visited[j] = true;
                }
            }
        }
        group.shrink_to_fit();
        groups.push_back(group);
    }
    return groups;
}

// 使用 groupGenesByFields 基因聚类
void groupGenesforSpeciesStrainset(Species_Strain_set& sset) {
    std::cout << "groupGenesforSpeciesStrainset(): " << sset.species << std::endl;
    g_logFile << "groupGenesforSpeciesStrainset(): " << sset.species << std::endl;
    // 1. 使用 groupGenesByName 来获取所有基因名对应的基因集
    groupGenesByName(sset);// 这个函数填充了sset的map<string, vector<Gene*>> gname_vec

    OutputAftergroupGenesByName(sset);

    // 2. 遍历 geneMap 的所有基因集（vector<Gene>）
    for (auto& [geneName, genes] : sset.gname_vec) {
        // 3. 使用 groupGenesByFields 来获得所有相同的Gene*集合
        sset.clusters[geneName] = groupGenesByFields(genes);
        g_logFile << "groupGenesByFields(): " << geneName << "  groups: " << sset.clusters.size()<<" sizes of one group: ";
        for (const auto& group : sset.clusters[geneName]) {
            g_logFile << group.size() << "  ";
        }
        g_logFile << std::endl;
    }
    // 减少内存开销
    /*
    for (auto& [_, vecs] : sset.clusters) {
        for (auto& vec : vecs) {
            vec.shrink_to_fit();
        }
        vecs.shrink_to_fit();
    }*/
}

// 查找并返回第一个非空字符的索引
size_t findFirstLetter(const string& line) {
    size_t i = 0;
    while (i < line.size() && isspace(line[i])) {
        ++i;
    }
    return i;
}

// 读取 KEGG 路径表并组织为 map 结构
void loadKeggPathwayTable(const string& tableFilePath,
    map<string, string>& pathwayDescriptions,
    map<string, string>& pathwayCategories) {
    ifstream inFile(tableFilePath);
    if (!inFile) {
        std::cerr << "无法打开文件: " << tableFilePath << std::endl;
        return;
    }

    string line;
    string currentBigCategory, currentSmallCategory;

    while (getline(inFile, line)) {
        // 跳过空行
        if (line.empty()) continue;

        // 检测标题行或路径描述行
        if (line[1] == '.' && isdigit(line[2])) {
            // 小类标题行
            size_t firstLetterIdx = findFirstLetter(line);
            currentSmallCategory = line.substr(firstLetterIdx);  // 提取小类名称
        }
        else if (line[1] == '.' && !isdigit(line[2])) {
            // 大类标题行
            size_t firstLetterIdx = findFirstLetter(line);
            currentBigCategory = line.substr(firstLetterIdx);  // 提取大类名称
        }
        else if (line[0] >= '0' && line[0] <= '9') {
            // 判断路径描述行，路径代码是五位数字
            string pathwayCode = line.substr(0, 5);  // 前五位是路径代码
            string description = line.substr(6);    // 描述从第六个字符开始

            pathwayDescriptions[pathwayCode] = description;
            // 将路径与大类和小类关联
            pathwayCategories[pathwayCode] = currentBigCategory + " > " + currentSmallCategory;
        }
    }

    inFile.close();
}

// 计算某个基因 (preferred_name) 在给定菌株集合中的存在比例
// 此函数暂时还没用过
double calculatePresenceRatio(const std::string& geneName,
    const std::vector<std::set<std::string>>& strainGeneSets)//每个set<string> 代表一个菌株的Geneset
{
    size_t count = 0;
    for (auto& oneStrainSet : strainGeneSets) {//遍历每个菌株
        if (oneStrainSet.find(geneName) != oneStrainSet.end()) {
            count++;//如果在当前菌株的Geneset中找到了geneName，则count++
        }
    }
    return double(count) / strainGeneSets.size();//计算出现率
}


void parseAllSpecies(const std::string& rootFolder,
    const std::vector<std::string>& species,
    std::vector<Species_Strain_set>& strainSets) {
    for (size_t i = 0; i < species.size(); ++i) {
        std::string speciesPath = rootFolder + "\\" + species[i];
        Species_Strain_set& currentSet = strainSets[i];
        currentSet.species = species[i];

        std::cout << "Directory: " << speciesPath <<" Species: "<< currentSet.species<< std::endl;
        g_logFile << "Directory: " << speciesPath << " Species: " << currentSet.species << std::endl;

        // 遍历该菌种文件夹下的所有菌株子文件夹
        for (const auto& strainEntry : fs::directory_iterator(speciesPath)) {
            if (!strainEntry.is_directory()) continue;

            std::string strainFolder = strainEntry.path().string();
            std::string fnaFile, annoFile;

            // 遍历菌株文件夹，寻找 .fna 和 .annotations 文件
            for (const auto& fileEntry : fs::directory_iterator(strainFolder)) {
                std::string fname = fileEntry.path().filename().string();
                if (fname.ends_with(".fna")) {
                    fnaFile = fileEntry.path().string();
                }
                else if (fname.ends_with(".annotations")) {
                    annoFile = fileEntry.path().string();
                }
            }

            g_logFile << ".fna file: "<< fnaFile<<std::endl;
            g_logFile << ".annotations file: " << annoFile << std::endl;
            std::cout<<".fna file: "<< fnaFile << std::endl;
            std::cout << ".annotations file: " << annoFile << std::endl;

            // 如果都找到了，就解析这个菌株
            if (!fnaFile.empty() && !annoFile.empty()) {
                Strain st = parseAnnotationsFile(annoFile, fnaFile, species[i]);
                currentSet.all_strains.push_back(std::move(st));
            }
        }
    }
    strainSets.shrink_to_fit();
}

struct HTMLinfo {
    std::vector<std::string> geneNames;
    std::vector<std::string> ECs;
    std::vector<std::string> kos;
    std::vector<std::string> reactions;
    HTMLinfo() {}
    HTMLinfo(vector<std::string> gn, vector<std::string> EC, vector<std::string> ko, vector<std::string> re) :geneNames(gn), ECs(EC), kos(ko), reactions(re) {}
};

//记录parseHtmlFile函数的情况
std::ofstream phout("parseHtmlFile_info.txt");
//解析html文件
void parseHtmlFile(const std::string& htmlFilePath,
    // std::string& pathwayDescription,
    std::vector<RectElement>& Rects,
    // std::vector<CircleElement>& Circles,//HTML中的方块与圆圈分开
    std::vector<std::string>& all_geneNamesHTML,
    std::vector<std::string>& all_kosHTML,
    std::vector<std::string>& all_ECsHTML,
    std::vector<std::string>& all_reactionsHTML,
    bool& file_open_error) {

    std::ifstream file(htmlFilePath);
    if (!file.is_open()) {
        std::cerr << "无法打开HTML文件: " << htmlFilePath << std::endl;
        file_open_error = true;
        return;
    }

    phout << "parseHtmlFile: " << htmlFilePath << std::endl;
    //读取新的HTML文件时，之前下面这些vector需要清空，以便存入新的信息，避免与之前的HTML信息混淆
    all_geneNamesHTML.clear();
    all_kosHTML.clear();
    all_ECsHTML.clear();
    all_reactionsHTML.clear();
    Rects.clear();
    // Circles.clear();

    std::string line;
    bool inImageMap = false;
    //bool descriptionCaptured = false;

    while (std::getline(file, line)) {
        /*
        if (!descriptionCaptured && line.find("<div id=\"description\" class=\"hidden\">") != std::string::npos) {
            std::getline(file, pathwayDescription);
            descriptionCaptured = true;
            phout << "description captured: " << pathwayDescription << endl;
            //std::cout << pathwayDescription << std::endl;
        }
        */
        if (line.find("<!-- pathway image start -->") != std::string::npos) {
            inImageMap = true;
            phout << "<!-- pathway image start -->  found" << endl;
            //std::cout << "<!-- pathway image start --> found" << std::endl;
        }
        else if (inImageMap && line.find("</map>") != std::string::npos) {
            phout << "</map>  found" << endl;//找到第一个</map>代表可以结束读取
            break;
        }
        else if (inImageMap && line.find("<area id") != std::string::npos) {
            phout << "<area id  found" << endl;
            std::regex idRegex(R"(id=\"(.*?)\")");
            std::regex shapeRegex(R"(shape=\"(.*?)\")");
            std::regex coordsRegex(R"(coords=\"(.*?)\")");
            std::regex titleRegex(R"(title=\"(.*?)\")");

            std::smatch idMatch, shapeMatch, coordsMatch, titleMatch;
            std::string area_id, shape, coords, title;

            if (std::regex_search(line, idMatch, idRegex)) area_id = idMatch[1];
            if (std::regex_search(line, shapeMatch, shapeRegex)) shape = shapeMatch[1];
            if (std::regex_search(line, coordsMatch, coordsRegex)) coords = coordsMatch[1];
            if (std::regex_search(line, titleMatch, titleRegex)) title = titleMatch[1];
            if (shape == "poly") continue;//读到poly直接跳过，poly不是所需类型
            //std::cout <<"title: "<< title << std::endl;
            //std::cout << "coords: " << coords << std::endl;
            //我在这里发现有的coords非常长，很可能是圆角矩形, 这种情况下不读取，这样的长coords也不是正常的Rect，可以加判断条件跳过
            if (shape == "rect") {
                if (title.find("map") != std::string::npos) continue;//如果有"map"在title中，则跳过
                auto elem = RectElement();
                elem.area_id = area_id;
                elem.shape = Shape::RECT;
                phout << "area id: " << area_id << "  coords: " << coords << "  title: " << title << endl;
                elem.parseCoords(coords);
                elem.parseTitle(title);
                Rects.push_back(elem);//改成RectElement的vector
                phout << "Rects.push_back(elem), current Rect: " << endl;
                phout << "area id: " << elem.area_id << " coords: " << elem.point1.first << "," << elem.point1.second << "," << elem.point2.first << "," << elem.point2.second << "  ";
                if (!elem.ko_to_geneName.empty()) {
                    for (const auto& [ko, gname] : elem.ko_to_geneName) {
                        phout << ko << "--" << gname << endl;
                    }
                }
                if (!elem.EC.empty()) {
                    phout << elem.EC << endl;
                }
                if (!elem.reaction.empty()) {
                    phout << elem.reaction << endl;
                }
                phout << endl;
                //以下是填充all_kosHTML, all_geneNamesHTML, all_ECsHTML, all_reactionsHTML
                for (const auto& [ko, gene] : elem.ko_to_geneName) {
                    // Add ko to all_kosHTML regardless of the condition
                    all_kosHTML.push_back(ko);

                    // Only add gene to all_geneNamesHTML if it is different from ko
                    if (ko != gene) {
                        all_geneNamesHTML.push_back(gene);// 我发现有时有这种现象：K01222 (E3.2.1.86A), K01223 (E3.2.1.86B)
                    }
                }
                all_ECsHTML.push_back(elem.EC);
                all_reactionsHTML.push_back(elem.reaction);
                //std::cout << "Rect read: "<<elem.area_id << std::endl;
            }
            /*
            else if (shape == "circle") {
                auto elem = CircleElement();
                elem.area_id = area_id;
                elem.shape = Shape::CIRCLE;
                elem.parseCoords(coords);
                elem.parseTitle(title);
                Circles.push_back(elem);
                //std::cout << "Circle read: " << elem.area_id << std::endl;
            }
            */
        }
    }
}

std::ofstream MarkPicLog("MarkPicLog.txt");

template <typename T>
void insertMatchingElementsToSet(
    const std::vector<T>& sourceContainer, // htmlinfo 中的元素容器
    const std::set<std::string>& targetElements, // gene 中的 ECs, KOs, 或 reactions
    std::set<std::string>& outputSet) { // 存储匹配的元素指针的集合

    MarkPicLog << "In gene container: "<<std::endl;
    for (const auto& element : targetElements) {
        MarkPicLog << element << "  ";
    }
    MarkPicLog << std::endl;

    for (const auto& element : targetElements) {
        auto it = std::find(sourceContainer.begin(), sourceContainer.end(), element);
        if (it != sourceContainer.end()) {
            outputSet.insert(*it);  // 将匹配元素的地址插入 set
            MarkPicLog << *it << " inserted" << std::endl;
        }
    }
}

void fillRectWithoutBlackText(cv::Mat& image, const cv::Rect& rect, const cv::Scalar& color) {
    // 复制ROI区域
    cv::Mat roi = image(rect);
    cv::Mat roiCopy = roi.clone();

    // 转为灰度图
    cv::Mat gray;
    cv::cvtColor(roi, gray, cv::COLOR_BGR2GRAY);

    // 创建mask，排除黑色文字（比如灰度值 <= 50 的像素）
    cv::Mat nonBlackMask;
    cv::inRange(gray, 51, 255, nonBlackMask);  // 灰度 > 50 认为不是文字

    // 创建纯色背景图（要填充的颜色）
    cv::Mat colorMat(roi.size(), roi.type(), color);

    // 只填非黑色区域
    colorMat.copyTo(roiCopy, nonBlackMask);

    // 贴回原图
    roiCopy.copyTo(image(rect));
}

void InsideMarkPic(const std::map<string, set<string>>& ECtoSpecies,
    const std::map<string, set<string>>& KOtoSpecies,
    const std::map<string, set<string>>& ReactoSpecies,
    const std::vector<RectElement>& Rects,
    const std::string& PNG) {

    // 首先读取PNG
    cv::Mat image = cv::imread(PNG);
    // 检查是否成功打开
    if (image.empty()) {
        std::cerr << "Image: " << PNG << " open failed" << endl;
    }
    else { // 成功打开图像后，进行标注
        MarkPicLog << PNG << " open succeeded" << std::endl;
        // 可以遍历ECtoSpecies, KOtoSpecies, ReactoSpecies，一个Rect会有多个KO，一个EC和一个Reaction，
        // 如果某个Rect的至少一个KO, EC, Reaction都在上面三个map中出现
        for (const auto& rect : Rects) {
            // 我的想法是：遍历Rects时，每个rect都有ko_to_geneName, EC, reaction
            // 一个rect的多个ko中可能会有一个或几个ko是在KOtoSpecies中能找到的，而且对应了一种或多种菌种
            // 然后rect的EC, Reaction都是string，把EC和Reaction也进行一样的操作，在ECtoSpecies和ReactoSpecies中查找
            // 然后同样rect的EC和Reaction在ECtoSpecies和ReactoSpecies中可能对应一个或多个菌种
            // 如果一个rect的全部或部分ko, EC, Reaction都对应了同一个或相同的多个菌种，那么这个rect可以在PNG中标注了，标注时依据菌种的数量在方框中填充颜色，如果只有一个菌种则填充一种颜色，多个菌种则把方框用数线划分，然后在分出的每份中填充不同颜色，注意一个颜色代表一种菌种
            std::set<std::string> speciesSet;  // 用于存储所有匹配的菌种
            bool KOfound = false,ECfound = false,Reacfound = false;
            MarkPicLog << "rect: " << rect.area_id << "  " << rect.EC << "  " << rect.reaction << std::endl;
            for (const auto& [ko,genename] : rect.ko_to_geneName) {
                MarkPicLog << ko << "->" << genename<<std::endl;
            }


            // 检查KO，EC和reaction是否能在相应的map中找到
            // 遍历每个KO，查看KO是否在KOtoSpecies中
            for (const auto& [ko, gname] : rect.ko_to_geneName) {
                if (KOtoSpecies.find(ko) != KOtoSpecies.end()) {
                    // 添加对应菌种到speciesSet
                    const auto& speciesList = KOtoSpecies.at(ko);
                    speciesSet.insert(speciesList.begin(), speciesList.end());
                    MarkPicLog << ko << " found in KOtoSpecies ";
                    for (const auto& one : speciesList) {
                        MarkPicLog << one << "  ";
                    }
                    MarkPicLog << std::endl;
                    KOfound = true;
                }
            }

            // 检查EC是否在ECtoSpecies中
            if (ECtoSpecies.find(rect.EC) != ECtoSpecies.end()) {
                // 添加对应菌种到speciesSet
                const auto& speciesList = ECtoSpecies.at(rect.EC);
                speciesSet.insert(speciesList.begin(), speciesList.end());
                MarkPicLog << rect.EC << " found in ECtoSpecies ";
                for (const auto& one : speciesList) {
                    MarkPicLog << one << "  ";
                }
                MarkPicLog << std::endl;
                ECfound = true;
            }

            // 检查reaction是否在ReactoSpecies中
            if (ReactoSpecies.find(rect.reaction) != ReactoSpecies.end()) {
                // 添加对应菌种到speciesSet
                const auto& speciesList = ReactoSpecies.at(rect.reaction);
                speciesSet.insert(speciesList.begin(), speciesList.end());
                MarkPicLog << rect.reaction << " found in ReactoSpecies ";
                for (const auto& one : speciesList) {
                    MarkPicLog << one << "  ";
                }
                MarkPicLog << std::endl;
                Reacfound = true;
            }

            // 如果找到了相应的菌种，那么就对该区域进行标记
            if (KOfound&& ECfound&&Reacfound && !speciesSet.empty()) {
                MarkPicLog <<std::endl<< "KOfound&& ECfound&&Reacfound, begin to mark" << std::endl;
                // 计算颜色
                cv::Scalar color1(0, 0, 0); // 默认黑色
                cv::Scalar color2(0, 0, 0);
                cv::Scalar color3(0, 0, 0);

                // 在switch语句外部定义
                int middleX = 0, middleX2 = 0;

                // 判断speciesSet中包含的菌种，选择颜色
                switch (speciesSet.size()) {
                    // 只有一个菌种的情况
                case 1: {
                      MarkPicLog << "One species";
                      if (speciesSet.count("Cris")) {
                          color1 = cv::Scalar(0, 255, 0);  // 绿色
                          MarkPicLog << "Cris" << std::endl;
                      }
                      else if (speciesSet.count("Iners")) {
                          color1 = cv::Scalar(230, 216, 173); // 浅蓝
                          MarkPicLog << "Iners" << std::endl;
                      }
                      else if (speciesSet.count("Gard")) {
                          color1 = cv::Scalar(255, 0, 255);  // 紫色
                          MarkPicLog << "Gard" << std::endl;
                      }
                      cv::Rect fillRect(rect.point1.first, rect.point1.second,
                          rect.point2.first - rect.point1.first,
                          rect.point2.second - rect.point1.second);
                      // 使用颜色填充方框
                      /*
                      cv::rectangle(image,
                          cv::Point(rect.point1.first, rect.point1.second),
                          cv::Point(rect.point2.first, rect.point2.second),
                          color1, cv::FILLED);  // 使用选定的颜色填充方框
                      */
                      fillRectWithoutBlackText(image, fillRect, color1);
                      break;
                }
                case 2: {
                    // 两个菌种的情况
                    MarkPicLog << "2 species" << std::endl;

                    // 默认分配绿色和蓝色
                    color1 = cv::Scalar(0, 255, 0);  // 绿色左半部分
                    color2 = cv::Scalar(230, 216, 173); // 浅蓝

                    // 检查speciesSet中包含的菌种，动态设置颜色
                    if (speciesSet.count("Cris") && speciesSet.count("Iners")) {
                        color1 = cv::Scalar(0, 255, 0);  // 绿色
                        color2 = cv::Scalar(230, 216, 173); //浅蓝
                        MarkPicLog << "Cris + Iners" << std::endl;
                    }
                    else if (speciesSet.count("Cris") && speciesSet.count("Gard")) {
                        color1 = cv::Scalar(0, 255, 0);  // 绿色
                        color2 = cv::Scalar(255, 0, 255);  // 紫色
                        MarkPicLog << "Cris + Gard" << std::endl;
                    }
                    else if (speciesSet.count("Iners") && speciesSet.count("Gard")) {
                        color1 = cv::Scalar(230, 216, 173); // 浅蓝
                        color2 = cv::Scalar(255, 0, 255);  // 紫色
                        MarkPicLog << "Iners + Gard" << std::endl;
                    }

                    // 绘制绿色和蓝色两部分
                    middleX = (rect.point1.first + rect.point2.first) / 2;
                    cv::Rect leftRect(rect.point1.first, rect.point1.second,
                        middleX - rect.point1.first,
                        rect.point2.second - rect.point1.second);
                    cv::Rect rightRect(middleX, rect.point1.second,
                        rect.point2.first - middleX,
                        rect.point2.second - rect.point1.second);

                    fillRectWithoutBlackText(image, leftRect, color1);
                    fillRectWithoutBlackText(image, rightRect, color2);
                    /*
                    cv::rectangle(image,
                        cv::Point(rect.point1.first, rect.point1.second),
                        cv::Point(middleX, rect.point2.second),
                        color1, cv::FILLED);  // 左半部分
                    cv::rectangle(image,
                        cv::Point(middleX, rect.point1.second),
                        cv::Point(rect.point2.first, rect.point2.second),
                        color2, cv::FILLED);  // 右半部分
                    */
                    break;
                }
                case 3: {

                    // 三个菌种的情况
                    MarkPicLog << "3 species" << std::endl;

                    // 默认分配绿色、蓝色和紫色
                    color1 = cv::Scalar(0, 255, 0);  // 绿色
                    color2 = cv::Scalar(230, 216, 173); //浅蓝
                    color3 = cv::Scalar(255, 0, 255);  // 紫色

                    // 如果speciesSet中包含所有三个菌种
                    if (speciesSet.count("Cris") && speciesSet.count("Iners") && speciesSet.count("Gard")) {
                        MarkPicLog << "Cris + Iners + Gard" << std::endl;
                    }
                    int x0 = rect.point1.first;
                    int x1 = (2 * rect.point1.first + rect.point2.first) / 3;
                    int x2 = (rect.point1.first + 2 * rect.point2.first) / 3;
                    int x3 = rect.point2.first;

                    cv::Rect part1(x0, rect.point1.second, x1 - x0, rect.point2.second - rect.point1.second);
                    cv::Rect part2(x1, rect.point1.second, x2 - x1, rect.point2.second - rect.point1.second);
                    cv::Rect part3(x2, rect.point1.second, x3 - x2, rect.point2.second - rect.point1.second);

                    fillRectWithoutBlackText(image, part1, color1);
                    fillRectWithoutBlackText(image, part2, color2);
                    fillRectWithoutBlackText(image, part3, color3);
                    /*
                    // 计算每部分的分割位置
                    middleX = (rect.point1.first + rect.point2.first) / 3;
                    middleX2 = 2 * middleX;

                    // 绘制绿色、蓝色和紫色三个部分
                    cv::rectangle(image,
                        cv::Point(rect.point1.first, rect.point1.second),
                        cv::Point(middleX, rect.point2.second),
                        color1, cv::FILLED);  // 左部分绿色
                    cv::rectangle(image,
                        cv::Point(middleX, rect.point1.second),
                        cv::Point(middleX2, rect.point2.second),
                        color2, cv::FILLED);  // 中间部分蓝色
                    cv::rectangle(image,
                        cv::Point(middleX2, rect.point1.second),
                        cv::Point(rect.point2.first, rect.point2.second),
                        color3, cv::FILLED);  // 右部分紫色
                    */
                    break;
                }
                case 4:
                    MarkPicLog << "4 species" << std::endl;
                    break;
                case 5:
                    // 五个菌种的情况（可以继续添加颜色并进行五分）
                    MarkPicLog << "5 species" << std::endl;
                    break;

                default:
                    // 如果speciesSet.size()不在1-5之间
                    MarkPicLog << "Unexpected number of species, maybe speciesSet.size()==0 or speciesSet.size()>5" << std::endl;
                    break;
                }
            }

        }


        // 提取文件名部分 (例如 map00010) 
        size_t startPos = PNG.find("map");  // 找到"map"的起始位置
        size_t endPos = PNG.find("@2x");   // 找到"@2x"的结束位置

        if (startPos != std::string::npos && endPos != std::string::npos) {
            std::string pathwayID = PNG.substr(startPos, endPos - startPos);  // 提取 "map00010"

            // 构造保存路径
            std::string savePath = "D:\\Bioana\\eggnogparsev6\\marked\\marked_" + pathwayID + ".png";
            // 保存标记后的图像
            cv::imwrite(savePath, image);
            std::cout << "Marked image saved at: " << savePath << std::endl;
        }
        else {
            std::cerr << "Invalid PNG path: " << PNG << std::endl;
        }
    }
}


void MarkPic(const set<string>& all_code,const vector<KeggPathwayGeneSet>& KPGSvec,const std::vector<Species_Strain_set>& AllSpecies) {
    // 遍历all_KEGGPathway_code，对每个code提取对应的png与源码，然后读取一下html的信息
    // 首先要读取KEGG_Pathway的png图片和html源码
    // 然后遍历KPGSvec，每个KeggPathwayGeneSet代表一个KEGG_Pathway的基因集和
    // 遍历一个KeggPathwayGeneSet的genes时，为每个基因所在的酶方框进行标记，标记前要检查：
    // 一个gene的EC, KO, Reaction是否在这个KEGGPathway的html中，即这个KEGGPathway的所有酶中有哪些是和这个gene对应的，找到对应的酶方框，有可能有多个方框
    // 这个（些）方框的EC, KO, 反应号各自在哪些菌种中有
    // 根据有相应EC, KO, 反应号的菌种，为方框标颜色，一个颜色代表一个菌种，多个菌种共有这个酶的情况，则标多个颜色
    for (const string& cur_code : all_code) {
        string HTML(KEGG_Pathway_HTML_PNG+"\\kegg_pathway_html\\"+cur_code+".html");
        string PNG(KEGG_Pathway_HTML_PNG + "\\kegg_pathway_images\\" + cur_code + "@2x.png");
        //读取HTML，组织Rects
        vector<RectElement> Rects;// Rects保存一个KEGGPathway的HTML源码中所有的酶方框
        HTMLinfo htmlinfo;// 这个结构体中保存一个HTML源码中出现过的所有酶方框的基因名
        bool file_open_error=false;
        parseHtmlFile(HTML, Rects, htmlinfo.geneNames, htmlinfo.kos, htmlinfo.ECs, htmlinfo.reactions, file_open_error);
        // 现在Rects与html已经被填充
        if (!file_open_error) { // 如果HTML源码成功打开
            // 先遍历KPGSvec，找到对应的KeggPathwayGeneSet
            for ( auto& KPGS : KPGSvec) { // 这里不是打算遍历，只是想在KPGSvec中找到对应当前cur_code的KPGS，因此后面搭配了break;
                if (KPGS.kegg_pathway == cur_code) { // 但是我觉得也可以不这样做，直接遍历KPGSvec，走到哪个编号就读哪个

                    //确认htmlinfo.geneNames, htmlinfo.kos, htmlinfo.ECs, htmlinfo.reactions被正确填充
                    MarkPicLog << "KEGGPathway: " << cur_code << std::endl<<std::endl<<"htmlinfo.geneNames:";
                    for (const auto& genename : htmlinfo.geneNames) {
                        MarkPicLog << genename << "  ";
                    }
                    MarkPicLog << std::endl <<"htmlinfo.kos:"<<std::endl;
                    for (const auto& ko : htmlinfo.kos) {
                        MarkPicLog << ko << "  ";
                    }
                    MarkPicLog << std::endl << "htmlinfo.ECs:" << std::endl;
                    for (const auto& EC : htmlinfo.ECs) {
                        MarkPicLog << EC << "  ";
                    }
                    MarkPicLog << std::endl << "htmlinfo.reactions" << std::endl;
                    for (const auto& reac : htmlinfo.reactions) {
                        MarkPicLog << reac << "  ";
                    }
                    MarkPicLog << std::endl;

                    map<string, set<string>> ECtoSpecies,KOtoSpecies,ReactoSpecies;// 这些map用于保存这个KEGGPathway中EC, KO, 反应号及它们对应的菌种，用于标记
                    // 开始遍历KPGS的genes
                    for ( auto& gene : KPGS.genes) {
                        // 首先要找到这个gene的ECs, KOs, Reac在htmlinfo中有哪些是对应的
                        set<std::string> ECinHTML, KOfoundinHTML, ReacfoundinHTML;

                        // 使用泛型函数处理 ECs, KOs 和 reactions
                        // ECinHTML，KOfoundinHTML，ReacfoundinHTML没有被正确填充
                        MarkPicLog << "insertMatchingElementsToSet()  EC"<<std::endl;
                        insertMatchingElementsToSet(htmlinfo.ECs, gene.ECs, ECinHTML);
                        MarkPicLog << "insertMatchingElementsToSet()  KO" << std::endl;
                        insertMatchingElementsToSet(htmlinfo.kos, gene.KOs, KOfoundinHTML);
                        MarkPicLog << "insertMatchingElementsToSet()  Reactions" << std::endl;
                        insertMatchingElementsToSet(htmlinfo.reactions, gene.reactions, ReacfoundinHTML);
                        MarkPicLog << "ECinHTML.size(): " << ECinHTML.size() << " KOfoundinHTML.size(): " << KOfoundinHTML.size() << " ReacfoundinHTML.size(): " << ReacfoundinHTML.size() << std::endl;

                        MarkPicLog << "Gene: " << gene.genename<<"  ";
                        for (const auto& EC : gene.ECs) {
                            MarkPicLog << EC << "  ";
                        }
                        for (const auto& KO : gene.KOs) {
                            MarkPicLog << KO << "  ";
                        }
                        for (const auto& Reac : gene.reactions) {
                            MarkPicLog << Reac << "  ";
                        }
                        MarkPicLog << std::endl;
                        // 确认当前的gene的EC, KO, Reactions在HTML中都出现过，这样才进入标记环节
                        if (ECinHTML.size() > 0 && KOfoundinHTML.size() > 0 && ReacfoundinHTML.size() > 0) {

                            MarkPicLog << "ECinHTML.size()>0&&KOfoundinHTML.size()>0&&ReacfoundinHTML.size()>0  true"<<std::endl;
                            // 然后把ECinHTML，KOfoundinHTML，ReacfoundinHTML的各元素在Species_Strain_set的ECtoStrains,KOtoStrains,ReactiontoStrains中寻找
                            // 如果找到而且菌株数足够大，则认为可以标记
                            MarkPicLog << "ECinHTML: ";
                            for (const auto& EC : ECinHTML) {
                                MarkPicLog << EC << "  ";
                            }
                            MarkPicLog << std::endl;
                            for (const auto& EC : ECinHTML) {
                                for (const Species_Strain_set& OneSpecies : AllSpecies) {
                                    // 确保找到了对应的 EC 并且它有足够的菌株
                                    auto it = OneSpecies.ECtoStrains.find(EC);
                                    if (it != OneSpecies.ECtoStrains.end() && it->second.size() > OneSpecies.all_strains.size()*0.5) {
                                        // 此时这个EC号是可以标记的
                                        ECtoSpecies[EC] .insert( OneSpecies.species);
                                        MarkPicLog << EC <<"->"<<OneSpecies.species<<" Pushed into ECtoSpecies"<<std::endl;
                                    }
                                }
                            }
                            MarkPicLog << "KOfoundinHTML: ";
                            for (const auto& KO : KOfoundinHTML) {
                                MarkPicLog << KO << "  ";
                            }
                            MarkPicLog << std::endl;
                            for (const auto& KO : KOfoundinHTML) {
                                for (const Species_Strain_set& OneSpecies : AllSpecies) {
                                    auto it = OneSpecies.KOtoStrains.find(KO);
                                    if (it != OneSpecies.KOtoStrains.end() && it->second.size() > OneSpecies.all_strains.size() * 0.5) {
                                        KOtoSpecies[KO].insert(OneSpecies.species);
                                        MarkPicLog << KO << "->" << OneSpecies.species << " Pushed into KOtoSpecies" << std::endl;
                                    }
                                }
                            }
                            MarkPicLog << "ReacfoundinHTML: ";
                            for (const auto& Reac : ReacfoundinHTML) {
                                MarkPicLog << Reac << "  ";
                            }
                            MarkPicLog << std::endl;
                            for (const auto& Reac : ReacfoundinHTML) {
                                for (const Species_Strain_set& OneSpecies : AllSpecies) {
                                    auto it = OneSpecies.ReactiontoStrains.find(Reac);
                                    if (it != OneSpecies.ReactiontoStrains.end() && it->second.size() > OneSpecies.all_strains.size() * 0.5) {
                                        ReactoSpecies[Reac].insert(OneSpecies.species);
                                        MarkPicLog << Reac << "->" << OneSpecies.species << " Pushed intoReactoSpecies" << std::endl;
                                    }
                                }
                            }
                        }
                    }

                    // 调用InsideMarkPic前可以遍历一下ECtoSpecies,KOtoSpecies,ReactoSpecies
                    MarkPicLog << std::endl<<"ECtoSpecies: "<<std::endl;
                    for (const auto& [EC, species] : ECtoSpecies) {
                        MarkPicLog << EC << "->";
                        for (const auto& one : species) {
                            MarkPicLog << one << "  ";
                        }
                        MarkPicLog << std::endl;
                    }
                    MarkPicLog << std::endl << "KOtoSpecies: " << std::endl;
                    for (const auto& [KO, species] : KOtoSpecies) {
                        MarkPicLog << KO << "->";
                        for (const auto& one : species) {
                            MarkPicLog << one << "  ";
                        }
                        MarkPicLog << std::endl;
                    }
                    MarkPicLog << std::endl << "ReactoSpecies: " << std::endl;
                    for (const auto& [reac, species] : ReactoSpecies) {
                        MarkPicLog << reac << "->";
                        for (const auto& one : species) {
                            MarkPicLog << reac << "  ";
                        }
                        MarkPicLog << std::endl;
                    }

                    MarkPicLog << std::endl << "Starting InsideMarkPic()..." << std::endl;
                    // 现在依据ECtoSpecies,KOtoSpecies,ReactoSpecies进行标记
                    // 需要一个函数，读取ECtoSpecies,KOtoSpecies,ReactoSpecies，Rects，PNG，在.png中标记
                    // 先在Rects中找到ECtoSpecies,KOtoSpecies,ReactoSpecies中每个EC, KO, Reaction对应的Rect，读取坐标，坐标*2，然后根据它们对应的菌种名在图中标注
                    InsideMarkPic(ECtoSpecies,KOtoSpecies,ReactoSpecies,Rects, PNG);

                    break;// 找到并处理完成后就可以break，无需继续遍历
                }// if (KPGS.kegg_pathway == cur_code)
            }// for ( auto& KPGS : KPGSvec)
        }// if (!file_open_error)
    }//for (const string& cur_code : all_code)
}// 函数结尾

int main(int argc, char* argv[]) {
    // 设定总文件夹路径
    const std::string rootFolder = "D:\\CrisIners";
    
    // 获取菌种名称列表
    std::vector<std::string> species = getSpeciesFromFolder(rootFolder);

    // 使用std::vector<Strain_set>来存储所有菌种的信息
    std::vector<Species_Strain_set> strainSets(species.size());

    // 然后需要读取各菌种文件夹
    // 需要一个函数，读取rootFolder之下每个文件夹，这些文件夹里有每个菌株的文件夹，菌株文件夹下会有一个.fna文件和一个.annotations文件
    // 这个函数利用parseAnnotationsFile函数为每个菌株生成一个Strain对象，然后把同一个菌种内的所有菌株的Strain组织到Species_Strain_set中
    // 这个函数应接受参数std::string rootFolder，std::vector<std::string> species，这两个参数用来拼接成菌种文件夹，还需要参数std::vector<Species_Strain_set>&，把strainSets传入，把生成的每个Species_Strain_set都放到strainSets中
    parseAllSpecies(rootFolder, species, strainSets);

    for (auto& sset : strainSets) {
        groupGenesforSpeciesStrainset(sset); // 对sset的gname_vec和clusters修改
        modifyKEGGtoGeneclusters(sset);// 修改sset的KEGGtoGeneclusters
    }
    // 我打算在这里添加日志输出的代码，查看一下sset的gname_vec, clusters, KEGGtoGeneclusters
    // 写一个函数，

    set<string> all_KEGGPathway_code;

    vector<KeggPathwayGeneSet> KPGSvec;
    // 5. 生成 KEGG_Pathway_GeneSet.txt 文件
    GenerateKeggPathwaySet(strainSets, all_KEGGPathway_code, KPGSvec);

    // 把strainSets中的每个Species_Strain_set都设置好
    ModifyECKOReac(strainSets);

    // 现在需要进行标注，只剩下这一个函数的问题需要解决了
    MarkPic(all_KEGGPathway_code, KPGSvec, strainSets);
    // InsideMarkPic这个函数有问题，遍历Rects时，每次遍历时都会进行一次标记，每次标记都会覆盖之前的标记图，最后留下来的图上什么标记都没有
    // 我在控制台输出中发现了很多Marked image saved at:，这就说明了这一点

    std::cout << "Finished processing all species and generated KEGG_Pathway_GeneSet.txt." << std::endl;
    return 0;
}
