#include <iomanip>
#include <memory>
#include <iostream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <memory>
#include <regex>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <chrono>
#include <list>
#include <filesystem>
#include "CLI/App.hpp"
#include "CLI/Validators.hpp"
#include "spdlog/spdlog.h"
#include "CLI/CLI.hpp"
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>    

namespace fs = std::filesystem;

struct gapped_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::gapped<seqan3::dna5>;
    using sequence_legal_alphabet = seqan3::gapped<seqan3::dna5>;
};

using seqloc = std::pair<std::string, unsigned int>;

//struct key_hash : public std::unary_function<seqloc, std::size_t>
//{
 //std::size_t operator()(const seqloc& k) const
 //{
   //return std::get<0>(k) ^ std::get<1>(k);
 //}
//};

int main(int argc, char** argv){
    /*
     *Command line parsing
     */
    CLI::App app{"RAMBO application"};
    //app.require_subcommand(1, 1);
    //std::vector<fs::path> input_files;
    fs::path query_file = "";
    fs::path truth_file = "";
    fs::path output_prefix("vm-results");
    bool verbose = false;
    unsigned int num_threads = 1;

    // Generic flags
    app.add_option(
            "query-file",
            query_file,
            "Query MSA file"
    )->required()->check(CLI::ExistingFile);
    app.add_option(
            "truth-file",
            truth_file,
            "Truth MSA file"
    //)->required()->check(CLI::ExistingFile);
    );
    app.add_flag(
            "-v,--verbose",
            verbose,
            "Display debug output"
    );
    app.add_option(
            "-p,--threads",
            num_threads,
            "Number of threads to use"
    );
    app.add_option(
            "-o,--output",
            output_prefix,
            "Prefix of output directory"
    );

    CLI11_PARSE(app, argc, argv);

    /*
     *Constants and assertions 
     */
    spdlog::set_level(verbose ? spdlog::level::debug : spdlog::level::info); // Set global log level to debug
    //omp_set_num_threads(num_threads);
    

    seqan3::sequence_file_input<gapped_traits> truth_file_h{truth_file};
    int curr_block = -1;
    std::vector<std::shared_ptr<std::set<seqloc>>> groups; 
    std::map<seqloc, std::shared_ptr<std::set<seqloc>>> loc_to_idx; 
    for (auto & [seq, id, qual] : truth_file_h)
    {
        //seqan3::debug_stream << "ID:     " << id << '\n';
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        int start, strand, block_id;
        if(std::regex_search(s, matches, rgx)) {
            seq_id = matches[1].str();
            block_id = std::stoi(matches[2].str());
            start = std::stoi(matches[3].str());
            strand = matches[5].str() == "-" ? -1 : 1;
        } else {
            std::cout << "Match not found\n";
            continue;
        }
        if (block_id != curr_block) {
            curr_block = block_id;
            groups = std::vector<std::shared_ptr<std::set<seqloc>>>(seq.size());
            for (auto i = 0; i < seq.size(); ++i) {
                groups[i] = std::make_shared<std::set<seqloc>>();
            }
        }
        unsigned int ng = 0;
        for (auto col_idx = 0; col_idx < seq.size(); ++col_idx) {
            if (seq[col_idx] == seqan3::gap{})
                continue;
            ++ng;
            seqloc sl = std::make_pair(seq_id, start + ng);
            groups[col_idx]->insert(sl);
            loc_to_idx[sl] = groups[col_idx];
            //seqan3::debug_stream << *loc_to_idx[std::make_pair(seq_id, start + strand*ng)] << std::endl;
        }
    }

    seqan3::sequence_file_input<gapped_traits> query_file_h{query_file};
    curr_block = -1;
    std::vector<std::shared_ptr<std::vector<std::set<seqloc>>>> block_groups;
    std::shared_ptr<std::vector<std::set<seqloc>>> query_groups; 
    for (auto & [seq, id, qual] : query_file_h)
    {
        //seqan3::debug_stream << "ID:     " << id << '\n';
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        int start, strand, block_id;
        if(std::regex_search(s, matches, rgx)) {
            seq_id = matches[1].str();
            block_id = std::stoi(matches[2].str());
            start = std::stoi(matches[3].str());
            strand = matches[5].str() == "-" ? -1 : 1;
        } else {
            std::cout << "Match not found\n";
            continue;
        }
        if (block_id != curr_block) {
            //score block
            //
            //
            if (curr_block != -1) {
                block_groups.push_back(query_groups);
            }
            //wipe block data structures
            curr_block = block_id;
            query_groups = std::make_shared<std::vector<std::set<seqloc>>>(seq.size());
        }
        unsigned int ng = 0;
        for (auto col_idx = 0; col_idx < seq.size(); ++col_idx) {
            if (seq[col_idx] == seqan3::gap{})
                continue;
            ++ng;
            (*query_groups)[col_idx].insert(std::make_pair(seq_id, start + ng));
        }
    }


    // Compute intersections
    unsigned int fp = 0, fn = 0, tp = 0;
    for (auto block_group : block_groups) {
        for (auto query_group : *block_group) {
            //seqan3::debug_stream << query_group << std::endl;
            for (seqloc sl : query_group) {
                auto truth_group = *loc_to_idx[sl];
                std::vector<seqloc> v_intersection;
                std::set_intersection(truth_group.begin(), truth_group.end(), query_group.begin(), query_group.end(), std::back_inserter(v_intersection));
                tp += v_intersection.size();
                v_intersection.clear();
                std::set_difference(truth_group.begin(), truth_group.end(), query_group.begin(), query_group.end(), std::back_inserter(v_intersection));
                fn += v_intersection.size();
                //if (!v_intersection.empty())
                    //seqan3::debug_stream << truth_group << std::endl;
                v_intersection.clear();
                std::set_difference(query_group.begin(), query_group.end(), truth_group.begin(), truth_group.end(), std::back_inserter(v_intersection));
                fp += v_intersection.size();
                //if (!v_intersection.empty())
                    //seqan3::debug_stream << truth_group << std::endl;
                v_intersection.clear();
            }
        }
    }
    spdlog::info("TP = {}, FN = {}, FP = {}", tp, fn, fp);
}
