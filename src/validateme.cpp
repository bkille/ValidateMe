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
#include <bitset>
#include <omp.h>
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
#include "tqdm.h"

namespace fs = std::filesystem;

struct gapped_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::gapped<seqan3::dna5>;
    using sequence_legal_alphabet = seqan3::gapped<seqan3::dna5>;
};

using seqloc = std::pair<std::size_t, unsigned int>;

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
    )->required()->check(CLI::ExistingFile);
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
    omp_set_num_threads(num_threads);
    
    // Get block count for both files
    std::ifstream inFile(truth_file); 
    unsigned int truth_seq_count = std::count(std::istreambuf_iterator<char>(inFile), 
        std::istreambuf_iterator<char>(), '>'); 
    inFile.close();
    inFile = std::ifstream(query_file);
    unsigned int query_seq_count = std::count(std::istreambuf_iterator<char>(inFile), 
        std::istreambuf_iterator<char>(), '>'); 
    inFile.close();
    seqan3::sequence_file_input<gapped_traits> truth_file_h{truth_file};
    unsigned int truth_blocks = 0;
    unsigned int query_blocks = 0;
    auto truth_sequences = std::set<std::string>{};
    auto seq_to_idx = std::map<std::string, int>{};
    spdlog::debug("Counting truth blocks...");
    tqdm bar;
    unsigned int idx = 0;
    #pragma omp parallel
    #pragma omp single
    {
    for (auto & [seq, id, qual] : truth_file_h)
    {
        idx += 1;
        #pragma omp task firstprivate(seq, id, qual, idx)
        {
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        if(std::regex_search(s, matches, rgx)) {
            auto truth_blocks_thread = std::stoi(matches[2].str());
            seq_id = matches[1].str();
            #pragma omp critical 
            {
            bar.progress(idx, truth_seq_count);
            truth_sequences.insert(seq_id);
            truth_blocks = truth_blocks_thread > truth_blocks 
                ? truth_blocks_thread : truth_blocks;
            }
        } else {
            spdlog::error("{} does not fit the header format", id);
        }
        }
    }
    }
    bar.finish();
    bar.reset();
    idx = 0;
    for (auto seq_id : truth_sequences) {
        spdlog::debug("{}: {}", idx, seq_id);
        seq_to_idx[seq_id] = idx++;
    }
    spdlog::debug("Truth blocks: {}", truth_blocks);
    spdlog::debug("Counting query blocks...");
    seqan3::sequence_file_input<gapped_traits> query_file_h{query_file};
    idx = 0;

    #pragma omp parallel
    #pragma omp single
    {
    for (auto & [seq, id, qual] : query_file_h)
    {
        idx += 1;
        bar.progress(idx, query_seq_count);
        #pragma omp task firstprivate(seq, id, qual)
        {
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        if(std::regex_search(s, matches, rgx)) {
            auto query_blocks_thread = std::stoi(matches[2].str());
            #pragma omp critical 
            {
            query_blocks = query_blocks_thread > query_blocks 
                ? query_blocks_thread : query_blocks;
            }
        } else {
            spdlog::error("{} does not fit the header format", id);
        }
        }
    }
    }
    bar.finish();
    bar.reset();
    spdlog::debug("Query blocks: {}", query_blocks);


    // Create homology groups and mappings for the truth set
    spdlog::debug("Creating truth data homology groups");
    truth_file_h = seqan3::sequence_file_input<gapped_traits>(truth_file);
    using GroupVec = std::vector<std::shared_ptr<std::set<seqloc>>>; 
    std::vector<GroupVec> truth_block_vec = std::vector<GroupVec>(truth_blocks);
    std::map<seqloc, std::shared_ptr<std::set<seqloc>>> loc_to_idx; 
    std::map<seqloc, unsigned int> loc_to_start; 
    unsigned int seq_count = 0;
    #pragma omp parallel
    #pragma omp single
    {
    for (auto & [seq, id, qual] : truth_file_h) {
        #pragma omp task firstprivate(seq, id, qual)
        {
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        unsigned int start, block_id, length;
        int strand;
        if(std::regex_search(s, matches, rgx)) {
            seq_id = matches[1].str();
            block_id = std::stoi(matches[2].str()) - 1;
            start = std::stoi(matches[3].str());
            strand = matches[5].str() == "-" ? -1 : 1;
            length = std::stoi(matches[6].str());
            if (strand == -1) {
                start = length - start + 1;
            }
        }
        if (truth_block_vec[block_id].empty()) {
            #pragma omp critical 
            {
            if (truth_block_vec[block_id].empty()) {
                truth_block_vec[block_id].resize(seq.size());
                for (auto i = 0; i < seq.size(); ++i) {
                    truth_block_vec[block_id][i] = std::make_shared<std::set<seqloc>>();
                }
            }
            }
        }
        unsigned int ng = 0;
        for (auto col_idx = 0; col_idx < seq.size(); ++col_idx) {
            if (seq[col_idx] == seqan3::gap{})
                continue;
            ++ng;
            seqloc sl = std::make_pair(seq_to_idx[seq_id], start + ng*strand);
            #pragma omp critical 
            {
            truth_block_vec[block_id][col_idx]->insert(sl);
            loc_to_idx[sl] = truth_block_vec[block_id][col_idx];
            loc_to_start[sl] = start;
            }
        }
        #pragma omp critical 
        {
        seq_count += 1;
        bar.progress(seq_count, truth_seq_count);
        }
        }
    }
    }
    bar.finish();
    bar.reset();

    // Create homologies for query
    spdlog::debug("Creating query data homology groups");
    query_file_h = seqan3::sequence_file_input<gapped_traits>(query_file);
    std::vector<GroupVec> query_block_vec = std::vector<GroupVec>(query_blocks);
    seq_count = 0;
    #pragma omp parallel
    #pragma omp single
    {
    for (auto & [seq, id, qual] : query_file_h) {
        #pragma omp task firstprivate(seq, id, qual)
        {
        std::string s = id;
        std::regex rgx("(.+?) ID=(\\d+?)\\s(\\d+?)\\s(\\d+?)\\s([+,-])\\s(\\d+*)");
        std::smatch matches;
        std::string seq_id;
        unsigned int start, block_id, length;
        int strand;
        if(std::regex_search(s, matches, rgx)) {
            seq_id = matches[1].str();
            block_id = std::stoi(matches[2].str()) - 1;
            start = std::stoi(matches[3].str());
            strand = matches[5].str() == "-" ? -1 : 1;
            length = std::stoi(matches[6].str());
            if (strand == -1) {
                start = length - start + 1;
            }
        }
        if (query_block_vec[block_id].empty()) {
            #pragma omp critical 
            {
            if (query_block_vec[block_id].empty()) {
                query_block_vec[block_id].resize(seq.size());
                for (auto i = 0; i < seq.size(); ++i) {
                    query_block_vec[block_id][i] = std::make_shared<std::set<seqloc>>();
                }
            }
            }
        }
        unsigned int ng = 0;
        for (auto col_idx = 0; col_idx < seq.size(); ++col_idx) {
            if (seq[col_idx] == seqan3::gap{})
                continue;
            ++ng;
            seqloc sl = std::make_pair(seq_to_idx[seq_id], start + ng*strand);
            #pragma omp critical 
            {
            query_block_vec[block_id][col_idx]->insert(sl);
            }
        }
        #pragma omp critical 
        {
        seq_count += 1;
        bar.progress(seq_count, query_seq_count);
        }
        }
    }
    }
    bar.finish();
    bar.reset();


    // Compute intersections
    spdlog::debug("Computing metrics...");
    idx = 0;
    unsigned int fp = 0, fn = 0, tp = 0;
    #pragma omp parallel for reduction(+:fp, fn, tp)
    for (auto block_group : query_block_vec) {
        for (auto query_group : block_group) {
            for (seqloc sl : *query_group) {
                auto truth_group = loc_to_idx[sl];
                std::vector<seqloc> v_intersection;
                std::set_intersection(
                        truth_group->begin(), truth_group->end(), 
                        query_group->begin(), query_group->end(), 
                        std::back_inserter(v_intersection));
                tp += v_intersection.size();
                v_intersection.clear();
                std::set_difference(
                        truth_group->begin(), truth_group->end(), 
                        query_group->begin(), query_group->end(), 
                        std::back_inserter(v_intersection));
                fn += v_intersection.size();
                v_intersection.clear();
                std::set_difference(
                        query_group->begin(), query_group->end(), 
                        truth_group->begin(), truth_group->end(), 
                        std::back_inserter(v_intersection));
                fp += v_intersection.size();
                v_intersection.clear();
            }
        }
#pragma omp critical 
        {
        bar.progress(idx, query_block_vec.size());
        idx += 1;
}
    }
    bar.finish();
    bar.reset();
    spdlog::info("TP = {}, FN = {}, FP = {}", tp, fn, fp);
}
