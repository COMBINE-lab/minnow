#include <cstdio>
#include <iostream>
#include <random>
#include <unistd.h>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <cstring>

#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "FASTAParser.hpp"
#include "ProgOpts.hpp"
#include "Transcript.hpp"
#include "MinnowUtil.hpp"
#include "macros.hpp"

FASTAParser::FASTAParser() {}
FASTAParser::FASTAParser(const std::string& fname) : fname_(fname) {}

void FASTAParser::updateGeneLevelIntron(
    std::vector<Transcript>& transcripts ,
    std::unordered_map<std::string, uint32_t>& geneMap,
    std::unordered_map<uint32_t, uint32_t>& gene2LastTrMap,
    std::unordered_set<uint32_t>& geneIntronMap
){
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  std::unordered_map<std::string, size_t> nameToID;
  std::string sepStr = "::" ;


  std::vector<std::string> readFiles{fname_};
  size_t maxReadGroup{1000}; // Number of files to read simultaneously
  // size_t concurrentFile{1};  // Number of reads in each "job"

  single_parser parser(readFiles, 1, 1, maxReadGroup);
  parser.start();

  auto rg = parser.getReadGroup() ;
  while(parser.refill(rg)){
    for(auto& read : rg){
      std::string& header = read.name ;
      std::string geneName = header.substr(0, header.find_first_of(sepStr)) ;
      auto it = geneMap.find(geneName) ;
      if(it != geneMap.end()){
        auto geneId = it->second ;
        // Get the last transcript
        auto it2 = gene2LastTrMap.find(geneId) ;
        // If not present then search for the insert 
        // last transcript 
        uint32_t lastTr{0} ;
        if(it2 != gene2LastTrMap.end()){
          lastTr = it2->second ; 
        }
        else{
          std::cerr << "ERROR !! in FASTAParser::updateGeneLevelIntron -> Line 56\n" ;
          std::exit(1) ;
        }
        geneIntronMap.insert(geneId) ;
        auto& tr = transcripts[lastTr] ;

        std::string& seq = read.seq; 
        tr.insertPolyASeq(seq) ;
      }
    }
  }

  parser.stop() ;
}


void FASTAParser::updateTranscriptLevelIntron(
  std::vector<Transcript>& transcripts ,
  std::unordered_map<std::string, uint32_t>& transcriptNameMap
){
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using std::string;
  using std::unordered_map;

  unordered_map<string, size_t> nameToID;
  std::string sepStr = "::" ;

  std::vector<std::string> readFiles{fname_};
  size_t maxReadGroup{1000}; // Number of files to read simultaneously
  // size_t concurrentFile{1};  // Number of reads in each "job"

  single_parser parser(readFiles, 1, 1, maxReadGroup);
  parser.start();

  // constexpr char bases[] = {'A', 'C', 'G', 'T'};
  // Create a random uniform distribution
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<> dis(0, 3);
  // uint64_t numNucleotidesReplaced{0};

  // All header names we encounter in the fasta file
  std::unordered_set<std::string> fastaNames;
  // size_t id{0} ;

  auto rg = parser.getReadGroup();
  while (parser.refill(rg)){
    for(auto& read : rg){
      std::string& header = read.name ;
      std::string concatTranscriptNames = header.substr(0, header.find_first_of(sepStr)) ;
      std::vector<std::string> transcriptNames ; 
      
      uint32_t chrEndPos = std::stoul(header.substr(header.find_first_of("-")+1));

      util::split(concatTranscriptNames, 
                  transcriptNames,
                  ",") ;

      for(auto& name: transcriptNames){

        auto& tr = transcripts[transcriptNameMap[name]] ;
        std::string& seq = read.seq;
        //size_t readLen = seq.length();

        //Replace non-ACGT bases
        //for (size_t b = 0; b < readLen; ++b) {
        //  seq[b] = ::toupper(seq[b]);
        //  int c = jellyfish::mer_dna::code(seq[b]);
        //  // Replace non-ACGT bases with pseudo-random bases
        //  if (jellyfish::mer_dna::not_dna(c)) {
        //    char rbase = bases[dis(eng)];
        //    c = jellyfish::mer_dna::code(rbase);
        //    seq[b] = rbase;
        //    ++numNucleotidesReplaced;
        //  }
        //}

        tr.insertIntronSeq(chrEndPos, seq) ;
      }

    }
  }

  parser.stop() ;

}


void FASTAParser::populateTargets(
  std::vector<Transcript>& refs,
  uint32_t readLength
) {
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using std::string;
  using std::unordered_map;

  unordered_map<string, size_t> nameToID;
  
  // Separators for the header (default ' ' and '\t')
  // If we have the gencode flag, then add '|'.
  // TODO: make an option
  bool gencodeRef{true} ;

  std::string sepStr = " \t";
  if (gencodeRef) {
    sepStr += '|';
  }

  std::vector<std::string> readFiles{fname_};
  size_t maxReadGroup{1000}; // Number of files to read simultaneously
  //size_t concurrentFile{1};  // Number of reads in each "job"

  single_parser parser(readFiles, 1, 1, maxReadGroup);
  parser.start();

  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  // Create a random uniform distribution
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<> dis(0, 3);
  uint64_t numNucleotidesReplaced{0};

  // All header names we encounter in the fasta file
  std::unordered_set<std::string> fastaNames;
  size_t id{0} ;

  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    for (auto& read : rg) {
      std::string& header = read.name;
      std::string name = header.substr(0, header.find_first_of(sepStr));
      fastaNames.insert(name);
      
        // std::string& seq = j->data[i].seq;
        std::string& seq = read.seq;
        size_t readLen = seq.length();
        if(readLen < readLength)
          continue ;

        // Replace non-ACGT bases
        for (size_t b = 0; b < readLen; ++b) {
          seq[b] = ::toupper(seq[b]);
          int c = jellyfish::mer_dna::code(seq[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            seq[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // store transcript information
        refs.emplace_back(
          id++,
          name,
          readLen
        );
        // allocate space for the new copy
        char* seqCopy = new char[seq.length() + 1];
        std::strcpy(seqCopy, seq.c_str());
        refs[id-1].setSequenceOwned(seqCopy);
        // seqCopy will only be freed when the transcript is destructed!
    }
    
  }

  parser.stop();

  if(numNucleotidesReplaced > 0)
    std::cout<<
      "replaced " << numNucleotidesReplaced << " non-ACGT nucleotides with random nucleotides";
}

void FASTAParser::populateIntronTargets(
  std::vector<Transcript>& refs,
  std::string& intronFileName,
  std::unordered_map<std::string, uint32_t>& transcriptNameMap,
  uint32_t readLength
) {
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using std::string;
  using std::unordered_map;

  unordered_map<string, size_t> nameToID;
  
  // Separators for the header (default ' ' and '\t')
  // If we have the gencode flag, then add '|'.
  // TODO: make an option
  bool gencodeRef{true} ;

  std::string sepStr = " \t";
  if (gencodeRef) {
    sepStr += '|';
  }

  std::vector<std::string> readFiles{intronFileName};
  size_t maxReadGroup{1000}; // Number of files to read simultaneously
  //size_t concurrentFile{1};  // Number of reads in each "job"

  single_parser parser(readFiles, 1, 1, maxReadGroup);
  parser.start();

  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  // Create a random uniform distribution
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<> dis(0, 3);
  uint64_t numNucleotidesReplaced{0};

  // All header names we encounter in the fasta file
  std::unordered_set<std::string> fastaNames;

  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    for (auto& read : rg) {
      std::string& header = read.name;
      std::string name = header.substr(0, header.find_first_of(sepStr));
      fastaNames.insert(name);
      
        // std::string& seq = j->data[i].seq;
        std::string& seq = read.seq;
        size_t readLen = seq.length();
        if(readLen < readLength)
          continue ;

        // Replace non-ACGT bases
        for (size_t b = 0; b < readLen; ++b) {
          seq[b] = ::toupper(seq[b]);
          int c = jellyfish::mer_dna::code(seq[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            seq[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // find out the sequence id for this transcript
        auto it = transcriptNameMap.find(name);
        if(it != transcriptNameMap.end()){
          auto id  = it->second ;
          // allocate space for the new copy
          char* seqCopy = new char[seq.length() + 1];
          std::strcpy(seqCopy, seq.c_str());
          refs[id].setIntronSequenceOwned(seqCopy);
        }
        // seqCopy will only be freed when the transcript is destructed!
    }
    
  }

  parser.stop();

  if(numNucleotidesReplaced > 0)
    std::cout<<
      "replaced " << numNucleotidesReplaced << " non-ACGT nucleotides with random nucleotides";
}
