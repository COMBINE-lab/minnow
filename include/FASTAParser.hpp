#ifndef FASTA_PARSER
#define FASTA_PARSER

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include "Transcript.hpp"
#include "ProgOpts.hpp"

class FASTAParser {
public:
  FASTAParser();
  FASTAParser(const std::string& fname);
  void populateTargets(std::vector<Transcript>& transcripts);
  void populateIntronTargets(
    std::vector<Transcript>& refs,
    std::string& intronFileName,
    std::unordered_map<std::string, uint32_t>& transcriptNameMap
  ) ;
  void updateTranscriptLevelIntron(std::vector<Transcript>& transcripts, 
    std::unordered_map<std::string, uint32_t>& transcriptNameMap
  ) ;
  void updateGeneLevelIntron(
    std::vector<Transcript>& trancripts,
    std::unordered_map<std::string, uint32_t>& geneMap,
    std::unordered_map<uint32_t, uint32_t>& gene2LastTrMap,
    std::unordered_set<uint32_t>& geneIntronMap
  ) ;

private:
  std::string fname_;
};

#endif // FASTA_PARSER
