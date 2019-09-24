#ifndef __MINNOW_PROG_OPTS_HPP__
#define __MINNOW_PROG_OPTS_HPP__
#include <cstdint>
#include <string>

class IndexOptions {
public:
  uint32_t k{31};
  uint32_t p{16};
  std::string gfa_file;
  std::string cfile;
  std::vector<std::string> rfile;
  std::string outdir;
  std::string header_sep{""};
  bool isSparse{false};
  bool keep_duplicates{false};
  uint32_t lossy_rate{5};
  int32_t filt_size{-1};
  std::string twopaco_tmp_dir{""};
};


class SimulateOptions {
public:
  bool alevinMode{false} ;

  bool splatterMode{false} ;
  bool normalMode{false} ;
  bool testUniqness{false} ;
  bool reverseUniqness{false} ;


  bool velocityMode{false} ;
  
  bool useDBG{false} ;
  
  bool binary{false} ;
  bool gencode{false} ;
  bool dupCounts{false} ;
  bool useWhiteList{false} ;
  bool generateNoisyCells{false} ;
  bool useEqClass{false} ;
  bool noDump{false} ;
  bool useWeibull{false} ;
  bool switchOnEffModel{false} ;

  bool samplePolyA{false} ;

  std::string clusterFile{""} ;
  std::string bfhFile{""} ;
  
  std::string modelDir{""} ;

  std::string intronFile{""} ;
  std::string exonLengthFile{""} ;
  std::string polyAsiteFile{""} ;
  std::string polyAsiteFractionFile{""} ;
  std::string matrixFile;
  std::string gene2txpFile{""} ;
  std::string gfaFile{""} ;
  std::string genomefile{""} ;
  std::string rspdFile{""} ;
  std::string geneProbFile{""} ;
  std::string countProbFile{""} ;
  std::string numMolFile{""} ;
  std::string uniquenessFile{""} ;
  std::string illuminaModelFile{""} ;

  std::string refFile ;
  std::string outDir ;
  uint32_t sampleCells{0} ;
  uint32_t numOfCells{10};
  uint32_t numOfTranscripts{100};
  uint32_t numOfPCRCycles{5};
  double errorRate{0.01};
  double mutationProb{0.01};
  uint32_t numThreads{2};
  uint32_t numOfDoublets{0};

  uint32_t CBLength{16} ;
  uint32_t UMILength{10} ;
  uint32_t ReadLength{100} ;
};

class EstimateOptions {
  // To estimate we need the exon junctions taken 
  // from a gtf file and gene name
  public:
    std::string gene2txpFile{""} ;
    std::string eqClassFolder{""} ;

    std::string bfhFile{""} ;
    std::string outDir{""} ;
    std::string clusterFile{""} ;
    
};


#endif 
