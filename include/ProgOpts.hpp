#ifndef __MINNOW_PROG_OPTS_HPP__
#define __MINNOW_PROG_OPTS_HPP__
#include <cstdint>
#include <string>

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
