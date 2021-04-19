#ifndef __MINNOW_PROG_OPTS_HPP__
#define __MINNOW_PROG_OPTS_HPP__
#include <cstdint>
#include <string>

enum class BarcodeEnd { FIVE = 5, THREE = 3 };
enum class Sequence { BARCODE, UMI };
enum class ProtocolType {
                         CHROMIUMV3 = 0,
                         CHROMIUM = 1,
                         DROPSEQ = 2,
                         CELSEQ = 3,
                         CELSEQ2 = 4,
                         CUSTOM = 5
};

class IndexOptions {
public:
  uint32_t k{31};
  // k-mer to build the de-Bruijn graph
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

  // ===============================================
  // common options shared by both modes
  // ===============================================
  std::string inputdir ;
  std::string outDir ;
  std::string protocol{""} ;
  std::string whitelistFile ;
  uint32_t numOfPCRCycles{5};
  double errorRate{0.01};
  size_t librarySize{100000};
  double mutationProb{0.01};
  std::string refFile ;
  uint32_t sampleCells{0} ;
  uint32_t numOfCells{10};
  uint32_t numThreads{2};
  bool gencode{false} ;
  bool useDBG{false} ;
  bool switchOnEffModel{false} ;
  // PCR model
  uint32_t CBLength{16} ;
  uint32_t UMILength{10} ;
  uint32_t ReadLength{100} ;
  std::string gene2txpFile{""} ;
  std::string gfaFile{""} ;
  std::string clusterFile{""} ;
  std::string bfhFile{""} ;
  std::string rspdFile{""} ;
  bool noDump{false} ;
  // intended for not dumping the output files
  std::string illuminaModelFile{""} ;
  bool generateNoisyCells{false} ;

  // ===============================================
  // options for running from alevinMode
  // ===============================================
  bool alevinMode{false} ;
  bool dupCounts{false} ;
  bool useWhiteList{false} ;
  bool useEqClass{false} ;

  // ===============================================
  // options for running from splatterMode
  // ===============================================
  bool customNames{false} ;
  bool splatterMode{false} ;
  bool normalMode{false} ;
  bool testUniqness{false} ;
  bool reverseUniqness{false} ;
  bool useWeibull{false} ;
  std::string uniquenessFile{""} ;
  std::string geneProbFile{""} ;
  std::string countProbFile{""} ;
  std::string numMolFile{""} ;
  std::string metadataDir{""} ;
  std::string duplicateFile{""};

  // ================================================
  // currently not used in main branch
  // ================================================
  bool binary{false} ;
  // intended for reading the binary mode of the file
  std::string modelDir{""} ;
  // intended for storing all model files for analysis
  // Velocity mode
  bool velocityMode{false} ;
  // to switch on the velocity mode
  std::string intronFile{""} ;
  // splicing rate vector
  std::string spliceVecFile{""};

  // intended for intron retention
  std::string exonLengthFile{""} ;
  // exon lengths are stored to help simulate intron
  // retained reads
  bool samplePolyA{false} ;
  // Sample polyA tails
  std::string polyAsiteFile{""} ;
  // File that stores polyA site files
  std::string polyAsiteFractionFile{""} ;
  // Fraction of polyA site again needed for intron retention
  std::string matrixFile;
  // We don't use it now as input directory is changed
  std::string genomefile{""} ;
  // We needed genome file for intron sequences
  uint32_t numOfTranscripts{100};
  // In case this was a transcript level matrix
  uint32_t numOfDoublets{0};
  // The number of doublets to be produced

};

class EstimateOptions {
  // To estimate we need the exon junctions taken
  // from a gtf file and gene name
  public:
    std::string gene2txpFile{""} ;
    std::string eqClassFolder{""} ;
    std::string refFile ;
    uint32_t ReadLength{100};

    std::string bfhFile{""} ;
    std::string outDir{""} ;
    std::string clusterFile{""} ;

};


#endif
