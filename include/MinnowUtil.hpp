#ifndef MINNOW_UTIL_HPP
#define MINNOW_UTIL_HPP


#include <cstdio>
#include <cmath>
#include <random>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <unistd.h>
#include "MinnowFS.hpp"

#include "string_view.hpp"
#include "macros.hpp"



namespace util{

  struct EditBlock{
    int editedInd;
    int editChar;
    bool active;
  };

  struct SequenceBlock{
    std::string sequence ;
    int count ;
  };

  struct TrRelPos{
        uint32_t eqId ;
        uint32_t start ;
        uint32_t end ;
        uint32_t end2 ;

    };

  //Main structural blocks that
  //makes up minnow and therefore structure
  //of the entire tree
  struct CellBarcodeUMIBasicBlock{
    int transcriptId; // trancriptId

    std::string cellBarCode ;
    std::string UMICode ;

    bool count ;

    uint32_t fragmentStart ;
    uint32_t fragmentEnd ;

    bool junction ;

    // next childBlock to generate in the next cycle
    // std::vector<CellBarcodeUMIChildBlock>::iterator childBlock ;
    // std::vector<EditBlock> childBlockChain ;
    // std::vector<EditBlock>::iterator childBlockChain;

    CellBarcodeUMIBasicBlock(
      int transcriptIdIn,
      std::string cellBarCodeIn,
      std::string UMICodeIn,
      bool countIn,
      uint32_t fragmentStartIn,
      uint32_t fragmentEndIn,
      bool junctionIn
    ) :
    transcriptId(transcriptIdIn),
    cellBarCode(std::move(cellBarCodeIn)),
    UMICode(std::move(UMICodeIn)),
    count(countIn),
    fragmentStart(fragmentStartIn),
    fragmentEnd(fragmentEndIn),
    junction(junctionIn)
    {}
  };

  // Made in the image of the
  // previous structure

  struct CellBarcodeUMISegmentBlock{
    uint32_t transcriptId; // trancriptId
    size_t segmentId ; // segment id
    uint32_t segmentStart ;
    uint32_t segmentEnd ;
    std::string cellBarCode ;
    std::string UMICode ;
    bool ore ;

    CellBarcodeUMISegmentBlock(
      uint32_t transcriptIdIn,
      uint32_t segmentIdIn,
      uint32_t segmentStartIn ,
      uint32_t segmentEndIn ,
      std::string cellBarCodeIn,
      std::string UMICodeIn,
      bool oreIn
    ) :
    transcriptId(transcriptIdIn),
    segmentId(segmentIdIn),
    segmentStart(segmentStartIn),
    segmentEnd(segmentEndIn),
    cellBarCode(std::move(cellBarCodeIn)),
    UMICode(std::move(UMICodeIn)),
    ore(oreIn)
    {}
  };


  struct CellBarcodeUMIChildBlock{
    std::string sequence ;
    int count ;

    CellBarcodeUMIChildBlock(
      std::string sequenceIn,
      int countIn
    ) :
    sequence(sequenceIn),
    count(countIn)
    {}

    inline void incrementCount() { count++ ; }
  };

  inline void split(const std::string& str, std::vector<std::string>& tokens, const std::string& delim){
    const std::string whiteSpace = " " ;
    size_t prev = 0, pos = 0;
    do
      {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty() and token != whiteSpace) tokens.push_back(token);
        prev = pos + delim.length();
      }
    while (pos < str.length() && prev < str.length());
  }

  inline std::string genRandomSeq(uint32_t length){
    std::vector<char> nuclVec = {'A', 'T', 'G', 'C'} ;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 3); // define the range

    std::string seq(length, 'A') ;
    for(size_t i = 0; i < length; ++i){
      seq[i] = nuclVec[distr(eng)] ;
    }
    return seq ;
  }


  template<typename T>
  inline T generateRandomNumber(T start, T end){
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<T> distr(start, end); // define the range

    return distr(eng);
  }

 // Follows the different start position distribution
  inline uint32_t generateStartPosition(){
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<double> distr(4.0, 6.0); // define the range

    return static_cast<uint32_t>(std::exp(distr(eng)));
  }

  inline std::vector<std::string> generateUMIList(){
    std::vector<std::string> whiteList ;
    std::cerr << "PRINTING DEBUG: POOL_SIZE " << POOL_SIZE << "\n\n" ;
    whiteList.resize(POOL_SIZE) ;
    for(int i=0 ; i < POOL_SIZE; ++i){
      whiteList[i] = genRandomSeq(UMI_LENGTH) ;
    }
    return whiteList ;
  }

  inline std::vector<std::string> generateCBList(int numCells){
    std::vector<std::string> CBList;
    CBList.resize(numCells) ;
    for(int i=0 ; i < numCells; ++i){
      CBList[i] = genRandomSeq(CB_LENGTH) ;
    }
    return CBList ;
  }


  inline std::vector<std::string> generate10XCBList(int numCells){


    //std::cout << "reading 10X specific list\n" ;
    std::vector<int> indexArray(TENX_WHITELIST_LENGTH) ;
    std::iota(indexArray.begin(), indexArray.end(), 0) ;

    std::vector<std::string> CBList ;
    CBList.resize(numCells) ;

    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 eng(seed);


    std::shuffle(indexArray.begin(), indexArray.end(), eng) ;

    std::vector<int> lines_to_select(indexArray.begin(), indexArray.begin() + numCells ) ;
    std::sort(lines_to_select.begin(), lines_to_select.end()) ;
    
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("/"));
    std::string dir_path_root = dir_path.substr(0, dir_path.rfind("/"));
    std::string whitelist_filename = dir_path_root + "/data/737K-august-2016.txt" ;

    std::cerr << "10X whitelist file " << whitelist_filename << "\n";
    
    //std::string whitelist_filename("/mnt/scratch1/hirak/minnow/data/737K-august-2016.txt") ;

    std::string line;
    if(! util::fs::FileExists(whitelist_filename.c_str())){
			std::cerr << "10X whitelist file does not exist\n" ;
			std::exit(1) ;
		}

    //std::cout << "here \n" ;
    std::ifstream file(whitelist_filename.c_str()) ;

    int line_number = 0;
    int list_number = 0;

    while (std::getline(file, line)) {
        if (line_number == lines_to_select[list_number]) {
          CBList[list_number] = line ;
          ++list_number ;
        }
      ++line_number;
    }

    return CBList ;
  }

   inline char complementt(char& c) {
    switch (c) {
    case 'A':
      c = 'T';
      return c;
    case 'T':
      c = 'A';
      return c;
    case 'C':
      c = 'G';
      return c;
    case 'G':
      c = 'C';
      return c;
    }
    return 'N';
  }

  inline std::string revcomp(std::string s) {
    int n = s.size();
    int halfLength = s.size() / 2;
    for (int i = 0; i < halfLength; i++) {
      char temp = complementt(s[i]);
      s[i] = complementt(s[n - 1 - i]);
      s[n - 1 - i] = temp;
    }
    if (s.size() % 2 != 0) {
      s[halfLength] = complementt(s[halfLength]);
    }
    return s;
  }

  inline uint32_t calculateParent(size_t index, size_t blockSize){
    auto bucketId = index / blockSize ;
    uint32_t prevPcrCycleNum = std::floor(std::log2(bucketId))  ;

    return (index % (blockSize * static_cast<size_t>(std::pow(2, prevPcrCycleNum)))) ;
  }


  inline void readIlluminaErrorModel(
    std::string& errorProfileFile, 
    std::vector<std::vector<std::vector<double>>>& errorModel
  ){

    if(! util::fs::FileExists(errorProfileFile.c_str())){
      std::cerr << "Error model file " << errorProfileFile << " does not exist \n" ;
      std::exit(1) ;
    }
    
    errorModel.resize(101) ;
    for (int i = 0; i < 101; i++){
        errorModel[i].resize(5);
        for (int j = 0; j < 5; j++)
        {
          errorModel[i][j].resize(5);
        }
    }


    std::string line;
    std::ifstream modelStream{errorProfileFile.c_str()} ;

    size_t lineCount = 0 ;
    size_t nuclCount = 0 ;

    while(std::getline(modelStream, line)){
      if(lineCount == 0){
        lineCount++ ;
        continue ;
      }

      std::vector<std::string> cellCountsString ; 
      util::split(line, cellCountsString, "\t") ;
      int pos = std::stoi(cellCountsString[6]) ;

      for(size_t fid = 1 ; fid < 6 ; ++fid){

        errorModel[pos][nuclCount%5][fid - 1] = std::stod(cellCountsString[fid]) ; 
      
      }
      nuclCount += 1 ;
    }
  }

}



#endif
