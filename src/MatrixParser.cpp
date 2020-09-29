
#include <set>
#include <memory>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "MatrixParser.hpp"
#include "MinnowFS.hpp"
#include "BFHClass.hpp"
#include "GFAReader.hpp"

#include <unordered_set>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <cassert>

#include "zstr.hpp"


#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)

void loadTSVFile(
                 std::string& tsvFile,
                 std::unordered_map<std::string, uint32_t>& geneMap,
                 std::vector<double>& fractionVector
                 ){
    if(! util::fs::FileExists(tsvFile.c_str())){
        std::cerr << "the tsv file does not exist\n" ;
        std::exit(1) ;
    }

    fractionVector.resize(geneMap.size(), 0.0) ;
    std::ifstream dataStream(tsvFile.c_str()) ;
    std::string line ;

    // fill the fraction vector
    while(std::getline(dataStream, line)){
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        std::vector<std::string> tokens ;
        util::split(line, tokens, "\t") ;
        if(tokens.size() == 2){
            auto it = geneMap.find(tokens[0]) ;
            if(it != geneMap.end()){
                fractionVector[it->second] = std::stod(tokens[1]) ;
            }
        }
    }
}

template <typename T>
void pupulateGeneCountMatrix(
    std::string& countFile,
    std::vector<std::vector<T>>& geneCount,
    uint32_t& numCells,
    size_t& numOfGenes
){
    if(! util::fs::FileExists(countFile.c_str())){
        std::cerr << "quants_mat.csv file does not exist\n" ;
        std::exit(1) ;
    }
    std::ifstream dataStream(countFile.c_str()) ;
    std::string line ;

    size_t cellId{0} ;
    while(std::getline(dataStream, line)){
        cellId += 1 ;

        std::vector<std::string> cellGeneCountsString ;
        util::split(line, cellGeneCountsString, ",") ;

        // should not have missing count
        assert(cellGeneCountsString.size() == numOfGenes) ;

        std::transform(cellGeneCountsString.begin(), cellGeneCountsString.end(), geneCount[cellId].begin(), [](const std::string& val)
        {
            return std::stod(val);
        });
        if (cellId == numCells){
            break ;
        }
        _verbose("\rNumber of cells read genes : %lu", cellId);
    }

}


template<typename T>
void populateGeneCountMatrix(
  std::string& countMatFilename,
  std::vector<std::vector<T>>& geneCount,
  uint32_t& numCells,
  size_t& numOfGenes,
  std::unordered_map<uint32_t, uint32_t>& original2whitelistMap,
  std::unordered_map<uint32_t, uint32_t>& original2NoisyMap,
  uint32_t& numOfOriginalCells
){


      if(! util::fs::FileExists(countMatFilename.c_str())){
        std::cerr << "quants_mat.gz file does not exist\n" ;
        std::exit(1) ;
      }

      bool predefinedCells = (original2whitelistMap.size() != 0) ;
      auto popcount = [](uint8_t n) {
        size_t count {0};
        while (n) {
          n &= n-1;
          ++count;
        }
        return count;
      };

      uint32_t zerod_cells {0};
      size_t numFlags = std::ceil(numOfGenes/8.0);
      std::vector<uint8_t> alphasFlag (numFlags, 0);
      size_t flagSize = sizeof(decltype(alphasFlag)::value_type);

      std::vector<float> alphasSparse;
      alphasSparse.reserve(numFlags/2);
      size_t elSize = sizeof(decltype(alphasSparse)::value_type);

      std::unique_ptr<std::istream> in =
        std::unique_ptr<std::istream> (
             new zstr::ifstream(countMatFilename.c_str(),
             std::ios::in | std::ios::binary)
        ) ;

      //std::unique_ptr<std::istream> in =
      //  std::unique_ptr<std::istream> (
      //                                 new std::ifstream(countMatFilename.c_str(),
      //                                 std::ios::in | std::ios::binary)
      //                                 ) ;
      // Loop over the binary matrix

      size_t cellCount{0} ;
      for (size_t cellId = 0 ; cellId < numOfOriginalCells ; ++cellId){

        // read the row in alphasSparse
        in->read(reinterpret_cast<char*>(alphasFlag.data()), flagSize * numFlags);
        size_t numExpGenes {0};

        std::vector<size_t> indices;
        for (size_t j=0; j < alphasFlag.size(); j++) {
          uint8_t flag = alphasFlag[j];
          size_t numNonZeros = popcount(flag);
          numExpGenes += numNonZeros;


          for (size_t i=0; i<8; i++){
            if (flag & (128 >> i)) {
              indices.emplace_back( i+(8*j) );
            }
          }
        }

        if (indices.size() != numExpGenes) {
          std::cerr<< "binary format reading error "  << indices.size() << numExpGenes << "\n";
          exit(84);
        }
        alphasSparse.clear();
        alphasSparse.resize(numExpGenes, 0.0);

        in->read(reinterpret_cast<char*>(alphasSparse.data()), elSize * numExpGenes);

        // the cell quant values are present in alphasSparse
        float readCount {0.0};
        readCount = std::accumulate(alphasSparse.begin(), alphasSparse.end(), 0.0);

        if(predefinedCells){
          if (original2whitelistMap.find(cellId) != original2whitelistMap.end()){
            uint32_t whitelistCellId = original2whitelistMap[cellId] ;
            if(whitelistCellId >= numCells)
              continue ;
            for(size_t i = 0 ; i < numExpGenes ; i++){
              //std::cerr << alphasSparse[i] << "\n" ;
              geneCount[whitelistCellId][indices[i]] = alphasSparse[i] ;
            }
            cellCount++ ;
            _verbose("\rNumber of cells read  : %lu", cellCount);

            if(cellCount == numCells){
              break ;
            }

          }else if(original2NoisyMap.find(cellId) != original2NoisyMap.end()){
            uint32_t noisyCellId = original2NoisyMap[cellId] ;

            for(size_t i = 0 ; i < numExpGenes ; i++){
              //std::cerr << alphasSparse[i] << "\n";
              geneCount[noisyCellId][indices[i]] = alphasSparse[i] ;
            }
            cellCount++ ;
            _verbose("\rNumber of cells read noisy  : %lu", cellCount);
            if (cellCount == numCells){
              break ;
            }
          }
        }else{

          cellCount++ ;

          for(size_t i = 0 ; i < numExpGenes ; i++){
            //std::cerr << alphasSparse[i] << "\n";
            geneCount[cellId][indices[i]] = alphasSparse[i] ;
          }

          _verbose("\rNumber of cells read  : %lu", cellCount);
          if (cellCount == numCells){
            break ;
          }

        }


        if (readCount == 0.0){
          zerod_cells += 1;
        }

      } // end-for each cell

      if (zerod_cells > 0) {
        std::cerr << "Found {} cells with 0 counts " <<  zerod_cells << "\n";
      }


      double totSum{0} ;
      for(size_t i = 0 ; i < geneCount.size() ; ++i){
        for(size_t j = 0 ; j < geneCount[0].size() ; ++j){
          totSum += geneCount[i][j] ;
        }
      }
      std::cerr << "\n\tTotal sum " << totSum << "\n" ;
}


template<typename T>
void populateGeneCountMatrix2(
    std::string& countFileBinary,
    std::vector<std::vector<T>>& geneCount,
    uint32_t& numCells,
    size_t& numOfGenes,
    std::unordered_map<uint32_t, uint32_t>& original2whitelistMap,
    std::unordered_map<uint32_t, uint32_t>& original2NoisyMap,
    uint32_t& numOfOriginalCells
)
{
    if(! util::fs::FileExists(countFileBinary.c_str())){
        std::cerr << "quants_mat.gz file does not exist\n" ;
        std::exit(1) ;
    }
    bool predefinedCells = (original2whitelistMap.size() != 0) ;
    std::cerr << "Num cells ::: " << numCells << "\n" ;
    std::cerr << "Gene count dim " << geneCount.size() << "\t" <<  geneCount[geneCount.size()-1].size() << "\n" ;
    std::cerr << "noisyMap size " << original2NoisyMap.size() << std::endl  ;
    //std::cerr << "\n DEBUG ==> reading form binary with numCells: " << numCells <<  " numGenes: " << numOfGenes << "\n"  ;
    //std::cerr << "\n DEBUG ==> size of geneCount "<< geneCount.size() << " " << geneCount[geneCount.size() - 1].size() << "\n" ;
    std::unique_ptr<std::istream> in =
        std::unique_ptr<std::istream> (new zstr::ifstream(countFileBinary.c_str(), std::ios::in | std::ios::binary)) ;

    size_t elSize = sizeof(typename std::vector<double>::value_type);
    size_t cellCount{0} ;
    if(predefinedCells){
        for(size_t cellId = 0; cellId < numOfOriginalCells; ++cellId){
            std::vector<T> cell(numOfGenes) ;
            in->read(reinterpret_cast<char*>(cell.data()), elSize * numOfGenes) ;
            if (original2whitelistMap.find(cellId) != original2whitelistMap.end()){
                uint32_t whitelistCellId = original2whitelistMap[cellId] ;
                // std::cerr <<" cellid: " << cellId << "\tWhitelistCellId: " << whitelistCellId << std::endl ;
                if(whitelistCellId >= numCells)
                    continue ;
                geneCount[whitelistCellId] = cell ;
                cellCount++ ;
                _verbose("\rNumber of cells read  : %lu", cellCount);
                if (cellCount == numCells){
                    break ;
                }
            }
            else if(original2NoisyMap.find(cellId) != original2NoisyMap.end()){
                uint32_t noisyCellId = original2NoisyMap[cellId] ;
                geneCount[noisyCellId] = cell ;
                cellCount++ ;
                _verbose("\rNumber of cells read noisy  : %lu", cellCount);
                if (cellCount == numCells){
                    break ;
                }
            }
        }
    }else{
        for(auto& cell : geneCount){
            cellCount += 1 ;
            //std::cerr << cellCount << "\n" ;
            _verbose("\rNumber of cells read  : %lu", cellCount);
            in->read(reinterpret_cast<char*>(cell.data()), elSize * numOfGenes) ;
            if (cellCount == numCells){
                break ;
            }
        }
        size_t numOfExpressedGenes{0} ;
        for(auto v : geneCount){
            for(auto c: v){
                if(c > 0){
                    numOfExpressedGenes += 1 ;
                }
            }
        }

    }

    assert(cellCount == numCells) ;

}

/*
**	This module loads data from alevin folder, it reads 
**	the matrix and other related information from the same 
**	directory in order to simulate reads 
*/

template<typename T>
void DataMatrix<T>::loadAlevinData(
    SimulateOptions& simOpts,
    Reference& refInfo
){

    // Load the values from simOpts
    // Basic Options
	auto alevinDir = simOpts.inputdir;
	auto sampleCells = simOpts.sampleCells;
	auto outDir = simOpts.outDir ;
	bool useDBG = (simOpts.gfaFile.size() != 0) ;
	auto gfaFile = simOpts.gfaFile ;
	auto bfhFile = simOpts.bfhFile ;
	bool fivePrime = true;

    // Advanced options
    bool samplePolyA = simOpts.samplePolyA;
    bool dupCounts = simOpts.dupCounts ;
    bool generateNoisyCells = simOpts.generateNoisyCells ;
    auto cellClusterFile = simOpts.clusterFile ;
    bool createDoublet = (simOpts.numOfDoublets == 0) ? false : true ;
    numOfDoublets = simOpts.numOfDoublets ;


    // load alevin related files
    if(! util::fs::DirExists(alevinDir.c_str())){
        consoleLog->error("Alevin directory does not exists") ;
        std::exit(1) ;
    }
    if(!simOpts.useWhiteList && generateNoisyCells){
        consoleLog->warn("--generateNoisyCells needs to be invoked in conjunction with --useWhiteList") ;
    }

    std::ifstream indata ;

    // Here we model the exact same duplicated counts for each
    // cell, this is proportional to the actually mapped the
    // reference. This is explicitely obtained fron alevin run

    if (dupCounts){
        consoleLog->info("Reading duplicated read numbers") ;
        bool alevin_updated{false} ;
        std::string dupCountFile_1 = alevinDir + "/MappedUmi.txt" ;
        std::string dupCountFile_2 = alevinDir + "/featureDump.txt" ;

        std::string dupCountFile ;

        if(! util::fs::FileExists(dupCountFile_2.c_str())){
        if(! util::fs::FileExists(dupCountFile_1.c_str())){
            consoleLog->error("Neither MappedUmi.txt nor featureDump.txt" 
                              "exists run without --dupCounts") ;
            std::exit(2) ;
        }else{
            dupCountFile = dupCountFile_1 ;
        }
        }else{
        alevin_updated = true ;
        dupCountFile = dupCountFile_2 ;
        }

        if(!alevin_updated){
            std::ifstream dupCountStream(dupCountFile) ;
            std::string line ;
            while(std::getline(dupCountStream, line)){
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                std::vector<std::string> tokens ;
                util::split(line, tokens, "\t") ;
                if (tokens.size() == 2){
                cellNamesDupCount[tokens[0]] = std::stoul(tokens[1]) ;
                }
            }
        }else{
            std::ifstream dupCountStream(dupCountFile) ;
            std::string line ;

            // throw away the first line
            std::getline(dupCountStream, line) ;

            while(std::getline(dupCountStream, line)){
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            std::vector<std::string> tokens ;
            util::split(line, tokens, "\t") ;
                if (tokens.size() > 2){
                    cellNamesDupCount[tokens[0]] = std::stoul(tokens[2]) ;
                }else{
                    consoleLog->error("{} is does not have enough columns", dupCountFile);
                    std::exit(1);
                }
            }
        }

        if(!dupCounts && (simOpts.numMolFile != "")){
            if(! util::fs::FileExists(simOpts.numMolFile.c_str())){
                std::cerr << simOpts.numMolFile << " does not exist\n" ;
            }else{
                std::ifstream dupCountStream(simOpts.numMolFile) ;
                std::string line ;
                while(std::getline(dupCountStream, line)){
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                std::vector<std::string> tokens ;
                util::split(line, tokens, "\t") ;
                if (tokens.size() == 2){
                    cellNamesDupCount[tokens[0]] = std::stoul(tokens[1]) ;
                }
                }

            }
        }
    }
      // End reading the duplicated counts


    // The map contains gid to tid map, where
    // the gid is created first time here as, we
    // get to know about them first here. The tids
    // are taken from the fasta file. Other transcripts
    // will be of no use as we don't know their
    // sequence.
    refInfo.updateGene2TxpMap() ;
    if(simOpts.duplicateFile.size()) { refInfo.updateDuplicateMap(simOpts.duplicateFile) ;}
    auto& geneMap = refInfo.geneMap ;
    consoleLog->info("Number of genes in the txp2gene file: {}", geneMap.size()) ;


    // This is an intermediate map that keeps a map
     // between gene index and the gene names from
    // the alevin produced file. This will be used to
    // map the index of the matrix to gene index in
    // the reference.

    std::unordered_map<uint32_t, std::string> alevinGeneIndex2NameMap ;

    // we might want to skip some genes
    std::unordered_set<uint32_t> skippedGenes ;
    std::unordered_set<std::string> skippedGeneNames ;
    uint32_t numOfSkippedGenes{0} ;
    size_t numOfOriginalGenes{0} ;


      consoleLog->info("====================Parsing Alevin Directory==========================") ;
      consoleLog->info("Start parsing Alevin Directory") ;
    consoleLog->info("Parsing {}/quants_mat_cols.txt",alevinDir) ;
    {
        std::string geneListFile = alevinDir + "/quants_mat_cols.txt" ;
        if(! util::fs::FileExists(geneListFile.c_str())){
            consoleLog->error("quants_mat_cols.txt file does not exist") ;
            std::exit(3) ;
        }
        std::ifstream geneNameStream(alevinDir + "/quants_mat_cols.txt") ;
        std::string line ;
        uint32_t geneIndex{0} ;
        while(std::getline(geneNameStream, line)){
            // strip new line
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            // search this gene name
            auto it = geneMap.find(line) ;
            if(it != geneMap.end()){
                alevin2refMap[geneIndex] = it->second ;
                alevinGeneIndex2NameMap[geneIndex] = line ;
                geneIndex += 1 ;
            }else{
                //consoleLog->error("quant_mat_cols.txt contains gene {} that is not present in the reference ",line)  ;
                //std::exit(1) ;
                skippedGeneNames.insert(line);
                skippedGenes.insert(numOfOriginalGenes) ;
                numOfSkippedGenes++ ;
            }
            numOfOriginalGenes++ ;
        }
    }

    //FIXME: One should not write the genes here. They can be truncated
    // std::string cellColFile = simOpts.outDir + "/alevin/quants_mat_cols.txt" ;
    // std::ofstream cellColStream(cellColFile.c_str()) ;
    // for(uint32_t i= 0; i < alevinGeneIndex2NameMap.size() ; ++i){
    // 	cellColStream << alevinGeneIndex2NameMap[i] << "\n" ;
    // }
    if (numOfSkippedGenes > 0){
        consoleLog->warn("Original number of genes: {}\tNumber of genes skipped: {}", numOfOriginalGenes, numOfSkippedGenes) ;
        consoleLog->warn("This means not all genes given in the input alevin matrix will be utilized");
        std::string skippedGenesFile = outDir + "/skipped_genes.txt" ;
        std::ofstream skippedGeneStream(skippedGenesFile.c_str()) ;
        
        for(auto skg : skippedGeneNames){
            skippedGeneStream << skg << "\n" ; 
        }

    }

    // consoleLog->info("Number of genes in the alevin produced files: {}",alevin2refMap.size()) ;

    auto& gene2transcriptMap = refInfo.gene2transcriptMap ;
    // alevin2refTranscriptMap, a map from columns of the
    // cell x transcript count matrix to be formed to the
    // transcript ids of the reference. This map is *very*
    // important since we need the reference id to get the
    // real sequences.

    std::map<uint32_t, uint32_t> alevinReverseMap ;
    {
        // create a transcript map
        uint32_t trId{0} ;
        for(auto it : alevin2refMap){
            auto gIt = gene2transcriptMap.find(it.second) ;
            if(gIt != gene2transcriptMap.end()){
                auto& trVec = gIt->second ;
                //std::cerr << "# of transcripts in this gene "<< gIt->first 
                //	      <<  ":\t" << trVec.size() << "\n" ;	
                for(auto tr : trVec){
                    alevin2refTranscriptMap[trId] = tr ;
                    alevinReverseMap[tr] = trId ;
                    trId += 1 ;
                }
            }else{
                consoleLog->error("This should not happen: gene {} not found",gIt->first) ;
                std::exit(5);
            }
        }
    }
      consoleLog->info("Parsing {}/quants_mat_rows.txt",alevinDir) ;
    // all cell names irrespective of whitelist
    std::map<std::string, uint32_t> allCellListMap ; // map cell-name -> id
    std::vector<std::string> allCellNames ; // {cell names}

    // track a cell by barcode
    std::string cellToTrack  = ""; //"GTATTCTGTAGTGAAT"; // "TGCGGGTCAGACGCAA" ;
    //uint32_t cellIdToTrack{0} ;
    // read cell names/barcodes that alevin produce
    {
        std::string cellListFile = alevinDir + "/quants_mat_rows.txt" ;
        if(! util::fs::FileExists(cellListFile.c_str())){
            consoleLog->error("quants_mat_rows.txt file does not exist") ;
            std::exit(1) ;
        }
        std::ifstream cellNameStream(alevinDir + "/quants_mat_rows.txt") ;
        std::string line ;
        while(std::getline(cellNameStream, line)){
            // strip new line
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            allCellNames.push_back(line) ;
            allCellListMap[line] = allCellNames.size() - 1 ;
        }
    }
    uint32_t numOfOriginalCells = allCellNames.size() ;
    consoleLog->info("Number of cells in the alevin produced files: {}",numOfOriginalCells) ;

    // Read the whitelist to know which cells are going to be finally part of the matrix.
    // If the whitelist is not present then treat the rows as whitelist
    std::string cellListFile = alevinDir + "/whitelist.txt" ;
    // When noisy cells are present and considered, then
    // allCells are the combination of whitelist cells and
    // noisy cells. These two maps would keep track of those
    // two groups and map original id of the cell list to
    // the whitelist and noisylist.

    // FIXME: If we sample from white-list, it might not be contiguous,
    // it is not contiguous create a second order mapping. in that case
    // we should first

    std::unordered_map<uint32_t, uint32_t> original2whitelistMap ;
    std::unordered_map<uint32_t, uint32_t> original2NoisyMap ;

    // Minnow would also keep a list of doublet cells, which will be
    // created later
    std::unordered_map<uint32_t, uint32_t> original2DoubletMap ;
    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> doubletContainer ;


    // Simulate noisy cell if asked by user.

    // If whitelist by alevin is not present we will treat all cells as whitelist
    // cells.
    // TODO: Given an in-built noise detector we can use to threshold
    // and call it noise that would be great. Not now.


    if((! util::fs::FileExists(cellListFile.c_str())) or !(simOpts.useWhiteList)){
        consoleLog->info("whitelist.txt file does not exist/ or will NOT be used");
        consoleLog->info("we need to assume that the rows are the whitelisted barcodes") ;
        // Copy the vector and map for now, and try something better
        // later
        cellNames = allCellNames ;
        cellWhiteListMap = allCellListMap ;
    }else{
        consoleLog->info("Parsing {}/whitelist.txt",alevinDir) ;
        {

            std::ifstream cellNameStream(alevinDir + "/whitelist.txt") ;
            std::string line ;
            while(std::getline(cellNameStream, line)){
                // strip new line
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                cellNames.push_back(line) ;
                cellWhiteListMap[line] = cellNames.size() - 1 ;
                // cellsToBeRead.insert(allCellListMap[line]) ;
                //original2whitelistMap[allCellListMap[line]] = cellWhiteListMap[line] ;
            }
        }
        // fill up the original2whitelistMap by going
        // over the cellWhiteListMap
        for(size_t i = 0 ; i < cellNames.size(); ++i){
            original2whitelistMap[allCellListMap[cellNames[i]]] = i ;
        }
        // NOTE: Noisy cells if neccessary will be appended
        // in the same vector of cell names. This makes
        // cellWhiteListMap deprecated for the case where
        // whitelist is provided. as original2whitelistMap and
        // allCellListMap together provide the same information
        
        // TODO: remove cellWhiteListMap in future release
        
        // NOTE: Not in use now
        consoleLog->info("Number of cells in whitelist file: {}", cellNames.size()) ;
        if(generateNoisyCells){
            consoleLog->info("Additionally reads from noisy cells will be generated too,"
                             "keeping track of noisy cells");

            if(numOfNoisyCells == 0){
                numOfNoisyCells = allCellNames.size() - cellNames.size() ;
            }
            uint32_t noisyCellId = cellNames.size() ;
            for(auto it : allCellListMap){
                if(cellWhiteListMap.find(it.first) == cellWhiteListMap.end()){
                    cellNames.push_back(it.first) ;
                    cellNoisyMap[it.first] = cellNames.size() - 1 ;
                    original2NoisyMap[it.second] =  cellNoisyMap[it.first];
                    noisyCellId++ ;
                    // TODO: insert a condition to break

                }
            }
            consoleLog->info("Additionally generate data from {} noisy experiments", noisyCellId) ;
        }
    }


    if ((sampleCells > 0) and (sampleCells < cellNames.size())){
        consoleLog->info("We will sample {} from {} cells",sampleCells, cellNames.size()) ;
    }else{
        sampleCells = cellNames.size() ;
    }

    numOfWhiteListedCells = cellWhiteListMap.size() ;
    numCells = sampleCells ;
    numOfTranscripts = alevin2refTranscriptMap.size() ;
    size_t numOfGenes = alevin2refMap.size() ;

    // We create a temporary gene count matrix which can serve the purpose of before
    // writing the counts to the data file
    std::vector<std::vector<T>> originalGeneCountMatrix ;

  // TODO: Depending of whether we use DBG we can get rid of either of
  // trueGeneCounts or data
  // NOTE: This matrix is to store the input of the alevin dataset directly
    originalGeneCountMatrix.resize(sampleCells, std::vector<T>(numOfOriginalGenes)) ;
  // NOTE: This matrix is when we need to have gene x transcripts
    data.resize(sampleCells + numOfDoublets, std::vector<T>(numOfTranscripts)) ;
  // NOTE: This matrix is made after multinomial sampling of the original matrix
    geneCounts.resize(sampleCells + numOfDoublets, std::vector<T>(numOfGenes)) ;
  // NOTE: This matrix is constructed in DBG mode
	trueGeneCounts.resize(sampleCells + numOfDoublets, std::vector<int>(numOfGenes)) ;

	// this matrix might have more genes
	for(auto& v : originalGeneCountMatrix)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : data)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : geneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : trueGeneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	std::string countFile = alevinDir + "/quants_mat.csv" ;
	std::string countFileBinary = alevinDir + "/quants_mat.gz" ;
  	bool binary{true} ;
  	if(! util::fs::FileExists(countFileBinary.c_str())){
  		binary = false ;
  	}


  	// TODO: Not sure if this is needed any more
	std::unordered_map<uint32_t, uint32_t> gene2LastTrMap ;
	{
		// Make a map from gene id to last last tid
		for(auto& gmpair : geneMap){
			auto trIds = refInfo.gene2transcriptMap[gmpair.second] ;
			auto lastTid = trIds[trIds.size() - 1] ;
			auto alevinLastTid = alevinReverseMap[lastTid] ;
			gene2LastTrMap[gmpair.second] = alevinLastTid ;
		}

	}

	// NOTE: This is important for intron retention
	// Decide whether to sample from the introns or not
	// depending on that use a flag for that cell and gene id
	// store a flag
	std::unordered_set<uint32_t> geneIntronMap ;
	if(samplePolyA){
		consoleLog->info("Sampling from polyA sites in between gene") ;
		auto& tsvFile = simOpts.polyAsiteFractionFile ; // fraction file 
		auto& polyAFastaFile = simOpts.polyAsiteFile ; // fasta file

		// load polyA fraction
		loadTSVFile(
			tsvFile,
			geneMap,
			fractionVector
		) ;
		// load the fasta file with the sequences 
		refInfo.updatePolyASequence(
			polyAFastaFile,
			geneMap,
			gene2LastTrMap,
			geneIntronMap
		) ;
	}


	bool useClusters{false} ;
	// TODO:
	if(useDBG){
		// This will read the dbg now
		if(!util::fs::FileExists(gfaFile.c_str())){
			consoleLog->error("GFA file {} does not exist EXITING !!", gfaFile);
			std::exit(4); 
		}

		auto rspdFile = simOpts.rspdFile ;
		if(rspdFile != ""){
			if(util::fs::FileExists(rspdFile.c_str())){
				rspdVec.resize(MAX_FRAGLENGTH) ;
				std::ifstream rspdStream(rspdFile.c_str()) ;
				std::string line ;
				while(std::getline(rspdStream, line)){
					std::vector<std::string> tokens ;
					util::split(line, tokens, ",") ;
					if(tokens.size() != 2)
						continue ;
					uint32_t pos = std::stoul(tokens[0]) ;
					uint32_t prob = std::stoul(tokens[1]) ;
					if(pos < rspdVec.size()){
						rspdVec[pos] = prob ;
					}
				}
			}else{
				consoleLog->info("RSPD file is not empty and doesn't exist, going with truncated sampling") ;
			}
		}

    	consoleLog->info("======================= Parsing GFA file {} ==========================",gfaFile) ;

		dbgPtr = new GFAReader(gfaFile, consoleLog) ;
		dbgPtr->parseFile(refInfo, outDir) ;


		// NOTE: stand alone call to bfh
		// ignore other eqclass stuff for now
	    // They might come handy later.
		// std::string eqFileDir = "dummy/dir" ;
    	consoleLog->info("======================= Parsing BFH/related file ==========================") ;
		eqClassPtr = new BFHClass(consoleLog) ;
		if((simOpts.countProbFile != "") || (simOpts.geneProbFile != "")){
			if(simOpts.countProbFile != ""){
				eqClassPtr->loadProbability(
					simOpts.countProbFile,
					refInfo,
					false
				) ;
			}else{
				eqClassPtr->loadProbability(
					simOpts.geneProbFile,
					refInfo,
					true
				) ;
			}
		}else if(bfhFile != ""){
      		// NOTE: This branch is currently
      		// active
			eqClassPtr->loadBFH(
				bfhFile,
				simOpts.clusterFile,
				refInfo,
				cellWhiteListMap,
				false,
				cellNoisyMap,
				simOpts.outDir
			) ;
		}else{
			// FIXME: Or load currently existing default file, this 
			// shoud be changed
			std::string geneProbFile = "../data/hg/geneLebelProb_pbmc_4k.txt" ;
			if(!util::fs::FileExists(geneProbFile.c_str())){
				consoleLog->error("alevin-mode is invoked with --dbg but neither bfh file"
								  " no probability files are produced"
				);
				std::exit(5) ;
			}
			eqClassPtr->loadProbability(
					geneProbFile,
					refInfo,
					true
			) ;
		}

		consoleLog->info("Parsed BFH/Prob file related information...") ;
		consoleLog->info("Genes in BFH: {}", numOfGenes) ;
		preCalculatedSegProb.resize(numOfGenes) ; // per gene
		preSegOreMapVector.resize(numOfGenes) ;
		geneSpecificTrVector.resize(numOfGenes) ;


		// fill segment based probability
		// gene to tr to prob
		// uint32_t trackGid{54819} ;
		for(uint32_t i = 0; i < numOfGenes; ++i){
			auto it = alevin2refMap.find(i) ;
			if(it != alevin2refMap.end()){
				auto originalGeneId = it->second ;
				auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

				std::unordered_map<size_t, uint32_t> localGeneProb ;
				std::unordered_map<size_t, bool> localSegOreMap ;
				std::unordered_map<size_t, std::vector<trInfo>> localTrVector ;

				for(auto tid : transcriptIds){
					refInfo.transcripts[tid].setGeneId(i);

					if(dbgPtr->trSegmentMap.find(tid) != dbgPtr->trSegmentMap.end()){
						auto segVec = dbgPtr->trSegmentMap[tid] ;
						bool everRC{false} ;
						for(auto seg : segVec){
							auto segCount = dbgPtr->eqClassMapCounts[seg] ;
							auto bfhCount = eqClassPtr->getGeneLevelProbCount(
								originalGeneId,
								segCount
							) ;

							if(bfhCount != 0){
								auto tcInfoVec = dbgPtr->eqClassMap[seg][tid] ;
								for(auto tInfo : tcInfoVec){
									if(
										(refInfo.transcripts[tid].RefLength - tInfo.eposInContig <= MAX_FRAGLENGTH) ||
										fivePrime
									){
										if(tInfo.eposInContig - tInfo.sposInContig < READ_LEN){
											consoleLog->error("encountered a contig shorter than read length",
															  "this is not permitted currently"
											);
											consoleLog->error("seg id: {} \t {}",tInfo.eposInContig,tInfo.sposInContig) ;
											std::exit(6) ;
										}
										localGeneProb[seg] = bfhCount ;
										localTrVector[seg].emplace_back(
											tid,
											tInfo.sposInContig,
											tInfo.eposInContig
										) ;
									}

									if(!tInfo.ore && !everRC){
										everRC = true ;
									}
								}
								if(everRC){
									localSegOreMap[seg] = false ;
								}else{
									localSegOreMap[seg] = true ;
								}
							}
						}
					}
				}
				preCalculatedSegProb[i] = localGeneProb ;
				preSegOreMapVector[i] = localSegOreMap ;
				geneSpecificTrVector[i] = localTrVector ;
			}
		}

		consoleLog->info("After loading bfh prob size {}",preCalculatedSegProb.size()) ;
		cellSegCount.resize(numCells + numOfDoublets) ;
	}


  	{

		size_t nonWhiteLisBarcodesSkipped{0} ;

		// The alevin data can be binary 
		// or can be non-binary, based on which
		// we decide it to read and fill the geneCount
		// matrix if useWhitelist is on we have to read
		// the entire matrix. 
		// Given N whitelistes and M noisy barcodes, if
		// asked minnow will generate (N+M) x numOfTranscripts
		// matrix in this order. 

		// if given in txt format.
		// the gene count matrix produced
		// by Alevin has the structure of 
		// cell to gene counts with a
		// comma as delimeter 
		// each line also ends with comma

		// we would simultaneously create the
		// cell x transcript matrix, it will have more 
		// columns presumably

		// Using a two level selection process for 
		// converting the the gene counts to transcript 
		// level counts. It is a multinomial distribution 
		// to begin with where gene counts are of type double	

		consoleLog->info("=======================Parsing the binary matrix file======================") ;
		if(!binary){
			pupulateGeneCountMatrix(countFile, geneCounts, numCells, numOfGenes) ;
		}else{
			//std::cout << "\n" ;
			consoleLog->info("After gathering all information about the matrix loading the binary Matrix") ;
			populateGeneCountMatrix(countFileBinary, originalGeneCountMatrix, numCells, numOfOriginalGenes, original2whitelistMap, original2NoisyMap, numOfOriginalCells) ;
		}

		consoleLog->info("Loaded the Matrix") ;

		// Make a truncated geneCounts file
		//T skippedCount{0} ;
		if(numOfSkippedGenes > 0){
			consoleLog->info("Truncating the matrix as not all genes are included") ;
			for(size_t cell_id = 0 ; cell_id < originalGeneCountMatrix.size(); ++cell_id){
				size_t geneCountsGeneIdx{0} ;
				for(size_t gene_id = 0 ; gene_id < numOfOriginalGenes ; ++ gene_id){
					if(skippedGenes.find(gene_id) == skippedGenes.end()){
						geneCounts[cell_id][geneCountsGeneIdx] = originalGeneCountMatrix[cell_id][gene_id] ;
						geneCountsGeneIdx++ ;
					}
				}
			}
		}else{
			geneCounts = originalGeneCountMatrix ;
		}


    	// NOTE: This feature is not tested yet
		// CREATE Doublets
		if(createDoublet){
			// Treat doublets as normal cells and put them in the list of all cells
			auto CB10XList = util::generate10XCBList(numOfDoublets,
                                               simOpts.whitelistFile,
                                               consoleLog) ;
            for(size_t i = 0; i < numOfDoublets; ++i){
                allCellNames.push_back(CB10XList[i]) ;
                allCellListMap[CB10XList[i]] = allCellNames.size() - 1 ;

                //std::cerr << "numOfWhiteListedCells "<< numOfWhiteListedCells << "\n" ;
                // generate doublet
                auto cell_id_1 = util::generateRandomNumber(0, static_cast<int>(numOfWhiteListedCells-1)) ;
                auto cell_id_2 = util::generateRandomNumber(0, static_cast<int>(numOfWhiteListedCells-1)) ;

                std::vector<T> doubletCount(numOfGenes) ;
                for(size_t j = 0; j < numOfGenes ; ++j){
                    geneCounts[sampleCells+i][j] =
            static_cast<T>((geneCounts[cell_id_1][j] + geneCounts[cell_id_2][j])/2) ;
                }
                // increase the cell numbers and update the maps
                cellNames.push_back(CB10XList[i]);
                cellDoubletMap[CB10XList[i]] = cellNames.size()-1 ;
                doubletContainer[cellDoubletMap[CB10XList[i]]] = {cell_id_1, cell_id_2} ;
            }

        }


        // By this stage the geneCount matrix should be
        // filled up and ready to be converted to the transcript level
        // counts

        //std::string diff_file_name = "diff_file.txt" ;
        //std::ofstream diffFileStream(diff_file_name) ;
        //diffFileStream << "l1" << "\t" << "l2" << "\t" << "sum_sampled" << "\t" << "sum_true" << "\t"
        //				<< "after" << "\t" << "before" << "\n" ;

        size_t cellId{0} ;
        size_t numOfExpressedGenesInput{0};
        //size_t numOfExpressedGenesOutput{0} ;

        consoleLog->info("We start to prepare Cell-Transcript matrix");
        for(auto& cellGeneCounts : geneCounts){

            // check if this cell Id is in cellWhiteListMap
            // std::cerr << "\nreading cell line\n" ;
            int droppedGeneExpression{0} ;
            uint32_t thisCellValidateGeneCount{0} ;
            //uint32_t thisCellValidateExonCount{0} ;
            uint32_t thisCellValidateTrCount{0} ;

            auto barcode = cellNames[cellId] ;
            uint32_t backupActualCellId ;

            // This step is unncessary when
            // whitelist is present, well sort of
            auto itw = cellWhiteListMap.find(barcode) ;
            auto itn = cellNoisyMap.find(barcode) ;
            auto itd = cellDoubletMap.find(barcode) ;
            if (
                (itw != cellWhiteListMap.end()) ||
                (itn != cellNoisyMap.end()) ||
                (itd != cellDoubletMap.end())
            ){

                bool trackTheCell{false} ;
                bool isNoisyCell{false} ;
                bool isDoublet{false} ;

                uint32_t actualCellId ;
                if(itn != cellNoisyMap.end()){
                    isNoisyCell = true ;
                    actualCellId = itn->second ;
                }else if(itw != cellWhiteListMap.end()){
                    actualCellId = itw->second ;
                }else{

                    actualCellId = itd->second ;
                }

                if(itw->first == cellToTrack){
                    trackTheCell = true ;
                }
                //auto actualCellId = it->second ;
                backupActualCellId = actualCellId ;


                // If clusterMap is present then
                uint32_t cluster_id{0} ;
                if(useClusters && !isNoisyCell){
                    // if it's doublet use one of the clusters
                    if(isDoublet){
                        auto parent_id = doubletContainer[actualCellId].first ;
                        cluster_id = eqClassPtr->cell2ClusterMap[parent_id] ;
                    }else{
                        cluster_id = eqClassPtr->cell2ClusterMap[actualCellId] ;
                    }

                    if(trackTheCell){ std::cerr << "the cluster id for " << barcode << " is " << cluster_id << "\n" ;}
                }

                int validateGeneCount{0} ;
                //int validateExonCount{0} ;
                //int validateTrCount{0} ;


                // treat the vector as a distribution
                std::vector<int> cellGeneCountSampled(numOfGenes,0) ;

                //double totGeneCountDouble = std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0) ;

                int totGeneCount = static_cast<int>(std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0)) ;

                //std::cerr << "[DEBUG]: Double and Int " << totGeneCountDouble << "\t" << totGeneCount << "\n" ;

                auto cellGeneCountsCopy = cellGeneCounts ;

                int toSubTract{0} ;
                for(uint32_t i = 0; i < cellGeneCounts.size() ; ++i){
                    if(cellGeneCounts[i] >= 1) {
                        --cellGeneCounts[i] ;
                        ++cellGeneCountSampled[i] ;
                        ++toSubTract ;
                    }
                }

                totGeneCount -= toSubTract ;
                // Sample form a multinomial dist with prob dist cellGeneCounts


                //std::random_device rdg;
                std::mt19937 geng(1);
                std::discrete_distribution<> dg(cellGeneCounts.begin(), cellGeneCounts.end()) ;
                for(int i = 0; i < totGeneCount ; ++i){
                    ++cellGeneCountSampled[dg(geng)] ; 
                }
                
                trueGeneCounts[cellId].assign(cellGeneCountSampled.begin(), cellGeneCountSampled.end()) ;	
                
                //diffFileStream << l1_diff << "\t" << std::fixed << std::setprecision(5) << "\t"
                //               << std::sqrt(l2_diff) << std::fixed << std::setprecision(5)  << "\t" 
                //			   << totGeneCountDouble << "\t" 
                //			   << totGeneCountRecheck <<  "\t" 
                //			   << numOfExpressedGenes3 << "\t" 
                //			   << numOfExpressedGenes2 << "\n" ;
                
                int numOfDroppedGenes{0} ;
                


                //totGeneCountDouble = std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0) ;
                //std::cerr << "DEBUG: After multinomial sampling number of expressed genes " << numOfExpressedGenes2 
                //          << "  Before number of expressed genes " << toSubTract <<  "\n" ;
                for(size_t i = 0 ; i < cellGeneCountSampled.size(); ++i){
                    // geneId
                    std::string geneName = alevinGeneIndex2NameMap[i] ;

                    int geneCount =  cellGeneCountSampled[i] ;

                    //int geneLevelTrCount{0} ;
                    //int geneLevelExCount{0} ;

                    validateGeneCount += geneCount ;
                    if (geneCount > 0){

                        bool dropThisGene{false} ;
                        numOfExpressedGenesInput++ ;

                        auto it = alevin2refMap.find(i) ;


                        if (it != alevin2refMap.end()){
                            auto originalGeneId = it->second ;
                            // fill geneCount
                            // Gene count is greater than 0 and gene count 
                            // Check if this gene is marked for polyA
                            int intronicCount{0} ;

                            if(samplePolyA){
                                auto it2 = geneIntronMap.find(originalGeneId) ;
                                if(it2 != geneIntronMap.end()){
                                    if(fractionVector[originalGeneId]){
                                         intronicCount = static_cast<int>(geneCount * fractionVector[originalGeneId]) ;
                                        geneCount = geneCount - intronicCount ;
                                    }
                                }
                            }

                            auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

                            // multinomial distribution with probabilities
                            // from the weibull distribution above
                            std::vector<double> probVec(transcriptIds.size()) ;
                            std::unordered_map<uint32_t, size_t> transcriptIdMap ;

                            for(size_t ttid = 0 ; ttid < transcriptIds.size() ; ++ttid){
                                transcriptIdMap[transcriptIds[ttid]] = ttid ;
                            }

                            // DBG based allocation
                            if(useDBG){
                                dropThisGene = false ;
                                auto segCountHist = preCalculatedSegProb[i] ;

                                if(segCountHist.size() == 0){
                                    droppedGeneExpression++ ;
                                    dropThisGene = true ;
                                    geneCount = 0 ;
                                }

                                if(!dropThisGene){
                                    std::vector<size_t> segIndex(segCountHist.size()) ;
                                    std::vector<uint32_t> segProbVec(segCountHist.size()) ;
                                    std::vector<int> segCounts(segCountHist.size(),0) ;

                                    size_t ind{0} ;
                                    for(auto it : segCountHist){
                                        segIndex[ind] = it.first ;
                                        segProbVec[ind] = it.second ;
                                        ind++ ;
                                    }

                                    std::random_device rdt4;
                                    std::mt19937 gent4(rdt4());

                                    std::unordered_map<size_t, int> segCountMap ;

                                    std::discrete_distribution<> dmt(segProbVec.begin(), segProbVec.end()) ;
                                    for(int j = 0 ; j < geneCount ; ++j){
                                        ++segCounts[dmt(gent4)] ;
                                    }
                                    for(size_t j = 0 ; j < segCounts.size() ; ++j){
                                        segCountMap[segIndex[j]] = segCounts[j] ;

                                    }

                                    cellSegCount[actualCellId][i] = segCountMap ;

                                }
                            }
                            else{
                                // In alevin mode if you are not using dbg file then by default
                                // the weibull distribution will be invoked.
                                // These are the learned distribution weibull paeameters from the paper
                                // https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005761
                                // Assign probability to individual transcripts

                                // do it only if gene count is more than 0
                                if(trackTheCell) {std::cerr << "Number of transcripts " << transcriptIds.size() << "\n";}

                                std::random_device rdt;
                                std::mt19937 gent(rdt());

                                std::weibull_distribution<> dwt(0.44,0.6);

                                std::random_device rd2;
                                std::mt19937 g2(rd2()) ;
                                std::vector<int> tIndex(transcriptIds.size()) ;
                                std::iota(tIndex.begin(), tIndex.end(), 0);
                                std::shuffle(tIndex.begin(), tIndex.end(), g2) ;

                                for(size_t j = 0; j < transcriptIds.size() ; ++j ){
                                    probVec[tIndex[j]] = dwt(gent) ;
                                }


                            }


                            //NOTE: If dbg is used we don't work with transcripts 
                            trueGeneCounts[cellId][i] = geneCount ;
                              // NOTE: In DBG mode we should be done by this
                              // point
                            if(useDBG){
                                continue ;
                            }


                            // NOTE: For weibull distribution the probability
                            // vector for each transcript could be 0 in which
                            // case we drop the gene altogether
                            auto sum_prob = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
                            if(sum_prob == 0){
                                dropThisGene = true ;
                            }

                            if (dropThisGene){
                                //cellExCount[actualCellId][i] = {} ;
                                numOfDroppedGenes++ ;
                                continue ;
                            }

                              // Otherwise we distribute the gene counts to
                              // individual transcripts
                            std::random_device rdt3;
                            std::mt19937 gent3(rdt3());
                            std::discrete_distribution<> dmt(probVec.begin(), probVec.end()) ;
                            // saple from this distribution geneCount times
                            std::vector<int> transcriptCounts(transcriptIds.size(), 0) ;
                            for(int j=0 ; j < geneCount ; ++j){
                                ++transcriptCounts[dmt(gent3)] ;
                            }

                            //auto beforeFuckingUp = std::accumulate(
                            //	transcriptCounts.begin(),
                            //	transcriptCounts.end(),
                            //	0
                            //);


                            // total gene count
                            T totCount{0} ;
                            uint32_t lastAlevinTid{0} ;
                            for(size_t j = 0; j < transcriptIds.size(); ++j){
                                auto tid = transcriptIds[j] ;
                                auto alevinTid = alevinReverseMap[tid] ;
                                data[actualCellId][alevinTid] = transcriptCounts[j] ;
                                totCount += data[actualCellId][alevinTid] ;
                                lastAlevinTid = alevinTid ;
                            }

                            thisCellValidateTrCount += totCount ;
                            thisCellValidateGeneCount += geneCount ;

                            // Add intronic count if there is a
                            // non-zero quantity to add
                            if (intronicCount > 0){
                                auto it3 = intronCountMap.find(actualCellId) ;
                                if (it3 != intronCountMap.end()){
                                    intronCountMap[actualCellId].push_back(std::make_pair(lastAlevinTid, intronicCount)) ;
                                }else{
                                    intronCountMap[actualCellId] = {{lastAlevinTid, intronicCount}} ;
                                }
                            }

                        }else{
                            consoleLog->error("This gene id: {} not found\n", i) ;
                        }
                    }

                }

                auto iVec = intronCountMap[backupActualCellId] ;
                int intronicSum = 0 ;
                for(auto i : iVec){
                    intronicSum += i.second ;
                }
                //std::cerr << barcode << "\t" << backupActualCellId << "\t" << cellId <<  "\t" <<  validateGeneCount << "\t" << intronicSum  << "\n" ;

            }else{
                //std::cerr << "barcode not found " << barcode << "\n" ;
                nonWhiteLisBarcodesSkipped++ ;
            }

            //if(useDBG){
      //  consoleLog->info("NOTE: no gene should be skipped die to continue block, this should not be printed");
            //	consoleLog->info("DBG::: total Gene Expression skipped {}", droppedGeneExpression) ;
            //}


            cellId += 1 ;
            _verbose("\rNumber of cells processed : %lu", cellId);

        }

        //std::cerr << "DEBUG: ==== Before dumping input expressedGeneCount " << numOfExpressedGenesInput << "\n" ;
        std::cout << "\n" ;
        consoleLog->info("{} cell barcodes are skipped as they are not in whitelist", nonWhiteLisBarcodesSkipped) ;
    }

    // FIXME: this should not be ad-hoc
    numCells = geneCounts.size() ;

    if(!useDBG){
        consoleLog->info("The transcript matrix is constructed, with dimension {} x {} \n" 
                        "\t\t\t\t\tfrom gene count with dimention {} x {}, geneCounts.size(): {}",
                        data.size(), data[0].size(), cellNames.size(), alevin2refMap.size(), geneCounts.size()) ;
    }

    std::string cellColFile = outDir + "/alevin/quants_mat_cols.txt" ;
    std::ofstream cellColStream(cellColFile.c_str()) ;
    
    for(size_t id = 0 ; id < alevinGeneIndex2NameMap.size() ; ++id){
        cellColStream << alevinGeneIndex2NameMap[id] << "\n" ; 
    }

}


// NOTE: Main module for reading splatter
template <typename T>
void readSplatterMatrix(
    std::string& countFile,
    std::vector<std::vector<T>>& geneCount,
    size_t numCells,
    uint32_t sampleCells,
    size_t numOfGenes
){
    if(! util::fs::FileExists(countFile.c_str())){
        std::cerr << "quants_mat.csv file does not exist\n" ;
        std::exit(1) ;
    }
    if (sampleCells > 0){
        numCells = sampleCells;
    } 

    std::ifstream dataStream(countFile.c_str()) ;
    std::string line ;
    std::cout << "Debug::matrix size " << geneCount.size() 
                << " x " << geneCount[0].size() << "\n";
    std::cout << "cell count " << numCells << "\n";

    size_t geneId{0} ;
    while(std::getline(dataStream, line)){

        std::vector<std::string> cellCountsString ; 
        util::split(line, cellCountsString, ",") ;

        for(size_t cellId = 0 ; cellId < numCells ; ++cellId){
            geneCount[cellId][geneId] = std::stoi(cellCountsString[cellId]) ; 
        }
        geneId += 1 ;
        _verbose("\rIn Splatter: Number of genes processed : %lu", geneId);
        if(geneId == numOfGenes){
            break ;
        }
    }

}


void readUniqueness(
     std::vector<std::pair<std::string, double>>& uniquenessInfo,
       std::unordered_map<std::string, uint32_t>& geneMap,
       std::string& uniquenessListFile,
     bool useReverse
){

    //std::string uniquenessListFile{"../data/hg/hg_stranded_gene_uniqueness.txt"} ;
    if(!util::fs::FileExists(uniquenessListFile.c_str())){
        std::cerr << "[ERROR] The gene uniqueness file:: " << uniquenessListFile << " does not exist, should exist, EXITING !!"; 
        std::exit(1) ;
    }

    std::ifstream uniqueFileStream(uniquenessListFile) ;

    std::string line ;
    // throw the fist line 
    std::getline(uniqueFileStream, line) ;

    while(std::getline(uniqueFileStream, line)){
        // strip new line
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        std::vector<std::string> tokens ;
        util::split(line, tokens, "\t") ;
        auto geneName = tokens[0] ;
        double score{0.0} ;
        if(std::stod(tokens[1]) != 0){
            score = std::stod(tokens[1]) / std::stod(tokens[2]) ;
        }
        if(geneMap.find(geneName) != geneMap.end()){
            uniquenessInfo.push_back({tokens[0], score}) ;
        }
    }

    using ptype = std::pair<std::string, double> ;

  // This should keep uniq genes at top 
  if(useReverse){
        std::sort(uniquenessInfo.begin(), uniquenessInfo.end(), [](ptype &left, ptype &right) {
            return left.second > right.second;
        });
    }else{
    // This should keep least uniq genes at top
        std::sort(uniquenessInfo.begin(), uniquenessInfo.end(), [](ptype &left, ptype &right) {
            return left.second < right.second;
        });
    }
}

void makeUniquenessBins(
    std::vector<std::unordered_set<std::string>>& uniqBins,
    std::unordered_map<std::string, uint32_t>& geneMap,
    std::string& uniquenessListFile
){


    //std::string uniquenessListFile{"../data/hg_stranded_gene_uniqueness.txt"} ;

    if(!util::fs::FileExists(uniquenessListFile.c_str())){
        std::cerr << uniquenessListFile << " does not exist, should exist, EXITING !!"; 
        std::exit(1) ;
    }

    std::ifstream uniqueFileStream(uniquenessListFile) ;
    std::vector<std::pair<std::string, double>> uniquenessInfo ;	
    uniqBins.resize(10) ;
    std::vector<double> scoreBins(11) ;

    int ind{0} ;
    for(double i = 0 ; i <= 1.0 ; i += 0.1){
        scoreBins[ind] = i ;
        ind++ ;
        //std::cerr << "score ::: " << ind << "\t" << scoreBins[ind-1] << "\n" ;
    }

    std::string line ;
    // throw the fist line 
    std::getline(uniqueFileStream, line) ;


    while(std::getline(uniqueFileStream, line)){
        // strip new line
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        std::vector<std::string> tokens ;
        util::split(line, tokens, "\t") ;
        auto geneName = tokens[0] ;
        double score{0.0} ;
        if(std::stod(tokens[1]) != 0){
            score = std::stod(tokens[1]) / std::stod(tokens[2]) ;
        }
        bool gotin{false} ;
        if(geneMap.find(geneName) != geneMap.end()){
            uniquenessInfo.push_back({tokens[0], score}) ;
            for(int i = 0; i < 10 ; ++i){
                //if(lookin) std::cout << scoreBins[i] << "\t" << scoreBins[i+1] << "\t" << score << "\n";
                if((score >= scoreBins[i]) && (score <= scoreBins[i+1])){
                    gotin = true ;
                    uniqBins[i].insert(tokens[0]) ;
                    break ;
                }
                if(score == 1){
                    uniqBins[9].insert(tokens[0]) ;	
                    gotin = true ;
                }
            }

            if(!gotin){
                std::cout <<"Score::: " << score << "\n" ;
                std::exit(1) ;
            }
        
        }
    }

    //for(size_t i = 0 ; i < uniqBins.size() ; ++i){
    //	std::cerr << scoreBins[i] << "\t" << uniqBins[i].size() << "\n" ;
    //}


}


template<typename T>
void DataMatrix<T>::loadSplatterData(
    SimulateOptions& simOpts,
    Reference& refInfo
){

	// Load data from the simOpts
	auto splatterDir = simOpts.inputdir ;
	size_t sampleCells = simOpts.sampleCells ;
	std::string outDir = simOpts.outDir ;
	std::string bfhFile = simOpts.bfhFile ;
	std::string gfaFile = simOpts.gfaFile ;
	std::string uniquenessListFile = simOpts.uniquenessFile ;
	bool fivePrime = true;
	//bool useEqClass = simOpts.useEqClass ;
	bool useWeibull = simOpts.useWeibull ;
	bool useDBG = simOpts.useDBG ;
	bool velocityMode = simOpts.velocityMode;


	std::string geneListFile = splatterDir + "/quants_mat_rows.txt" ;
	std::string cellListFile = splatterDir + "/quants_mat_cols.txt" ;

	// this is gene to cell matrix
	std::string countFile = splatterDir + "/quants_mat.csv" ;

	if(! util::fs::FileExists(geneListFile.c_str())){
		consoleLog->error("{} file does not exist", geneListFile) ;
		std::exit(1) ;
	}
	if(! util::fs::FileExists(cellListFile.c_str())){
		consoleLog->error("{} file  not exist",cellListFile) ;
		std::exit(1) ;
	}

	if(! util::fs::FileExists(countFile.c_str())){
		consoleLog->error("{} file does not exist",countFile) ;
		std::exit(1) ;
	}

  // Check command combinations first, since all combinations
  // are not allowed in splatter


    // Splatter does not specify gene names
    // We would randomly sample genes from from
    // the reference provided. Everything will be treated
    // as whitelisted cells.
    // update txp to gene


    refInfo.updateGene2TxpMap() ;

    auto& geneMap = refInfo.geneMap ;

  //if (geneMap.find("ENSG00000207552.1") != geneMap.end()){
  //  std::cout << "We have it !!!!!" ;
  //  std::exit(1) ;
  //}else{
  //  std::cout << "FUCK\n" ;
  //  std::exit(1) ;
  //}


    consoleLog->info("Number of genes in the txp2gene file: {}", geneMap.size()) ;
    
    consoleLog->info(
        "Parsing {}/quants_mat_cols.txt",splatterDir) ;	
    
    std::map<std::string, uint32_t> allCellListMap ;
    std::vector<std::string> allCellNames ;
    std::string cellToTrack  = ""; 
    //uint32_t cellIdToTrack{0} ;
    {
        std::ifstream cellNameStream(cellListFile) ;
        std::string line ;
        while(std::getline(cellNameStream, line)){
            // strip new line
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            allCellNames.push_back(line) ;
            allCellListMap[line] = allCellNames.size() - 1 ;
        } 
    }
    // In splatter there is no whitelist but there can be clusters 
  std::cout<<"=======================Reading Splatter Matrix=====================\n" ;
    consoleLog->info("{} cells are present ", allCellNames.size()) ;
    consoleLog->info("Start parsing Splatter output") ;
    consoleLog->info(
        "Parsing {}/quants_mat_rows.txt",splatterDir) ;
    
    std::unordered_map<std::string, std::string> splatterGeneMap ; //splatter gene name -> real gene name
    std::vector<std::string> splatterGeneNames ;
    //std::unordered_map<uint32_t, std::string> aleviGeneIndex2NameMap ;

    {
        std::ifstream geneNameStream(geneListFile) ;
        std::string line ;
        while(std::getline(geneNameStream, line)){
            // strip new line
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            splatterGeneNames.push_back(line) ;
        }
    }


    if(simOpts.numMolFile != ""){
        if(! util::fs::FileExists(simOpts.numMolFile.c_str())){
            consoleLog->warn("--numMolFile {} does not exist",simOpts.numMolFile) ;
        }else{
            std::ifstream dupCountStream(simOpts.numMolFile) ;
            std::string line ;
      size_t line_counter{0};
            while(std::getline(dupCountStream, line)){
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                std::vector<std::string> tokens ;
                util::split(line, tokens, "\t") ;
                if (tokens.size() == 2){
                    cellNamesDupCount[tokens[0]] = std::stoul(tokens[1]) ;
                }else{
          consoleLog->warn("{} of {} has more than 2 entries skipping", line_counter, simOpts.numMolFile) ;
          cellNamesDupCount.clear() ;
          break ;
        }
        line_counter += 1 ;
            }

        }
    }

    size_t numOfGenes = splatterGeneNames.size() ;
    size_t numOfOriginalGenes = numOfGenes ;


    // Let's assume no whitelist
    {
        cellNames = allCellNames ;
        cellWhiteListMap = allCellListMap ;
    }

    if ((sampleCells > 0) and (sampleCells < cellNames.size())){
        consoleLog->info("According to spefication We will sample {} from {} cells",sampleCells,cellNames.size());
    }else{
        sampleCells = cellNames.size() ;
    }

    numCells = sampleCells ;

    std::vector<std::vector<T>> originalGeneCountMatrix ;


    originalGeneCountMatrix.resize(sampleCells, std::vector<T>(numOfGenes)) ;
    consoleLog->info("Debug:: resizing {}",originalGeneCountMatrix.size());
    //geneCounts.resize(sampleCells, std::vector<T>(numOfGenes)) ;

    //trueGeneCounts.resize(sampleCells, std::vector<int>(numOfGenes)) ;

    std::unordered_set<uint32_t> skippedGenes ;

    for(auto& v : originalGeneCountMatrix)
        std::memset(&v[0], 0, sizeof(v[0]) * v.size());


    consoleLog->info("Debug:: reading splatter");
    readSplatterMatrix(
        countFile,
        originalGeneCountMatrix,
        allCellNames.size(),
        numCells,
        splatterGeneNames.size()
    ) ;



  std::cout<<"==================Done Parsing Splatter Matrix==================\n" ;

  std::string gene_name_to_track = ""; // "ENSG00000001084.13";
  size_t geneIdToTrack = std::numeric_limits<size_t>::max() ;

    consoleLog->info(
                        "Splatter matrix is read, with dimension {} x {}", 
                        originalGeneCountMatrix.size(), 
                        originalGeneCountMatrix[0].size()
                    ) ;
    std::unordered_map<uint32_t, std::string> alevinGeneIndex2NameMap ;
    size_t numOfSkippedGenes{0} ;
  std::set<size_t> validGeneIds;
  std::set<std::string> invalidGeneNames;

    {
        // go over gene to transcript map, include a
        // gene if contains a transcript that is present
        // in the reference fasta file

        bool testUniqueness{simOpts.testUniqness} ;
    bool revUniqness{simOpts.reverseUniqness} ;

        if(testUniqueness || revUniqness){
            // read uniqueness file
      std::cout << "\n !!!!!!!!!!!!!!!!!! IN TEST UNIQUENESS MODE !!!!!!!!!!!!!!!!!!!!!!!\n" ;

            std::vector<std::pair<std::string, double>> uniquenessInfo ;
            std::vector<std::pair<size_t, int>> countInfo ;

            readUniqueness(uniquenessInfo, refInfo.geneMap, uniquenessListFile, revUniqness) ;
            size_t mostUniqGeneId{0} ;
            size_t leastUniqGeneId{uniquenessInfo.size() - 1} ;

            for(size_t geneId = 0 ; geneId < splatterGeneNames.size() ; geneId++){
                int gcount{0} ;
                for(size_t cellId = 0 ; cellId < allCellNames.size() ; cellId++){
                    if(originalGeneCountMatrix[cellId][geneId] > 0){
                        gcount++ ;
                    }
                }

        // std::cout << "\n gcount " << gcount << "\n" ;
                countInfo.push_back({geneId, gcount}) ;
                _verbose("\rIn Splatter: Number of gene sum processed : %lu", geneId);
            }
            using ptype = std::pair<size_t, int> ;

            std::sort(countInfo.begin(), countInfo.end(), [](ptype &left, ptype &right) {
                return left.second > right.second;
            });

            std::cerr << "\nIn Splattermode countInfo.size() " << countInfo.size() 
                      << "\t uniquenessInfo.size() " << uniquenessInfo.size() 
                      << "\n" ;

            for(auto it : countInfo){
                //search until you find a proper gene
                while( geneMap.find(uniquenessInfo[mostUniqGeneId].first) == geneMap.end()){
                    mostUniqGeneId++ ;
                }

                //std::cerr << "found " << uniquenessInfo[mostUniqGeneId].first << "\n" ;
                if (mostUniqGeneId <= leastUniqGeneId){
                    auto splatterGeneId = it.first ;
                    auto spgName = splatterGeneNames[splatterGeneId] ;

                    auto realGeneName = uniquenessInfo[mostUniqGeneId].first ;
                    auto refGeneId = geneMap[realGeneName] ;

                    splatterGeneMap[realGeneName] = spgName ; // keep the real name

                    alevin2refMap[splatterGeneId] = refGeneId ; // gene id to real gene id
                    // this one is extra
                    alevinGeneIndex2NameMap[splatterGeneId] = realGeneName ; // keep mapping from real gene name

                }
                mostUniqGeneId++ ;
            }
        }else if(simOpts.normalMode){

      std::cout << "\n !!!!!!!!!!!!!!!!!! IN NORMAL MODE !!!!!!!!!!!!!!!!!!!!!!!\n" ;

            size_t splatterGeneId{0} ;
            while(splatterGeneId < splatterGeneNames.size()){


                auto spgName = splatterGeneNames[splatterGeneId] ;

                if(geneMap.find(spgName) != geneMap.end()){
                    auto selectedGeneName = spgName ;
                    auto selectedGeneId = geneMap[spgName] ;
                    splatterGeneMap[selectedGeneName] = spgName ;
                    alevin2refMap[splatterGeneId] = selectedGeneId ;
                    alevinGeneIndex2NameMap[splatterGeneId] = selectedGeneName ;
                }else{
                    //consoleLog->error("Possible FATAL ERROR: The reference file contains genes which are not present in the reference transcript to gene map") ;
                    //consoleLog->error("If you think that is not the error post an issue in n https://github.com/COMBINE-lab/minnow/issues with transcript reference") ;
                    skippedGenes.insert(splatterGeneId) ;
                    numOfSkippedGenes++ ;
                }


                splatterGeneId++ ;
    
                _verbose("\rIn Splatter: Number of genes processed : %lu", splatterGeneId);

            }

        }
        else if(useDBG){

      std::cout << "\n !!!!!!!!!!!!!!!!!! IN DBG MODE !!!!!!!!!!!!!!!!!!!!!!!\n" ;

      if(!util::fs::FileExists(gfaFile.c_str())){
        consoleLog->error("--dbg is invoked but is not passed with --gfa, exiting");
        std::exit(2) ;
      }

      // Load transcripts that are present in the gfa file 
      dbgPtr = new GFAReader(gfaFile, consoleLog) ;
      dbgPtr->parseFile(refInfo, outDir) ;


      // Take only genes that are present in trSegmentMap ;
      // std::set<size_t> validGeneIds ;
      for(auto& segObj : dbgPtr->trSegmentMap){
        auto& validTid = segObj.first ;
        auto geneId = refInfo.transcript2geneMap[validTid] ;
        validGeneIds.insert(geneId) ;
        // get corresponding gene id
      }
      consoleLog->info("The size of the gene id pool {}", validGeneIds.size()) ;

      // construct a reverse gene map
      std::unordered_map<uint32_t, std::string> reverseGeneMap ;
      for(auto geneMapIt: geneMap){
        reverseGeneMap[geneMapIt.second] = geneMapIt.first ;
      }

      size_t splatterGeneId{0} ;
      size_t truncatedGeneId{0} ;
      bool customNames = simOpts.customNames ;
      if(customNames){
        // gene names are already provided
        // fill up invalid gene ids
        while(splatterGeneId < splatterGeneNames.size()){
          auto spgName = splatterGeneNames[splatterGeneId] ;
          if(geneMap.find(spgName) != geneMap.end()){
            auto geneId = geneMap[spgName] ;
            if(validGeneIds.find(geneId) != validGeneIds.end()){
              splatterGeneMap[spgName] = spgName ;
              alevin2refMap[truncatedGeneId] = geneId ;
              alevinGeneIndex2NameMap[truncatedGeneId] = spgName ;
            //   if (spgName == gene_name_to_track){
            // 	  std::cout << "===========================================================gene id missed " << geneId << "\n"; 
            // 	  geneIdToTrack = truncatedGeneId ;
            //   }
              truncatedGeneId++ ;
            }else{
              //consoleLog->warn("the gene name {} is not present in de-Bruijn graph",spgName) ;
              skippedGenes.insert(splatterGeneId) ;
              invalidGeneNames.insert(spgName + "\t" + "Absent-De-Bruijn");
              numOfSkippedGenes++ ;
            }
          }else{
            //consoleLog->warn("the gene name {} is not present in reference",spgName) ;
            skippedGenes.insert(splatterGeneId) ;
            invalidGeneNames.insert(spgName + "\t" + "Absent-Gene-Txp");
            numOfSkippedGenes++ ;
          }

          splatterGeneId++ ;
          _verbose("\rIn Splatter: Number of genes processed : %lu", splatterGeneId);

        }
      }else{
        auto geneIdIt = validGeneIds.begin() ;
        while(splatterGeneId < splatterGeneNames.size()){
          if(geneIdIt != validGeneIds.end()){
            auto selectedGeneId = *geneIdIt ;
            auto selectedGeneName = reverseGeneMap[selectedGeneId] ;
            auto spgName = splatterGeneNames[splatterGeneId] ;
            splatterGeneMap[selectedGeneName] = spgName ;
            alevin2refMap[truncatedGeneId] = selectedGeneId ;
            alevinGeneIndex2NameMap[truncatedGeneId] = selectedGeneName ;
            truncatedGeneId++ ;
            ++geneIdIt ;
          }else{
            skippedGenes.insert(splatterGeneId) ;
            numOfSkippedGenes++ ;
          }

          splatterGeneId++ ;

          _verbose("\rIn Splatter: Number of genes processed : %lu", splatterGeneId);
        }
      }


        }else{

      size_t splatterGeneId{0} ;
            // Assign a random gene to the columns of the matrix  
            // Go over geneMap assign a gene that exists
            // and has at least one transcript that has more than 
            // READ_LEN 
            std::vector<std::unordered_set<std::string>> uniqBins ;
            std::vector<int> uniqHist ;

            makeUniquenessBins(uniqBins, geneMap, uniquenessListFile) ;
            uniqHist.resize(uniqBins.size()) ;

            //std::cerr << "Made uniqueness bins " << uniqBins.size() << "\t" << uniqHist.size() << "\n"; 
            size_t totalGenes{0} ;

            for(size_t i = 0; i < uniqBins.size(); ++i){
                uniqHist[i] = uniqBins[i].size() ;
                std::cout << "Bin " << i << " " << uniqHist[i] << " " << uniqBins[i].size() << "\n" ;
                totalGenes +=uniqHist[i] ; 
            }

            //std::cout << "Total genes " << totalGenes << "\t" << geneMap.size() << "\n";

            while(splatterGeneId < splatterGeneNames.size()){
                std::random_device rd;
                std::mt19937 gen(rd());	

                std::discrete_distribution<> d(uniqHist.begin(), uniqHist.end());
                int ind = d(gen) ;
                //while(!uniqHist[ind]){
                //	for(size_t i = 0; i < uniqBins.size(); ++i){
                //		uniqHist[i] = uniqBins[i].size() ;
                //		//std::cout << "Bin " << i << " " << uniqHist[i] << "\n" ;
                //	}	
                //	ind = d(gen) ;
                //	//std::exit(1) ;
                //}

                auto randIt = uniqBins[ind].begin() ;
                auto selectedGeneName = *randIt ;
                auto selectedGeneId = geneMap[selectedGeneName] ;

                uniqBins[ind].erase(randIt) ;	
                uniqHist[ind] = uniqHist[ind] - 1 ;

                auto spgName = splatterGeneNames[splatterGeneId] ;
                splatterGeneMap[selectedGeneName] = spgName ;
                alevin2refMap[splatterGeneId] = selectedGeneId ;
                alevinGeneIndex2NameMap[splatterGeneId] = selectedGeneName ;

                splatterGeneId++ ;	

                _verbose("\rIn Splatter: Number of genes processed : %lu", splatterGeneId);
            }

    }
        if(splatterGeneMap.size() > splatterGeneNames.size()){
            consoleLog->error("The splatter matrix contains more genes than provided in reference, truncating") ;
            consoleLog->error("splatterGeneMap.size() {} splatterGeneNames.size() {} geneMap.size() {}", 
                                splatterGeneMap.size(), 
                                splatterGeneNames.size(),
                                geneMap.size()
                            ) ;
            std::exit(2) ;
        }
    }

	// Make a truncated geneCounts file
	//T skippedCount{0} ;
	geneCounts.resize(sampleCells, std::vector<T>(numOfGenes - numOfSkippedGenes)) ;
	trueGeneCounts.resize(sampleCells, std::vector<int>(numOfGenes - numOfSkippedGenes)) ;
	for(auto& v : geneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : trueGeneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());


	if(numOfSkippedGenes > 0){
		if(!useDBG){
		consoleLog->warn("Skipping {} genes, either they are short or absent in reference",numOfSkippedGenes) ;
		}else{
		consoleLog->warn("Skipping {} genes, gene pool size of de-Bruijn graph {}",
						numOfSkippedGenes, validGeneIds.size()) ;
		}

		for(size_t cell_id = 0 ; cell_id < originalGeneCountMatrix.size(); ++cell_id){
			size_t geneCountsGeneIdx{0} ;
			for(size_t gene_id = 0 ; gene_id < numOfOriginalGenes ; ++ gene_id){
				if(skippedGenes.find(gene_id) == skippedGenes.end()){
					geneCounts[cell_id][geneCountsGeneIdx] = originalGeneCountMatrix[cell_id][gene_id] ;
					geneCountsGeneIdx++ ;
				}
			}
		}
		// write down the skipped genes if there is any
		if (invalidGeneNames.size() > 0){
			std::string invalidGeneFile = outDir + "/invalidnames.txt" ;
			std::ofstream invalidGeneNameStream(invalidGeneFile);
			for(auto& n : invalidGeneNames){
				invalidGeneNameStream << n << "\n" ;
			}
		}


		if (geneCounts[0].size() == 0){
			consoleLog->error("Skipped all genes, please check if gene names are consistent"
			                  " across transcript to gene mapping file and others"
			);
			std::exit(13);
		}

		consoleLog->info("Truncated the matrix to dimension {} x {}",geneCounts.size(), geneCounts[0].size()) ;
	}else{
		geneCounts = originalGeneCountMatrix ;
	}


	numOfGenes = geneCounts[geneCounts.size()-1].size() ;

	auto& gene2transcriptMap = refInfo.gene2transcriptMap ;
	std::map<uint32_t, uint32_t> alevinReverseMap ;
	{
		// create a transcript map
		uint32_t trId{0} ;
		for(auto it : alevin2refMap){
			auto gIt = gene2transcriptMap.find(it.second) ;
			if(gIt != gene2transcriptMap.end()){
				auto& trVec = gIt->second ;
				//std::cerr << "# of transcripts in this gene "<< gIt->first 
				//	      <<  ":\t" << trVec.size() << "\n" ;	
				for(auto tr : trVec){
					alevin2refTranscriptMap[trId] = tr ;
					alevinReverseMap[tr] = trId ;
					trId += 1 ;
				}
			}else{
				consoleLog->error("This should not happen: gene {} not found",gIt->first) ;
			}
		}
	}


	numOfTranscripts = alevin2refTranscriptMap.size() ;


	// If we use dbg then the transcripts used above won't be needed  
	if(useDBG){

		auto rspdFile = simOpts.rspdFile ;
		if(rspdFile != ""){
			if(util::fs::FileExists(rspdFile.c_str())){
				std::cerr << "RSPD file " << rspdFile << "\n" ;
				rspdVec.resize(MAX_FRAGLENGTH) ;
				std::ifstream rspdStream(rspdFile.c_str()) ;
				std::string line ;
				while(std::getline(rspdStream, line)){
					std::vector<std::string> tokens ;
					util::split(line, tokens, ",") ;
					if(tokens.size() != 2)
						continue ;
					uint32_t pos = std::stoul(tokens[0]) ;
					uint32_t prob = std::stoul(tokens[1]) ;
					if(pos < rspdVec.size()){
						rspdVec[pos] = prob ;
					}
				}
			}else{
				consoleLog->warn("RSPD file is not empty and doesn't exist, going with truncated sampling\n") ;
			}
		}else{
			consoleLog->warn("No RSPD file provided \n") ;
		}


		
      // stand alone call to bfh
      // ignore other eqclass stuff
      // FIXME we don't call equivalence
      // class folder any more
      //std::string eqFileDir = "/mnt/scratch1/hirak/minnow/metadata/hg/" ;
      eqClassPtr = new BFHClass(consoleLog) ;
      if(bfhFile != ""){
        eqClassPtr->loadBFH(
                            bfhFile, 
                            simOpts.clusterFile, 
                            refInfo, 
                            cellWhiteListMap, 
                            false, 
                            cellNoisyMap,
                            simOpts.outDir
                            ) ;
        consoleLog->info("Loaded the bfh.txt file") ;
      }else if(useDBG){
        //load default file
        std::string countProbFile{simOpts.countProbFile} ;
        if(countProbFile == ""){
          consoleLog->warn("picking hard coded count prob file") ;
          countProbFile = "../data/hg/countProb_pbmc_4k.txt" ;
        }
        if(util::fs::FileExists(countProbFile.c_str())){
            eqClassPtr->loadProbability(
                                        countProbFile,
                                        refInfo,
                                        false
                                        ) ;
          }else{
            consoleLog->error("Invoked with DBG mode but --countProb file is not present") ;
            consoleLog->error("Add --countProb option with the file countProb_pbmc_4k.txt") ;
            std::exit(10) ;
          }
      }


      consoleLog->info("The size of probability Vector {}",eqClassPtr->countProbability.size()) ;
	  // the size should not be zero
	  if (eqClassPtr->countProbability.size() == 0){
		  consoleLog->error("Recheck the bfh file, the probability vector size is zero");
		  std::exit(11);
	  }

		preCalculatedSegProb.resize(numOfGenes) ; // per gene
		preSegOreMapVector.resize(numOfGenes) ;
		geneSpecificTrVector.resize(numOfGenes) ;

		// fill segment based probability
		// gene to tr to prob mapping

		for(uint32_t i = 0; i < numOfGenes; ++i){
			auto it = alevin2refMap.find(i) ;
			if(it != alevin2refMap.end()){
				auto originalGeneId = it->second ;
				// if (i == geneIdToTrack)
				// 	std::cerr << "Tracking gene id " << i << " original gene id " << originalGeneId << "\n" ; 
				//std::cerr << "original gene id "  << originalGeneId << "\n" ;

				auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

				std::unordered_map<size_t, uint32_t> localGeneProb ;
				std::unordered_map<size_t, bool> localSegOreMap ;
				std::unordered_map<size_t, std::vector<trInfo>> localTrVector ;

				//std::cerr << "tr size " << transcriptIds.size() << "\n" ;
				size_t shortTranscriptLength{0};
				bool shortLength{false} ;
				uint32_t shortTid{0};
				bool tqVecZero{false} ;
				size_t numTids = transcriptIds.size() ;

				
				for(auto tid : transcriptIds){
					refInfo.transcripts[tid].setGeneId(i);
					if(dbgPtr->trSegmentMap.find(tid) != dbgPtr->trSegmentMap.end()){
						auto segVec = dbgPtr->trSegmentMap[tid] ;
						bool everRC{false} ;
						// check if this transcript has enough transcripts
						bool eligibleSegments{false};
						for(auto seg: segVec){
							auto segCount = dbgPtr->eqClassMapCounts[seg] ;
							//auto bfhCount = eqClassPtr->getGeneLevelProbCount(
							//	originalGeneId,
							//	segCount 
							//) ;
							auto bfhCount = eqClassPtr->countProbability[segCount] ;
							if(bfhCount != 0){
								eligibleSegments = true ;
								break ;
							}
						}
						
						double uniformProb{0.0};
						if(!eligibleSegments && (segVec.size() > 0)){
							uniformProb = static_cast<double>(1.0) / static_cast<double>(segVec.size());
						}

						for(auto seg : segVec){
							auto segCount = dbgPtr->eqClassMapCounts[seg] ;
							//auto bfhCount = eqClassPtr->getGeneLevelProbCount(
							//	originalGeneId,
							//	segCount 
							//) ;
							auto bfhCount = eqClassPtr->countProbability[segCount] ;
							// if this multimapping is not seen in count probability
							// vector then just 
							if (bfhCount == 0 && !eligibleSegments){
								bfhCount = uniformProb ;
							}else if (bfhCount == 0){
								continue;
							}

							//if(bfhCount != 0){
								// old one
								// localTrVector[seg].emplace_back(tid) ;

							auto tcInfoVec = dbgPtr->eqClassMap[seg][tid] ;
							tqVecZero = (tcInfoVec.size() == 0) ;
							for(auto tInfo : tcInfoVec){
								if (!fivePrime){
									if(
										(refInfo.transcripts[tid].RefLength - tInfo.eposInContig <= MAX_FRAGLENGTH)
									){
										localGeneProb[seg] = bfhCount ;
										localTrVector[seg].emplace_back(
											tid,
											tInfo.sposInContig,
											tInfo.eposInContig
										) ;
									}else{
										shortTid = tid ;
										shortLength = true ;
										shortTranscriptLength = refInfo.transcripts[tid].RefLength - tInfo.eposInContig;
									}
								}else{
									if(
										tInfo.sposInContig <= MAX_FRAGLENGTH - READ_LEN
									){
										localGeneProb[seg] = bfhCount ;
										localTrVector[seg].emplace_back(
											tid,
											tInfo.sposInContig,
											tInfo.eposInContig
										) ;
									}
								}

								if(!tInfo.ore && !everRC){
									everRC = true ;
								}
							}
							if(everRC){
								localSegOreMap[seg] = false ;
							}else{
								localSegOreMap[seg] = true ;
							}
							//}
						}
					}else{
						if(numTids == 1){
						consoleLog->warn("There is one transcript {} for this gene, length: {}, name: {}",
											tid, refInfo.transcripts[tid].RefLength,
											refInfo.transcripts[tid].RefName);
						}
          			}
                }

        if(localGeneProb.size() == 0){
          consoleLog->error("The gene got skipped, it should not") ;
          consoleLog->error("shortLength {}, length {}, tqVecZero {} numTids {} transcript name {} transcript id {}",
                            shortLength,shortTranscriptLength,tqVecZero,numTids, refInfo.transcripts[shortTid].RefName,
                            shortTid) ;
          std::exit(3) ;
        }
                preCalculatedSegProb[i] = localGeneProb ;
                preSegOreMapVector[i] = localSegOreMap ;
                geneSpecificTrVector[i] = localTrVector ;
                // if(i == geneIdToTrack){
                // 	std::cerr << " Corresponding probability sizes for gene id " << geneIdToTrack 
                //               << " " << localGeneProb.size() << "\t"
                // 			  << localSegOreMap.size() << "\t"
                // 			  << localTrVector.size() << "\n" ;
                // }

            }
        }

        std::cerr << "SPLATTER MODE: After loading bfh prob size " << preCalculatedSegProb.size() << "\n" ;
        cellSegCount.resize(numCells + numOfDoublets) ;
    }else{

        data.resize(sampleCells, std::vector<T>(numOfTranscripts)) ;
        for(auto& v : data)
                std::memset(&v[0], 0, sizeof(v[0]) * v.size());
    }


    size_t numOfDroppedGenes{0} ;
    size_t numOfBadCells{0} ;
    // create transcript level matrix 
    {
        size_t cellId{0} ;
        for(auto&  cellGeneCounts : geneCounts){
            size_t droppedGeneExpression{0} ;
            size_t cellSum{0} ;

              size_t cellSumCheck{0} ;
            auto barcode = cellNames[cellId] ;
            for(size_t i = 0 ; i < cellGeneCounts.size() ; ++i){
                std::string geneName = alevinGeneIndex2NameMap[i] ;

                int geneCount = cellGeneCounts[i] ;
                // if(geneName == gene_name_to_track && barcode == "Cell1"){
                // 	std::cerr << "gene id " << i << "\t"
                // 			  << "barcode " << barcode << "\t"
                //               << "gene count " << geneCount << "\n";
                // }

                // std::cerr << "\n new gcount " << geneCount << "\n" ;
                //int geneLevelTrCount{0} ;
                //int geneLevelExCount{0} ;

                if(geneCount > 0){

                    cellSum += geneCount ;

                    bool dropThisGene{false} ;

                    auto originalGeneId = alevin2refMap[i] ;
                    auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

                    std::vector<double> probVec(transcriptIds.size()) ;
                    std::unordered_map<uint32_t, size_t> transcriptIdMap ;

                    for(size_t ttid = 0 ; ttid < transcriptIds.size() ; ++ttid){
                        transcriptIdMap[transcriptIds[ttid]] = ttid ; 
                    }


                    // follow weibull 
                    if(useWeibull){
                        std::random_device rdt;
                        std::mt19937 gent(rdt());

                        // this is the learned distribution paeameters from the paper 
                        // https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005761
                        std::weibull_distribution<> dwt(0.44,0.6);

                        std::random_device rd2;
                        std::mt19937 g2(rd2()) ;
                        std::vector<int> tIndex(transcriptIds.size()) ;
                        std::iota(tIndex.begin(), tIndex.end(), 0);
                        std::shuffle(tIndex.begin(), tIndex.end(), g2) ;

                        for(size_t j = 0; j < transcriptIds.size() ; ++j ){
                            probVec[tIndex[j]] = dwt(gent) ;
                        }

                    }else if(useDBG){

                        dropThisGene = false ;

                        auto segCountHist = preCalculatedSegProb[i] ;

                        // if(i == geneIdToTrack  && barcode == "Cell1"){
                        // 	std::cerr << "segCountHist Size " << segCountHist.size() << "\t"
                        // 			  << "id to track " << geneIdToTrack << "\n" ;
                        // }	

                        if(segCountHist.size() == 0){
                            std::cerr << "No seg hist for this gene " << i << "\n" ; 
                              std::exit(1) ;
                              droppedGeneExpression += geneCount ;
                            dropThisGene = true ;
                            geneCount = 0 ;
                        }

                        if(!dropThisGene){
                            std::vector<size_t> segIndex(segCountHist.size()) ;
                            std::vector<uint32_t> segProbVec(segCountHist.size()) ;
                            std::vector<int> segCounts(segCountHist.size(),0) ;
                            
                            size_t ind{0} ;
                            for(auto it : segCountHist){
                                segIndex[ind] = it.first ;
                                segProbVec[ind] = it.second ;
                                ind++ ;
                            }

                            std::random_device rdt4;
                            std::mt19937 gent4(rdt4());

                            std::unordered_map<size_t, int> segCountMap ;
                            
                            std::discrete_distribution<> dmt(segProbVec.begin(), segProbVec.end()) ;
                            for(int j = 0 ; j < geneCount ; ++j){
                                ++segCounts[dmt(gent4)] ;
                            }
                            for(size_t j = 0 ; j < segCounts.size() ; ++j){
                                segCountMap[segIndex[j]] = segCounts[j] ;
                            }
                            // make sure the count is allocated
                            if((i == geneIdToTrack) && (barcode == "Cell1")){
                                std::cerr << "Tracking for gene id " << geneName << "\n"; 
                                int totNum = 0; 
                                for(auto& s : segCountMap){
                                    if(s.second != 0){
                                        totNum += s.second;
                                    }
                                }
                                if(totNum == geneCount){
                                    std::cerr << "We are good\n" ;
                                }
                            }

                            cellSegCount[cellId][i] = segCountMap ;

                        }
                    }
                    else{
                        consoleLog->error("This version takes either of --useEqClass/--weibull/--dbg") ;
                        std::exit(1) ;
                    }


                    trueGeneCounts[cellId][i] = geneCount ;
                    
                    if(useDBG){
                        continue ;	
                    }

                    auto sum_prob = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
                    if(sum_prob == 0){
                        dropThisGene = true ;
                    }

                    if (dropThisGene){
                        numOfDroppedGenes++ ;
                        continue ;
                    }
                    
                    // fill up the cell to transcriptome matrix 
                    std::random_device rdt3;
                    std::mt19937 gent3(rdt3());

                    std::discrete_distribution<> dmt(probVec.begin(), probVec.end()) ;
                    // saple from this distribution geneCount times
                    std::vector<int> transcriptCounts(transcriptIds.size(), 0) ; 
                    for(int j=0 ; j < geneCount ; ++j){
                        ++transcriptCounts[dmt(gent3)] ;
                    }
                    // total gene count 
                    // This should not be needed
                    T totCount{0} ;
                    //uint32_t lastAlevinTid{0} ;
                    for(size_t j = 0; j < transcriptIds.size(); ++j){
                        auto tid = transcriptIds[j] ;
                        auto alevinTid = alevinReverseMap[tid] ;
                        data[cellId][alevinTid] = transcriptCounts[j] ;
                        totCount += data[cellId][alevinTid] ;

                        cellSumCheck += data[cellId][alevinTid] ;
                        //lastAlevinTid = alevinTid ;
                    }
                }
            }

      //std::cerr << "\n cell sum before running " << cellSum << " cell sum after running " << cellSumCheck << "\n" ;
            cellId++ ;
            if(numOfDroppedGenes > numOfGenes - 100){
                numOfBadCells++ ;
            }
            numOfDroppedGenes = 0 ;
            _verbose("\rIn Splatter: Number of cells processed : %lu", cellId);
        
            if(useDBG && droppedGeneExpression > 0){
                std::cerr << " Num of genes dropped due to dbg " << droppedGeneExpression 
                << " out of " <<  cellSum << "\n" ;

            }
        }

    }
    if(numOfBadCells > 0)
        consoleLog->warn("{} cells have less than 100 genes present",numOfBadCells);
    
    numOfWhiteListedCells = cellWhiteListMap.size() ;
    //dumpSplatterMaps(outDir) ;
    //dump the generated cell cols 
    //which are genes 
    std::string cellColFile = outDir + "/alevin/quants_mat_cols.txt" ;
    std::ofstream cellColStream(cellColFile.c_str()) ;
    
    for(size_t splatterGeneId = 0 ; splatterGeneId < alevinGeneIndex2NameMap.size() ; ++splatterGeneId){
        cellColStream << alevinGeneIndex2NameMap[splatterGeneId] << "\n" ; 
    }
}


template<typename T>
void DataMatrix<T>::dumpGeneCountMatrix(
    std::string& fileName,
    std::unordered_set<uint32_t>& emptyCells
){
    std::ofstream geneMatrixStream(fileName) ;

    size_t numOfExpressedGenes{0} ;

    if(trueGeneCounts.size() != 0){
        for(size_t cell_id = 0 ; cell_id < trueGeneCounts.size() ; ++cell_id){
            if(emptyCells.find(cell_id) != emptyCells.end()){
                continue ;
            }
            auto v = trueGeneCounts[cell_id] ;
            for(auto c : v){
                geneMatrixStream << c << "," ;
                if(c > 0){
                    numOfExpressedGenes++ ;
                }
            }
            geneMatrixStream << "\n" ;
        }

    }else{
        for(size_t cell_id = 0 ; cell_id < geneCounts.size() ; ++cell_id){
            if(emptyCells.find(cell_id) != emptyCells.end()){
                continue ;
            }
            auto v = geneCounts[cell_id] ;
            for(auto c : v){
                geneMatrixStream << c << "," ;
                if(c > 0){
                    numOfExpressedGenes++ ;
                }
            }
            geneMatrixStream << "\n" ;
        }
    }

    //std::cerr<< "DEBUG: ==== At the end Expressed Gene count "<< numOfExpressedGenes << "\n" ;
}


template<typename T>
void DataMatrix<T>::dumpTrueCellNames(
    std::string& fileName
){
    std::ofstream cellNameStream(fileName) ;

    for(auto cn : cellNames){
        cellNameStream << cn << "\n" ;
    }

    //std::cerr<< "DEBUG: ==== At the end Expressed Gene count "<< numOfExpressedGenes << "\n" ;
}

template<typename T>
void DataMatrix<T>::dumpGeneCountMatrixUnmolested(
    std::string& fileName
){
    std::ofstream geneMatrixStream(fileName) ;
    for(auto& v : geneCounts){
        for(auto c : v){
            geneMatrixStream << c << "," ;
        }
        geneMatrixStream << "\n" ;
    }
}

template<typename T>
void DataMatrix<T>::dumpIntronCount(
    //std::string& fileName
){
    //std::ofstream geneMatrixStream(fileName) ;
    uint32_t intronCountSum{0} ;
    for(auto& v : intronCountMap){
        
        if(cellWhiteListMap["GTATTCTGTAGTGAAT"] == v.first){
            for(auto& i : v.second)
                intronCountSum += i.second ;
        }
        
    }
    std::cout << "BEFORE INTRON COUNT FOR GTATTCTGTAGTGAAT:  " << intronCountSum << "\n" ;
}

template<typename T>
DataMatrix<T>::DataMatrix(
    SimulateOptions& simOpts,
    Reference& refInfo,
    std::shared_ptr<spdlog::logger>& consoleLogIn
){

    consoleLog = consoleLogIn ;

    if (simOpts.alevinMode){
        loadAlevinData(simOpts, refInfo) ;
    }else if(simOpts.splatterMode || simOpts.normalMode){
        loadSplatterData(simOpts, refInfo) ;
    }else{
        consoleLog->info("Please specify either of --splatter-mode or --alevin-mode ") ;
    }
}

template class DataMatrix<int> ;
template class DataMatrix<double> ;


//template 
//void pupulateGeneCountFromTxt<int>(
//	std::string& countFile, 
//	std::vector<std::vector<int>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes	
//); 
//template 
//void pupulateGeneCountFromTxt<double>(
//	std::string& countFile, 
//	std::vector<std::vector<double>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes	
//); 
//
//template
//void populateGeneCountFromBinary<int>(
//	std::string& countFileBinary, 
//	std::vector<std::vector<int>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes
//);
//
//template
//void populateGeneCountFromBinary<double>(
//	std::string& countFileBinary, 
//	std::vector<std::vector<double>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes
//);
