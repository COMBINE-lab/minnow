#include "BFHClass.hpp"
#include "MinnowUtil.hpp"
#include <string>
#include <cstdlib>
#include <limits>

// <gene name> # number of keys <exon> <num> <transcripts> 

void BFHClass::dumpClusterHistoGram(std::string& file_name){
    std::ofstream histStream(file_name.c_str()) ;
    histStream << clusterCountHistogram.size() << "\n" ;
    for(auto& it : clusterCountHistogram){
        auto cluster_id = it.first ;
        auto size_of_map = it.second.size() ;
        histStream << cluster_id << "\n" ;
        // Number of Genes 
        histStream << size_of_map << "\n" ;
        for(auto it2 : it.second){
            // Gene id 
            histStream << it2.first << "\t" ;
            // Histogram for that gene id 
            for(auto it3 : it2.second){
                histStream << it3.first << "\t" << it3.second << "\t" ;
            }
            histStream << "\n" ;
        }
        //histStream << "\n" ;
    }
    histStream.close() ;
}


void BFHClass::loadBFHLite(
    std::string& bfhFile,
    Reference& refInfo,
    std::string& outDir
){
   if(! util::fs::FileExists(bfhFile.c_str())){
      consoleLog->error("{} does not exists", bfhFile) ;
      std::exit(1) ;
	}

    consoleLog->info("Reading BFH file {}",bfhFile);
	std::ifstream dataStream(bfhFile.c_str()) ;
	std::string line ;

    std::getline(dataStream, line) ;
    uint32_t numTranscripts = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numCells = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numEqClasses = std::stoul(line) ;

    consoleLog->info("numTranscripts: {}", numTranscripts) ;
    consoleLog->info("numCells: {}", numCells) ;
    consoleLog->info("numEqClasses: {}", numEqClasses) ;

    std::vector<std::string> trNames(numTranscripts) ;
    std::vector<std::string> CBNames(numCells) ;
    for(size_t i  = 0; i < numTranscripts; ++i){
        std::getline(dataStream, line) ;
        trNames[i] = line ;
    }

    consoleLog->info("Transcripts read ") ;
    for(size_t i  = 0; i < numCells; ++i){
        std::getline(dataStream, line) ;
        CBNames[i] = line ;
    }


    consoleLog->info("Cell names read ") ;

    // read equivalence classes now
    uint32_t tot_reads{0} ;
    std::unordered_map<uint32_t, uint32_t> countHistogram ;
    auto& transcript2geneMap = refInfo.transcript2geneMap ;

    while(std::getline(dataStream, line)){
        std::vector<std::string> tokens ;
        util::split(line, tokens, "\t") ;
        uint32_t num_labels = std::stoul(tokens[0]) ;
        std::unordered_set<uint32_t> geneIds ;
        for(uint32_t i = 0 ; i < num_labels; ++i){
            auto tname = trNames[std::stoul(tokens[1+i])] ;
            auto it = refInfo.transcriptNameMap.find(tname) ;
            if(it != refInfo.transcriptNameMap.end()){
                auto tid = it->second ;
                if(transcript2geneMap.find(tid) != transcript2geneMap.end()){
                  auto gid = transcript2geneMap[tid] ;
                  geneIds.insert(gid) ;
                }else{
                  consoleLog->error("transcript is in the list but no corresponding gene found "
                                    "this signifies that the BFH and the annotation does not belong "
                                    "to the same annotation (e.g. gencode version) or same organism"
                  ) ;
                  std::exit(2) ;
                }
            }
        }

        uint32_t num_reads = std::stoul(tokens[num_labels+1]) ;

        for(auto gid : geneIds){
          geneCountHistogram[gid][num_labels] += num_reads ;
        }

        auto it = countHistogram.find(num_labels) ;
        if(it != countHistogram.end()){
            countHistogram[num_labels] += num_reads ;
        }else{
            countHistogram[num_labels] = num_reads ;
        }
        tot_reads += num_reads ;
    } 

    consoleLog->info("countHistogram.size(): {}", countHistogram.size());

    auto x = std::max_element( countHistogram.begin(), countHistogram.end(),
        [](const std::pair<uint32_t, uint32_t>& p1, const std::pair<uint32_t, uint32_t>& p2) {
        return p1.first < p2.first; });

    consoleLog->info("Histogram statistics "
                     "max value: {}, max freq: {}, total reads {}",
                     x->first, x->second, tot_reads
    );

    consoleLog->info("Converting histogram to a probablity vector");
    countProbability.resize(x->first + 1, 0.0) ;
    for(auto it: countHistogram){
      if(it.first >= countProbability.size()){
        consoleLog->error("[DEBUG] out of memory {} >= {}",it.first,countProbability.size()) ;
      }
      countProbability[it.first] = static_cast<double>(it.second)/static_cast<double>(tot_reads) ; 
    }

    {
        std::string geneCountHistogramFile = outDir + "/geneLevelProb.txt" ;
        std::cerr << "DEBUG:  " << geneCountHistogramFile << "\n" ;
        std::ofstream probStream(geneCountHistogramFile.c_str()) ;
        probStream << geneCountHistogram.size() << "\n" ;
        std::cerr << "DEBUG: " << geneCountHistogram.size() << " " << geneCountHistogramFile << "\n" << std::flush;
        for(auto it : refInfo.geneMap){
            if(geneCountHistogram.find(it.second) != geneCountHistogram.end()){
                // print gene name
                probStream << it.first << "\n" ;
                probStream << geneCountHistogram[it.second].size() << "\n" ;
                for(auto it2: geneCountHistogram[it.second]){
                    probStream << it2.first << "\t" << it2.second << "\n" ;
                }
            }
        }

    }

    {
        std::string countProbabilityFile = outDir + "/countProb.txt" ;
        std::ofstream probStream(countProbabilityFile.c_str()) ;
        probStream << countProbability.size() << "\n" ;
        for(auto prob : countProbability){
            // print gene name
            probStream << prob << "\n" ;
        }
    }

    consoleLog->info("Exiting after reading BFH ") ; 

}


void BFHClass::loadBFH(
    std::string& bfhFile,
    std::string& cellClusterFile,
    Reference& refInfo,
    std::map<std::string, uint32_t>& cellWhiteListMap,
    bool generateNoiseProfile,
    std::unordered_map<std::string, uint32_t>& cellNoisyMap,
    std::string& outDir,
    bool dump
){



    if(! util::fs::FileExists(bfhFile.c_str())){
      consoleLog->error("{} does not exists", bfhFile) ;
      std::exit(1) ;
	  }

    bool createClusterLevelHist{false} ;
    if(cellClusterFile != ""){
        consoleLog->info("Feeding cell cluster file {}", cellClusterFile) ;
        if(! util::fs::FileExists(cellClusterFile.c_str())){
            consoleLog->error("Cell cluster file {} does not exist", cellClusterFile) ;
            std::exit(1) ;
        }
        createClusterLevelHist = true ;
    }

    if(createClusterLevelHist){
        std::ifstream clusterStream(cellClusterFile.c_str()) ;
        std::string line ;

        while(std::getline(clusterStream, line)){
            std::vector<std::string> tokens ;
            util::split(line, tokens, "\t") ;
            auto cell_name = tokens[0] ;
            uint32_t cluster_id = std::stoul(tokens[1]) ;
            auto it = cellWhiteListMap.find(cell_name);
            if(it != cellWhiteListMap.end()){
                auto cell_id = it->second ;
                cell2ClusterMap[cell_id] = cluster_id ;
            }else{
                consoleLog->error("Avoiding {} cluster", cluster_id) ;
            }
        }
        consoleLog->info("read cluster file with size {}",cell2ClusterMap.size());
    }

    consoleLog->info("Reading BFH file {}",bfhFile);
	std::ifstream dataStream(bfhFile.c_str()) ;
	std::string line ;

    std::getline(dataStream, line) ;
    uint32_t numTranscripts = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numCells = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numEqClasses = std::stoul(line) ;

    consoleLog->info("numTranscripts: {}", numTranscripts) ;
    consoleLog->info("numCells: {}", numCells) ;
    consoleLog->info("numEqClasses: {}", numEqClasses) ;

    std::vector<std::string> trNames(numTranscripts) ;
    std::vector<std::string> CBNames(numCells) ;
    for(size_t i  = 0; i < numTranscripts; ++i){
        std::getline(dataStream, line) ;
        trNames[i] = line ;
    }

    consoleLog->info("Transcripts read ") ;
    for(size_t i  = 0; i < numCells; ++i){
        std::getline(dataStream, line) ;
        CBNames[i] = line ;
    }


    consoleLog->info("Cell names read ") ;

    // read equivalence classes now
    uint32_t tot_reads{0} ;
    std::unordered_map<uint32_t, uint32_t> countHistogram ;
    auto& transcript2geneMap = refInfo.transcript2geneMap ;

    while(std::getline(dataStream, line)){
        std::vector<std::string> tokens ;
        util::split(line, tokens, "\t") ;
        uint32_t num_labels = std::stoul(tokens[0]) ;
        std::unordered_set<uint32_t> geneIds ;
        for(uint32_t i = 0 ; i < num_labels; ++i){
            auto tname = trNames[std::stoul(tokens[1+i])] ;
            auto it = refInfo.transcriptNameMap.find(tname) ;
            if(it != refInfo.transcriptNameMap.end()){
                auto tid = it->second ;
                if(transcript2geneMap.find(tid) != transcript2geneMap.end()){
                  auto gid = transcript2geneMap[tid] ;
                  geneIds.insert(gid) ;
                }else{
                  consoleLog->error("transcript is in the list but no corresponding gene found "
                                    "this signifies that the BFH and the annotation does not belong "
                                    "to the same annotation (e.g. gencode version) or same organism"
                  ) ;
                  std::exit(2) ;
                }
            }
        }

        uint32_t num_reads = std::stoul(tokens[num_labels+1]) ;

        for(auto gid : geneIds){
          geneCountHistogram[gid][num_labels] += num_reads ;
        }

        if(createClusterLevelHist){
            uint32_t numOfCells = std::stoul(tokens[num_labels+2]) ;
            size_t idx = num_labels + 3 ; // index of first CB
            //std::cerr << "Size of tokens " << tokens.size() << " idx " << idx << "\n"
            //          << " numCells " << numCells << "\n" 
            //          << " line "  << line << "\n" ;

            for(uint32_t i = 0; i < numOfCells ; ++i){
                size_t internalCellBarcodeId = std::stoul(tokens[idx]) ;

                // cell specific information
                size_t numOfUmis = std::stoul(tokens[idx+1]) ;
                uint32_t umiCountSum{0} ;

                for(size_t umi_id = 0 ; umi_id < numOfUmis ; ++umi_id){
                    size_t umi_cnt_idx = idx + 1 + ((umi_id * 2) + 2) ;
                    //std::cerr << " umi_cnt_idx " << umi_cnt_idx << "\n" ;
                    umiCountSum += std::stoul(tokens[umi_cnt_idx]) ;
                }

                std::string cell_name = CBNames[internalCellBarcodeId] ;
                // for now concentrate on cell id that are whitelisted
                auto it = cellWhiteListMap.find(cell_name) ;
                if(it != cellWhiteListMap.end()){
                    auto globalCellBarcodeId = it->second ;
                    auto cluster_id = cell2ClusterMap[globalCellBarcodeId] ;

                    for(auto gid : geneIds){
                        clusterCountHistogram[cluster_id][gid][num_labels] += umiCountSum ;
                    }

                }else if(generateNoiseProfile){
                    auto it2 = cellNoisyMap.find(cell_name) ;
                    if(it2 != cellNoisyMap.end()){
                        for(auto gid : geneIds){
                            noisyGeneCountHistogram[gid][num_labels] += umiCountSum ;
                        }
                    }
                }
                // advance idx 
                idx += 1 + std::stoul(tokens[idx+1]) * 2 + 1 ;
            }
        }
 
        auto it = countHistogram.find(num_labels) ;
        if(it != countHistogram.end()){
            countHistogram[num_labels] += num_reads ;
        }else{
            countHistogram[num_labels] = num_reads ;
        }
        tot_reads += num_reads ;
    }


    consoleLog->info("countHistogram.size(): {}", countHistogram.size());



    auto x = std::max_element( countHistogram.begin(), countHistogram.end(),
        [](const std::pair<uint32_t, uint32_t>& p1, const std::pair<uint32_t, uint32_t>& p2) {
        return p1.first < p2.first; });

    consoleLog->info("Histogram statistics "
                     "max value: {}, max freq: {}, total reads {}",
                     x->first, x->second, tot_reads
    );

    consoleLog->info("Converting histogram to a probablity vector");
    countProbability.resize(x->first + 1, 0.0) ;
    for(auto it: countHistogram){
      if(it.first >= countProbability.size()){
        consoleLog->error("[DEBUG] out of memory {} >= {}",it.first,countProbability.size()) ;
      }

        countProbability[it.first] = static_cast<double>(it.second)/static_cast<double>(tot_reads) ; 
    }


    if(dump)
    {
        std::string geneCountHistogramFile = outDir + "/geneLevelProb.txt" ;
        std::cerr << "DEBUG:  " << geneCountHistogramFile << "\n" ;
        std::ofstream probStream(geneCountHistogramFile.c_str()) ;
        probStream << geneCountHistogram.size() << "\n" ;
        std::cerr << "DEBUG: " << geneCountHistogram.size() << " " << geneCountHistogramFile << "\n" << std::flush;
        for(auto it : refInfo.geneMap){
            if(geneCountHistogram.find(it.second) != geneCountHistogram.end()){
                // print gene name
                probStream << it.first << "\n" ;
                probStream << geneCountHistogram[it.second].size() << "\n" ;
                for(auto it2: geneCountHistogram[it.second]){
                    probStream << it2.first << "\t" << it2.second << "\n" ;
                }
            }
        }

    }

    if(dump)
    {
        std::string countProbabilityFile = outDir + "/countProb.txt" ;
        std::ofstream probStream(countProbabilityFile.c_str()) ;
        probStream << countProbability.size() << "\n" ;
        for(auto prob : countProbability){
            // print gene name
            probStream << prob << "\n" ;
        }
    }

    consoleLog->info("Exiting after reading BFH ") ;

}

void BFHClass::loadProbability(std::string& file, Reference& refInfo, bool geneLevel){
    std::ifstream fileStream(file.c_str()) ;
    if(geneLevel){
        std::string line ;
        std::getline(fileStream, line) ;
        size_t numOfGenes = std::stoul(line) ;

        for(size_t i = 0 ; i < numOfGenes ; ++i){
            std::getline(fileStream, line) ;
            std::string geneName = line ;
            uint32_t geneId{0} ;
            bool skipThisGene{false} ;
            if(refInfo.geneMap.find(geneName) != refInfo.geneMap.end()){
                geneId = refInfo.geneMap[geneName] ;
            }else{
                skipThisGene = true ;
            }

            std::getline(fileStream, line) ;
            size_t numOfEntries = std::stoul(line) ;
            for(size_t j =0 ; j < numOfEntries; ++j){
                std::getline(fileStream, line) ;
                
                if(skipThisGene)
                    continue ;

                std::vector<std::string> tokens ;
                util::split(line, tokens, "\t") ;
                if(tokens.size()!= 2){
                    std::cerr << "FATAL ERROR!!! The probability file is ill-formed\n" ;
                    std::exit(1) ;
                }
                uint32_t value = std::stoul(tokens[0]) ;
                uint32_t count = std::stoul(tokens[1]) ;

                geneCountHistogram[geneId][value] = count ;
                 
            }
        }
    }else{
        std::string line ;
        std::getline(fileStream, line) ;
        size_t numOfEntries = std::stoul(line) ;
        countProbability.resize(numOfEntries) ;
        for(size_t i = 0 ; i < numOfEntries; ++i){
            std::getline(fileStream, line) ;
            countProbability[i] = std::stod(line) ;
        }
    }
}
