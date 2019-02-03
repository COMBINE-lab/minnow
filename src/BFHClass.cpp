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


void BFHClass::loadBFH(
    std::string& bfhFile,
    std::string& cellClusterFile,
    Reference& refInfo,
    std::map<std::string, uint32_t>& cellWhiteListMap,
    bool generateNoiseProfile,
    std::unordered_map<std::string, uint32_t>& cellNoisyMap
){
    
    if(! util::fs::FileExists(bfhFile.c_str())){
		std::cerr << bfhFile << " does not exists \n" ;
		std::exit(1) ; 
	}

    bool createClusterLevelHist{false} ;
    std::cout << "cell Clust file " << cellClusterFile << "\n" ;
    if(cellClusterFile != ""){
        if(! util::fs::FileExists(cellClusterFile.c_str())){
            std::cerr << cellClusterFile << " is not empty and does not exist \n" ;
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
                std::cerr << "Avoiding this cluster\n" ;
            }
        }
    }
    std::cerr << "DEBUG : read cluster file with size " << cell2ClusterMap.size() << "\n" ;

	std::ifstream dataStream(bfhFile.c_str()) ;
	std::string line ; 

    std::getline(dataStream, line) ;
    uint32_t numTranscripts = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numCells = std::stoul(line) ;
    std::getline(dataStream, line) ;
    uint32_t numEqClasses = std::stoul(line) ;

    std::cerr << "DEBUG: numEqClasses: " << numEqClasses << "\n" ;

    std::vector<std::string> trNames(numTranscripts) ;
    std::vector<std::string> CBNames(numCells) ;
    for(size_t i  = 0; i < numTranscripts; ++i){
        std::getline(dataStream, line) ;
        trNames[i] = line ;
    }
    
    for(size_t i  = 0; i < numCells; ++i){
        std::getline(dataStream, line) ;
        CBNames[i] = line ;
    }

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
                auto gid = transcript2geneMap[tid] ;
                geneIds.insert(gid) ;
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

    

    auto x = std::max_element( countHistogram.begin(), countHistogram.end(),
        [](const std::pair<uint32_t, uint32_t>& p1, const std::pair<uint32_t, uint32_t>& p2) {
        return p1.first < p2.first; });

    countProbability.resize(x->first, 0.0) ;
    for(auto it: countHistogram){
        countProbability[it.first] = static_cast<double>(it.second)/static_cast<double>(tot_reads) ; 
    }
    
}
