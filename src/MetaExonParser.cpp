#include "MetaExonParser.hpp"
#include "MinnowUtil.hpp"
#include <string>
#include <cstdlib>
#include <limits>

// <gene name> # number of keys <exon> <num> <transcripts> 

void ExonEqClass::loadEqClassInfo(Reference& refInfo){
    
    {

        if(! util::fs::FileExists(eqClassCountFile.c_str())){
            std::cerr << eqClassCountFile << " does not exists \n" ;
            std::exit(1) ; 
        }
        std::ifstream dataStream(eqClassCountFile.c_str()) ;
        std::string line ;
        
        std::getline(dataStream, line) ;
        uint32_t numOfEqClasses = std::stoul(line) ;
        eqCountVec.resize(numOfEqClasses) ;
        while(std::getline(dataStream, line)){
            std::vector<std::string> tokens ;
            util::split(line, tokens, "\t") ;
            if (tokens.size() == 3){
                uint32_t id = std::stoul(tokens[0]) ;
                uint32_t count = std::stoul(tokens[2]) ;
                eqCountVec[id] = count ;
            }else{
                std::cerr << eqClassCountFile << " is not formatted properly should contain 3 fields within a file\n" ;
                std::exit(1) ;
            }
        }
    }

    std::cerr << "DEBUG: ==== in ExonEqClass::loadEqClassInfo eqCountVec.size() =  "<< eqCountVec.size() << "\n" ; 
    {
        if(! util::fs::FileExists(tr2EqClassFile.c_str())){
            std::cerr << tr2EqClassFile << " does not exists \n" ;
            std::exit(1) ; 
        }
        size_t skippedExonJuntions{0} ;

        std::ifstream dataStream(tr2EqClassFile.c_str()) ;
        std::string line ;
        uint32_t numOfTranscriptsSkipped{0} ;

        while(std::getline(dataStream, line)){
            std::vector<std::string> tokens ;
            util::split(line, tokens, "\t") ;
            // check sanity 
            if(tokens.size()%4 != 2){
                std::cerr << eqClassCountFile << " is not formatted properly should contain 3 fields within a file\n" ;
                std::exit(1) ;
            }

            uint32_t tid{std::numeric_limits<uint32_t>::max()} ;
            uint32_t refLength{0} ;
            auto trName = tokens[0] ;
            auto it = refInfo.transcriptNameMap.find(trName) ;
            bool atleastOneAdded{false} ;
            if(it != refInfo.transcriptNameMap.end()){
                // this will be index to the reference
                // not minnow matrix
                tid = it->second ;
                refLength = refInfo.transcripts[tid].RefLength ;
            }else{
                //std::cerr << "Trancript not present while reading equivalence classes\n" ;
                numOfTranscriptsSkipped++ ;
            }
            std::vector<uint32_t> numTokens ; 
            std::transform(tokens.begin()+1, tokens.end(), std::back_inserter(numTokens), [](const std::string& val)
    		{
    			return std::stoul(val);
    		});

            if(tid != std::numeric_limits<uint32_t>::max()){

                uint32_t numOfClasses = numTokens[0] ;
                //std::cerr << "Number of classes for " << trName << "\t" << numOfClasses << "\n" ; 
                for(uint32_t i = 1; i < numOfClasses*4; i+=4){
                    // Decide if you want to consider this exon junction or not 
                    // If the distance from start of Exon_1 is less then READ_LEN
                    // then we ignore that equivalence class
                    
                    uint32_t farthestEnd = std::max(numTokens[i+2], numTokens[i+3]) ; 
                    int distance = static_cast<int>(farthestEnd - numTokens[i+1] + 1) ;
                    if(distance >= READ_LEN){
                        transcript2EqMap[tid][numTokens[i]] =
                            {
                                numTokens[i],
                                numTokens[i+1],
                                numTokens[i+2],
                                numTokens[i+3]
                            } ;
                        atleastOneAdded = true ;
                    }else{
                        skippedExonJuntions++ ;
                    }
                }
            }
            //if(!atleastOneAdded){
            //    std::cerr << " None of this transcript " << refInfo.transcripts[tid].RefName << " added \t"  
            //              << " Reference length " << refInfo.transcripts[tid].RefLength 
            //              << " line in file " << line << "\n" ;       

            //    //std::exit(1) ;
            //}
        }

       std::cerr << "Num of transcripts skipped " << numOfTranscriptsSkipped << "\n" ;
       std::cerr << "DEBUG: ==== skippedExonJuntions " << skippedExonJuntions 
                 << " skipped transcripts as a whole " << refInfo.transcripts.size() - transcript2EqMap.size()
                 << " Out of " << refInfo.transcripts.size() << "\n" ;        

    }

    std::cerr << "DEBUG: ==== in ExonEqClass::loadEqClassInfo transcript2EqMap.size() =  "<< transcript2EqMap.size() << "\n" ; 
}


void ExonEqClass::dumpClusterHistoGram(std::string& file_name){
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

void ExonEqClass::getTrLevelDistribution(
    uint32_t& cluster_id,
    uint32_t& geneId,
    std::vector<uint32_t>& tids,
    std::unordered_map<uint32_t, uint32_t>& exonTrMap,
    std::unordered_map<uint32_t, std::vector<double>>& clusterProbVec,
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& clusterExonLevelProb
){
    std::unordered_map<uint32_t, uint32_t> trExMap ;
    //probVec.resize(tids.size(),0.0) ;
    std::unordered_set<uint32_t> eqClassIds ;
    for(auto tid : tids){
        auto exonList = transcript2EqMap[tid] ;
        for(auto ex : exonList){
            if (exonTrMap.find(ex.first) == exonTrMap.end()){
                exonTrMap[ex.first] = tid ;
                trExMap[tid] = ex.first ;
            }
            eqClassIds.insert(ex.first) ;
        }
    }
    std::unordered_map<uint32_t,std::unordered_map<uint32_t, double>> clusterTrLevelProb ;

    double totProb{0.0} ;

    for(auto exonId : eqClassIds){
        auto count = eqCountVec[exonId] ; 
        clusterExonLevelProb[cluster_id][exonId] = static_cast<double>(clusterCountHistogram[cluster_id][geneId][count]) ;
        totProb += clusterExonLevelProb[cluster_id][exonId] ; 
    }

    for(auto it2: clusterExonLevelProb[cluster_id]){
        auto tid = exonTrMap[it2.first] ;
        clusterTrLevelProb[cluster_id][tid] += it2.second / totProb ;
    }


    clusterProbVec[cluster_id] = std::vector<double>(tids.size(), 0.0) ;
    for(size_t i = 0 ; i < tids.size(); ++i){
        clusterProbVec[cluster_id][i] = clusterTrLevelProb[cluster_id][tids[i]] ;
    }
    //std::cerr << "In getTrLevelDist " << cluster_id << "\n" ;

}

void ExonEqClass::getTrLevelDistribution(
    uint32_t& geneId,
    std::vector<uint32_t>& tids,
    std::unordered_map<uint32_t, uint32_t>& exonTrMap,
    std::unordered_map<uint32_t, std::vector<double>>& clusterProbVec,
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& clusterExonLevelProb
){

    //std::cerr << "GET TR level distribution \n" ;

    std::unordered_map<uint32_t, uint32_t> trExMap ;
    //probVec.resize(tids.size(),0.0) ;
    std::unordered_set<uint32_t> eqClassIds ;
    for(auto tid : tids){
        auto exonList = transcript2EqMap[tid] ;
        for(auto ex : exonList){
            if (exonTrMap.find(ex.first) == exonTrMap.end()){
                exonTrMap[ex.first] = tid ;
                trExMap[tid] = ex.first ;
            }
            eqClassIds.insert(ex.first) ;
        }
    }
    std::unordered_map<uint32_t,std::unordered_map<uint32_t, double>> clusterTrLevelProb ;

    for(auto it : cell2ClusterMap){
        auto cluster_id = it.second ;
        double totProb{0.0} ;

        for(auto exonId : eqClassIds){
            auto count = eqCountVec[exonId] ; 
            clusterExonLevelProb[cluster_id][exonId] = static_cast<double>(clusterCountHistogram[cluster_id][geneId][count]) ;
            totProb += clusterExonLevelProb[cluster_id][exonId] ; 
        }

        for(auto it2: clusterExonLevelProb[cluster_id]){
            auto tid = exonTrMap[it2.first] ;
            clusterTrLevelProb[cluster_id][tid] += it2.second / totProb ;
        }


        clusterProbVec[cluster_id] = std::vector<double>(tids.size(), 0.0) ;
        for(size_t i = 0 ; i < tids.size(); ++i){
            clusterProbVec[cluster_id][i] = clusterTrLevelProb[cluster_id][tids[i]] ;
        }
    }


    //totProb = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
    //if (totProb > 0.5){
    //    std::cout << "\n Sum " << totProb << "\n" ;
    //    for(auto eit: exonLevelProb){
    //        std::cerr << "ExonId " << eit.first 
    //        << "\ttid " << exonTrMap[eit.first] 
    //        << "\t" << " Prob " << eit.second << "\n" ;
    //    }
    //    std::exit(1) ;
    //}
    //std::cerr << "DONE WITH  TR level distribution \n" ;
    
}


void ExonEqClass::getNoisyTrLevelDistribution(
    uint32_t& geneId,
    std::vector<uint32_t>& tids,
    std::unordered_map<uint32_t, uint32_t>& exonTrMap,
    std::vector<double>& probVec,
    std::unordered_map<uint32_t, double>& exonLevelProb
){

    //std::cerr << "GET TR level distribution \n" ;

    std::unordered_map<uint32_t, uint32_t> trExMap ;
    probVec.resize(tids.size(),0.0) ;
    std::unordered_set<uint32_t> eqClassIds ;
    for(auto tid : tids){
        auto exonList = transcript2EqMap[tid] ;
        for(auto ex : exonList){
            if (exonTrMap.find(ex.first) == exonTrMap.end()){
                exonTrMap[ex.first] = tid ;
                trExMap[tid] = ex.first ;
            }
            eqClassIds.insert(ex.first) ;
        }
    }
    std::unordered_map<uint32_t, double> trLevelProb ;
    double totProb{0.0} ;

    for(auto exonId : eqClassIds){
        auto count = eqCountVec[exonId] ; 
        exonLevelProb[exonId] = static_cast<double>(noisyGeneCountHistogram[geneId][count]) ;
        totProb += exonLevelProb[exonId] ; 

    }

    for(auto it: exonLevelProb){
        auto tid = exonTrMap[it.first] ;
        if (trLevelProb.find(tid) != trLevelProb.end()){
            trLevelProb[tid] += it.second / totProb ;
        }else{
            trLevelProb[tid] = it.second / totProb ; 
        }
    }

    for(size_t i = 0 ; i < tids.size(); ++i){
        if (trLevelProb.find(tids[i]) != trLevelProb.end()){
            probVec[i] = trLevelProb[tids[i]] ;
        }  
    }
    //totProb = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
    //if (totProb > 0.5){
    //    std::cout << "\n Sum " << totProb << "\n" ;
    //    for(auto eit: exonLevelProb){
    //        std::cerr << "ExonId " << eit.first 
    //        << "\ttid " << exonTrMap[eit.first] 
    //        << "\t" << " Prob " << eit.second << "\n" ;
    //    }
    //    std::exit(1) ;
    //}
    //std::cerr << "DONE WITH  TR level distribution \n" ;
    
}

void ExonEqClass::getTrLevelDistribution(
    uint32_t& geneId,
    std::vector<uint32_t>& tids,
    std::unordered_map<uint32_t, uint32_t>& exonTrMap,
    std::vector<double>& probVec,
    std::unordered_map<uint32_t, double>& exonLevelProb,
    uint32_t& numCells
){

    //std::cerr << "GET TR level distribution \n" ;

    std::unordered_map<uint32_t, uint32_t> trExMap ;
    probVec.resize(tids.size(),0.0) ;
    std::unordered_set<uint32_t> eqClassIds ;
    for(auto tid : tids){
        auto exonList = transcript2EqMap[tid] ;
        for(auto ex : exonList){
            if (exonTrMap.find(ex.first) == exonTrMap.end()){
                exonTrMap[ex.first] = tid ;
                trExMap[tid] = ex.first ;
            }
            eqClassIds.insert(ex.first) ;
        }
    }
    std::unordered_map<uint32_t, double> trLevelProb ;
    double totProb{0.0} ;

    for(auto exonId : eqClassIds){
        auto count = eqCountVec[exonId] ; 
        if(numCells < 10){
            exonLevelProb[exonId] = countProbability[count] ;
            totProb += countProbability[count] ; 

        }else{
            exonLevelProb[exonId] = static_cast<double>(geneCountHistogram[geneId][count]) ;
            totProb += exonLevelProb[exonId] ; 

        }
    }

    for(auto it: exonLevelProb){
        auto tid = exonTrMap[it.first] ;
        if (trLevelProb.find(tid) != trLevelProb.end()){
            trLevelProb[tid] += it.second / totProb ;
        }else{
            trLevelProb[tid] = it.second / totProb ; 
        }
    }

    for(size_t i = 0 ; i < tids.size(); ++i){
        if (trLevelProb.find(tids[i]) != trLevelProb.end()){
            probVec[i] = trLevelProb[tids[i]] ;
        }  
    }
    //totProb = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
    //if (totProb > 0.5){
    //    std::cout << "\n Sum " << totProb << "\n" ;
    //    for(auto eit: exonLevelProb){
    //        std::cerr << "ExonId " << eit.first 
    //        << "\ttid " << exonTrMap[eit.first] 
    //        << "\t" << " Prob " << eit.second << "\n" ;
    //    }
    //    std::exit(1) ;
    //}
    //std::cerr << "DONE WITH  TR level distribution \n" ;
    
}


void ExonEqClass::loadBFH(
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

                }else{
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


void ExonEqClass2::fillGeneLevelMap(
    Reference& refInfo
){
    if(! util::fs::FileExists(fname_.c_str())){
		std::cerr << fname_ << " does not exists \n" ;
		std::exit(1) ; 
	}
	std::ifstream dataStream(fname_.c_str()) ;
	std::string line ; 

	size_t geneId{0} ;
	while(std::getline(dataStream, line)){
		geneId += 1 ;

		std::vector<std::string> tokens ; 
		util::split(line, tokens, "\t") ;

        std::string geneName = tokens[0] ;
        int numKeys = std::stoi(tokens[1]) ;

        size_t ind = 2 ;
        for(int i = 0; i < numKeys; ++i){
            std::string exon_name = tokens[ind++] ;
            int trCount = std::stoi(tokens[ind++]) ;
            size_t j = 0 ;
            std::vector<uint32_t> trList(trCount) ;

            while(j < trCount){
                std::string trName = tokens[ind++] ;
                auto it = refInfo.transcriptNameMap.find(trName) ;
                if(it != refInfo.transcriptNameMap.end()){
                    // this will be index to the reference
                    // not minnow matrix
                    auto trId = it->second ;
                    trList[j] = trId ;
                }else{
                    std::cerr << "Trancript not present while reading equivalence classes\n" ;
                }
                ++j ;
            }
            // The tr list is completely read Now 
            // we can update the map 
            geneLevelMap[geneName].emplace_back(
                exon_name,
                trList,
                trList.size()
            ); 
        }


	}
}


void ExonEqClass2::loadDistributionVector(
    std::string& geneName,
    std::unordered_map<uint32_t, size_t>& trIdMap,
    std::vector<double> probVec  
){
    if(geneLevelMap.size() == 0){
        std::cerr << "this function should be called after fillGeneLevelmap " ;
    }

    auto it = geneLevelMap.find(geneName) ;
    double totMass{0} ;
    probVec.resize(trIdMap.size(), 0.0) ;

    if(it != geneLevelMap.end()){
        std::vector<ExonLevel> exEqClass = it->second ;
        for(auto ex : exEqClass){
            auto tid = ex.trList[ex.labelLength - 1] ;
            double prob = countProbability[ex.labelLength] ;
            auto it = trIdMap.find(tid) ;
            if(it != trIdMap.end()){
                probVec[trIdMap[tid]] += prob ;
            }else{
                std::cerr << "trList.Size " << ex.trList.size() << "\n" ;
                std::cerr << geneName << "\n" ;
                for(auto tr: ex.trList){
                    std::cerr << tr << "\n" ;
                }
                std::cerr << "trIdMap size " << trIdMap.size() << "\n" ;
                for(auto it2 : trIdMap){
                    std::cerr << it2.first << "\t" << it2.second << "\n" ;
                }
                std::cerr << "The exon belongs to a transcript outside gene \n" ;
                std::exit(1) ; 
            }
            totMass += prob ;
        }
    }
    for(size_t i = 0; i < probVec.size() ; i++){
        probVec[i] = probVec[i] / totMass ;
    }

}