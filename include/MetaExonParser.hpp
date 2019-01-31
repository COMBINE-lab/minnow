#ifndef BFH_CLASS_HPP
#define BFH_CLASS_HPP

#include <cstdio>
#include <cmath>
#include <random>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include "MinnowFS.hpp"

#include "string_view.hpp"
#include "ReferenceInfo.hpp"
#include "MinnowUtil.hpp"
#include "macros.hpp"

class BFHClass{
    public:
    BFHClass(
        std::string& bfhFileIn
    ){
        bfhFile = bfhFileIn ;
    }

    void loadBFH(
        std::string& bfhFile,
        std::string& cellClusterFile,
        Reference& refInfo,
        std::map<std::string, uint32_t>& cellWhiteListMap,
        bool generateNoiseProfile,
        std::unordered_map<std::string, uint32_t>& cellNoisyMap
    ) ;
    void dumpClusterHistoGram(std::string& file_name) ;
    
    inline uint32_t getGeneLevelProbCount(uint32_t geneId, uint32_t count){
        return geneCountHistogram[geneId][count] ;
    }


    std::string bfhFile  ;
    std::vector<double> countProbability ;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> geneCountHistogram ; // Gene id -> (EqClass_Length -> Numebr)
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> noisyGeneCountHistogram ; // Gene id -> (EqClass_Length -> Numebr)
    std::unordered_map<uint32_t,
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>> clusterCountHistogram ; // Cluster id -> (Gene id -> (EqClass Length -> Number))
    std::unordered_map<uint32_t, uint32_t> cell2ClusterMap ;

};

#endif