#ifndef MATRIX_EXON_HPP
#define MATRIX_EXON_HPP


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

class ExonEqClass{
    public:
    ExonEqClass(
        std::string& exon_info_dir_in
    ){
        exon_info_dir = exon_info_dir_in ;
        eqClassCountFile = exon_info_dir + "/eqclass_count.txt" ;
        tr2EqClassFile = exon_info_dir + "/transcript_to_eqclass.txt" ;
    }

    

    void loadEqClassInfo(Reference & refInfo) ;
    //void loadTranscriptEqMap(Reference& refInfo) ;
    void getTrLevelDistribution(
        uint32_t& geneId,
        std::vector<uint32_t>& tids,
        std::unordered_map<uint32_t, uint32_t>& exonTrMap,
        std::vector<double>& probVec,
        std::unordered_map<uint32_t, double>& exonLevelProb,
        uint32_t& numCells
    ) ;

    void getTrLevelDistribution(
        uint32_t& geneId,
        std::vector<uint32_t>& tids,
        std::unordered_map<uint32_t, uint32_t>& exonTrMap,
        std::unordered_map<uint32_t, std::vector<double>>& clusterProbVec,
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& clusterExonLevelProb
    ) ;
    //void loadBFH(std::string& bfhFileName, Reference& refInfo) ;
    
    void getTrLevelDistribution(
        uint32_t& cluster_id,
        uint32_t& geneId,
        std::vector<uint32_t>& tids,
        std::unordered_map<uint32_t, uint32_t>& exonTrMap,
        std::unordered_map<uint32_t, std::vector<double>>& clusterProbVec,
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& clusterExonLevelProb
    );
    
    void getNoisyTrLevelDistribution(
        uint32_t& geneId,
        std::vector<uint32_t>& tids,
        std::unordered_map<uint32_t, uint32_t>& exonTrMap,
        std::vector<double>& probVec,
        std::unordered_map<uint32_t, double>& exonLevelProb
    );

    void loadBFH(
        std::string& bfhFile,
        std::string& cellClusterFile,
        Reference& refInfo,
        std::map<std::string, uint32_t>& cellWhiteListMap,
        bool generateNoiseProfile,
        std::unordered_map<std::string, uint32_t>& cellNoisyMap
    ) ;
    void dumpClusterHistoGram(std::string& file_name) ;
    
    inline bool getTrEqExon(uint32_t tid, uint32_t eqId, util::TrRelPos& exTrPos){
        
        if(transcript2EqMap.find(tid) != transcript2EqMap.end()){
            if(transcript2EqMap[tid].find(eqId) != transcript2EqMap[tid].end()){
                exTrPos = transcript2EqMap[tid][eqId] ;
                return true ;
            }else{
                return false ;
            }
        }else{
            return false ;
        }
    }

    inline uint32_t getGeneLevelProbCount(uint32_t geneId, uint32_t count){
        return geneCountHistogram[geneId][count] ;
    }


    std::string exon_info_dir ;
    std::string eqClassCountFile ;
    std::string tr2EqClassFile ;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, util::TrRelPos>> transcript2EqMap ;
    std::vector<uint32_t> eqCountVec ;
    std::vector<double> countProbability ;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> geneCountHistogram ; // Gene id -> (EqClass_Length -> Numebr)
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> noisyGeneCountHistogram ; // Gene id -> (EqClass_Length -> Numebr)
    std::unordered_map<uint32_t,
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>> clusterCountHistogram ; // Cluster id -> (Gene id -> (EqClass Length -> Number))
    std::unordered_map<uint32_t, uint32_t> cell2ClusterMap ;

};

class ExonEqClass2{
    public: 
        ExonEqClass2(
            std::string& fnameIn
        ){
            fname_ = fnameIn ;
            
        }

        struct ExonLevel{
            std::string junction_name ;
            std::vector<uint32_t> trList ;
            size_t labelLength ; 

            ExonLevel(
                std::string& junction_nameIn,
                std::vector<uint32_t>& trListIn,
                size_t labelLengthIn 
            ){
                junction_name = junction_nameIn ;
                trList = trListIn ;
                labelLength = labelLengthIn ;
            }
        };

        void fillGeneLevelMap(Reference& refInfo) ;
        void loadDistributionVector(
            std::string& geneName, 
            std::unordered_map<uint32_t, size_t>& trIdMap,
            std::vector<double> probVec 
        ) ;

        
        //void loadBFH(std::string& bfhFileName) ;


        std::string fname_ ;
        std::unordered_map<std::string, std::vector<ExonLevel> > geneLevelMap ;
        std::vector<double> countProbability ;
};

#endif