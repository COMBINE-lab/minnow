#ifndef GFA_READER_HPP
#define GFA_READER_HPP

#include "MinnowUtil.hpp"
#include "ReferenceInfo.hpp"


class GFAReader{
    public:

    GFAReader(
        std::string gfaFileIn
    ){
        gfaFileName_ = gfaFileIn ;
    }

    void parseFile(
        Reference& refInfo
    ) ;

    void readUnitigs() ;

    void updateEqClass(
      std::string& transcriptName, 
      std::vector<std::pair<size_t, bool>>& contigVec,
      Reference& refInfo,
      size_t overlap
    );
    
    struct posInConfig{
        uint32_t sposInContig ;
        uint32_t eposInContig ;
        bool ore ;

        posInConfig(){} ;

        posInConfig(
            uint32_t sposInContigIn ,
            uint32_t eposInContigIn ,
            bool oreIn 
        ):
        sposInContig(sposInContigIn),
        eposInContig(eposInContigIn),
        ore(oreIn)
        {}
    }; 

    using trContigMap = std::unordered_map<uint32_t, std::vector<posInConfig>> ;
    std::unordered_map<size_t, trContigMap> eqClassMap ; // segment -> {}
    std::unordered_map<size_t, uint32_t> eqClassMapCounts ; 

    std::unordered_map<uint32_t, std::vector<size_t>> trSegmentMap ; // tr -> segmenteqclass
    
    std::string gfaFileName_ ;
    std::unique_ptr<std::ifstream> file ; 
    std::unordered_map<size_t, std::string> unitigMap ;



}; 

#endif 
