#ifndef PCR_HEADER
#define PCR_HEADER

#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <cmath>
#include <type_traits>
#include <memory>
#include <mutex>
#include <thread>
#include <algorithm> 
#include <sstream>
#include <map>

#include "MinnowUtil.hpp"
#include "macros.hpp"

using edit_pair = std::pair<uint32_t, uint8_t> ;


class PCRClass{
public:
        
    PCRClass(
        size_t numOfUniqueMoleculeIn, 
        size_t numOfPCRCyclesIn,
        double mutationProbabilityIn,
        double captureProbabilityIn,
        double errorProbabilityIn,
        bool switchOnEffModelIn = false,
        double muIn = 0.45,
        double sigmaIn = 0.2 
    ): 
    numOfUniqueMolecules(numOfUniqueMoleculeIn),
    numOfPCRCycles(numOfPCRCyclesIn),
    mutationProbability(mutationProbabilityIn),
    captureProbability(captureProbabilityIn),
    errorProbability(errorProbabilityIn),
    switchOnEffModel(switchOnEffModelIn),
    mu(muIn),
    sigma(sigmaIn)
    {
        //std::cerr << "PCR :::: errorProbability " << errorProbability << "\n";
        totNum =  std::pow(2, numOfPCRCycles) * numOfUniqueMolecules ;
        changeBlockVector.resize(totNum) ;
    }   

    struct changeBlock{
        std::vector<edit_pair> editChain ; // vector if changed position and char
        bool active{false} ; // this variable is important while creating, not while decoding
        bool exist{false} ; // if it really exists (effective for a dropped out molecule)
        bool sameAsParent{false} ;

        changeBlock(){} // default
        void setActive(){active = true ;} // captured in this cycle 
        void setExist(){exist=true ;} // this block exists 
        void setSameAsParent(){sameAsParent=true ;} // this block exists 

    };

    void runPCR(
       std::map<uint32_t, std::string>& sequenceMap
    ){
        if(!switchOnEffModel){
            //std::cerr << "PCR:::: Going to call normal \n" ;
            runPCRNormal(
                sequenceMap
            ) ;
        }else{
            //std::cerr << "PCR:::: Going to call PCR model6  \n " ;
            runPCRInheritedModel(
                sequenceMap,
                mu,
                sigma 
            ) ;
        }
    }

    void runPCRNormal(
        std::map<uint32_t, std::string>& sequenceMap
    ) ;

    void runPCRInheritedModel(
        std::map<uint32_t, std::string>& sequenceMap,
        double mu,
        double sigma
    ) ;

    std::string getEditedSequence(
        uint32_t& index,
        std::string& parentSeq 
    ) ;

    bool constructOrCacheSequence(
        uint32_t ind,
        std::map<uint32_t, std::string>& sequenceMap
    ) ;

    inline void checkChainLength(){
        size_t numOfNontrivalChains{0} ;
        for(size_t i = 0 ; i < changeBlockVector.size() ; ++i){
            if(changeBlockVector[i].editChain.size() > 1){
                numOfNontrivalChains++ ;
            }
        }
        //std::cerr << "Number of more than one element chains " << numOfNontrivalChains << "\n" ;
    };

    inline bool mutateSeq(
        edit_pair& pc, 
        size_t length
    ){
        // Seed with a real random value, if available
        std::random_device r;
        std::mt19937 generator(r()) ;
        // to change ?
        std::uniform_real_distribution<double> distribution(0.0,1.0);

        // if there will be an error or not 
        double toss = distribution(generator) ;
        bool isChanged{false} ;
        if (toss < errorProbability){
            // which nucliotide to change 
            std::uniform_int_distribution<size_t> index_distribution(0,length-1);
            // what to change to 
            std::uniform_int_distribution<uint8_t> snp_distribution(0, 2);

            uint32_t randind = static_cast<uint32_t>(index_distribution(generator)) ;
            isChanged = true ;
            pc = std::make_pair(randind, static_cast<uint8_t>(snp_distribution(generator))) ;
        }    
        return isChanged ;
    };

    std::map<char , std::vector<char>> charMap = {
        {'A', {'T', 'G', 'C'}}, 
        {'T', {'A', 'G', 'C'}},
        {'G', {'A', 'T', 'C'}},
        {'C', {'A', 'G', 'T'}} 
    };

    size_t numOfUniqueMolecules ; 
    size_t numOfPCRCycles ;
    size_t totNum ;
    double mutationProbability ;
    double captureProbability ;
    double errorProbability ;
    bool switchOnEffModel ;
    double mu;
    double sigma ;
    std::vector<changeBlock> changeBlockVector ;

};

#endif // PCR_HEADER