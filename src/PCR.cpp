#include "PCR.hpp"

/*
This idea of this PCR model is directly taken from the paper, and mentioned as model 6 in the paper  
although the $\mu$  and $\sigma$ can be given by user to change the model. The implementation is 
is no longer available in the github repository mentioned in the paper. The hope is my implementation
would not daviate too much from the original one. 

Best, Katharine et al. “Computational analysis of stochastic heterogeneity in PCR amplification 
efficiency revealed by single molecule barcoding” Scientific reports vol. 5 14629. 13 Oct. 2015, 
doi:10.1038/srep14629
*/

void PCRClass::runPCRInheritedModel(
    std::map<uint32_t, std::string>& sequenceMap,
    double mu = 0.45,
    double sigma = 0.2
){


    //std::cerr << "Running PCR amplification bias from Katharine et al. \n" ;

    std::vector<double> efficiencyVector(numOfUniqueMolecules) ;
    std::random_device norm{};
    std::mt19937 gen{norm()};
    if (mu > 1.0 || sigma > 1.0){
        mu = 0.45 ;
        sigma = 0.2 ;
    }

    std::normal_distribution<> norm_dist{mu, sigma};

    for(size_t i = 0 ; i < numOfUniqueMolecules ; ++i){
        changeBlockVector[i].setExist() ;
        changeBlockVector[i].setActive() ;
        efficiencyVector[i] = norm_dist(gen) ;

    }

    std::random_device rd; 
    std::random_device rd_right;

    std::mt19937 generator_left(rd()) ;
    std::mt19937 generator_right(rd_right()) ;

    std::uniform_real_distribution<double> distribution(0.0,1.0);
    bool mutateBoth{false} ;

    size_t totalMismatches{0} ;
    size_t cumIndex = numOfUniqueMolecules ;

    for(uint32_t numCycle = 0 ; numCycle < numOfPCRCycles ; ++numCycle){
        // go over the entire vector
        for(uint32_t j = 0 ; j < cumIndex ; ++j){
            uint32_t grandParent_right = (j + cumIndex) % numOfUniqueMolecules ;
            uint32_t grandParent_left = j % numOfUniqueMolecules ;
            // for now we are always active
            if(changeBlockVector[j].exist){
                // check whether to silence left or right child 
                
                double toss_left = distribution(generator_left) ;
                double toss_right = distribution(generator_right) ;

                if (toss_right < efficiencyVector[grandParent_right]){
                    // child will exist 
                    changeBlockVector[j + cumIndex].setExist() ;
                    changeBlockVector[j + cumIndex].setActive() ;
                
                    // will mutate (?)
                    std::pair<uint32_t, uint8_t> mutateInfo ;
                    // NOTE: this step is required just to get seq size
                    size_t fragmentLength = sequenceMap[grandParent_right].size() ;

                    bool isChanged = mutateSeq(
                        mutateInfo,
                        fragmentLength 
                    ) ;

                    if (isChanged) {
                        totalMismatches += 1 ;
                        changeBlockVector[j + cumIndex].editChain.emplace_back(mutateInfo) ;
                    }
                }else{
                    changeBlockVector[j + cumIndex].exist = false ;
                    changeBlockVector[j + cumIndex].active = false ;
                }
                // child won't exist, but we might want to 
                // give this a pass
                if(toss_left < efficiencyVector[grandParent_left]) {
                    // will mutate (?)
                    if(!mutateBoth){
                        continue ;
                    }
                    std::pair<uint32_t, uint8_t> mutateInfo ;
                    // NOTE: this step is required just to get seq size
                    size_t fragmentLength = sequenceMap[grandParent_left].size() ;

                    bool isChanged = mutateSeq(
                        mutateInfo,
                        fragmentLength 
                    ) ;

                    if (isChanged) {
                        totalMismatches += 1 ;
                        changeBlockVector[j].editChain.emplace_back(mutateInfo) ;
                    }
                }else{
                    changeBlockVector[j].active = false ;
                } 
            }
        }
        cumIndex += cumIndex ;
    }
}


void PCRClass::runPCRNormal(
    std::map<uint32_t, std::string>& sequenceMap
){
    for(size_t i = 0 ; i < numOfUniqueMolecules ; ++i){
        changeBlockVector[i].setExist() ;
        changeBlockVector[i].setActive() ;
    }

    std::random_device rd; 
    std::random_device rd_right;

    std::mt19937 generator_left(rd()) ;
    std::mt19937 generator_right(rd_right()) ;

    std::uniform_real_distribution<double> distribution(0.0,1.0);
    bool mutateBoth{false} ;

    size_t totalMismatches{0} ;
    size_t cumIndex = numOfUniqueMolecules ;

    for(uint32_t numCycle = 0 ; numCycle < numOfPCRCycles ; ++numCycle){
        // go over the entire vector
        for(uint32_t j = 0 ; j < cumIndex ; ++j){
            // for now we are always active
            if(changeBlockVector[j].active){
                // check whether to silence left or right child 
                
                double toss_left = distribution(generator_left) ;
                double toss_right = distribution(generator_right) ;

                if (toss_right < captureProbability){
                    // child will exist 
                    changeBlockVector[j + cumIndex].setExist() ;
                    changeBlockVector[j + cumIndex].setActive() ;
                
                    // will mutate (?)
                    std::pair<uint32_t, uint8_t> mutateInfo ;
                    // NOTE: this step is required just to get seq size
                    uint32_t grandParent = (j + cumIndex) % numOfUniqueMolecules ;
                    size_t fragmentLength = sequenceMap[grandParent].size() ;

                    bool isChanged = mutateSeq(
                        mutateInfo,
                        fragmentLength 
                    ) ;

                    if (isChanged) {
                        totalMismatches += 1 ;
                        changeBlockVector[j + cumIndex].editChain.emplace_back(mutateInfo) ;
                    }
                }else{
                    changeBlockVector[j + cumIndex].exist = false ;
                    changeBlockVector[j + cumIndex].active = false ;
                }
                // child won't exist, but we might want to 
                // give this a pass
                if(toss_left < captureProbability) {
                    // will mutate (?)
                    if(!mutateBoth){
                        continue ;
                    }
                    std::pair<uint32_t, uint8_t> mutateInfo ;
                    // NOTE: this step is required just to get seq size
                    uint32_t grandParent = j % numOfUniqueMolecules ;
                    size_t fragmentLength = sequenceMap[grandParent].size() ;

                    bool isChanged = mutateSeq(
                        mutateInfo,
                        fragmentLength 
                    ) ;

                    if (isChanged) {
                        totalMismatches += 1 ;
                        changeBlockVector[j].editChain.emplace_back(mutateInfo) ;
                    }
                }else{
                    changeBlockVector[j].active = false ;
                } 
            }
            else{
                // TODO: introduce efficiency 
            }
        }
        cumIndex += cumIndex ;
    }
}

std::string PCRClass::getEditedSequence(
    uint32_t& index,
    std::string& parentSeq 
){
    // start from editChain and make the changes 
    if(changeBlockVector[index].editChain.size() > 0){
        auto& editChain = changeBlockVector[index].editChain ; 
        for(auto ep : editChain){
            parentSeq[ep.first] = charMap[parentSeq[ep.first]][ep.second] ;
        }
    }
    return parentSeq ;
}

bool PCRClass::constructOrCacheSequence(
    uint32_t ind,
    std::map<uint32_t, std::string>& sequenceMap
){

    bool validBlock{true} ;
    
    //check if this molecule exists or not
    if(!changeBlockVector[ind].exist){
        validBlock = false ;
        return validBlock ;
    }

    size_t blockSize{numOfUniqueMolecules} ;
    auto grandParentId = ind % blockSize ;

    uint32_t immediateParentId = util::calculateParent(ind, blockSize) ;
    auto it = sequenceMap.find(immediateParentId) ;
    if(it != sequenceMap.end()){
        // this baby is cached
        auto parentSeq = it->second ;
        sequenceMap[ind] = getEditedSequence(ind, parentSeq) ;
    }else{
        // cache this baby
        // search the tree (recusrsively)
        std::vector<uint32_t> nodeStack ;
        
        nodeStack.push_back(ind) ;
        nodeStack.push_back(immediateParentId) ;


        // push in stack 
        while(true){
            immediateParentId = util::calculateParent(immediateParentId, blockSize) ;
            if(immediateParentId > ind){
                std::cerr  << "Should not happen !!!\n" ;
                std::exit(1) ;
            }

            it = sequenceMap.find(immediateParentId) ;
            if(it != sequenceMap.end() || immediateParentId < blockSize){
                break ;
            }else{
                nodeStack.push_back(immediateParentId) ;
            }
        }

        // get parent seq
        std::string parentSeq ;
        if (it != sequenceMap.end()){
            parentSeq = it->second ;
        }else{
            std::cerr << "Should not happen\n" ;
            std::exit(1) ;
        }        
        // pop from stack
        for(int j = nodeStack.size() - 1 ; j >= 0; --j){
            sequenceMap[nodeStack[j]] = getEditedSequence(nodeStack[j], parentSeq) ; 
        }
    }
    
    return validBlock ;

}

