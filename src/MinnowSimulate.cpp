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
#include <functional>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "ProgOpts.hpp"
#include "MinnowUtil.hpp"
#include "MatrixParser.hpp"
#include "Transcript.hpp"
#include "ReferenceInfo.hpp"
#include "xxhash.h"
#include "FASTAParser.hpp"
#include "SpinLock.hpp"
#include "ScopedTimer.hpp"
#include "MinnowFS.hpp"
#include "string_view.hpp"
#include "concurrentqueue.h"
#include "PCR.hpp"
#include "GFAReader.hpp"

#include "zstr.hpp"

#if defined __APPLE__
using SpinLockT = SpinLock;
#else
using SpinLockT = std::mutex;
#endif

#define MAX_FRAGLENGTH 1000
#define MAX_QUEUE_SIZE 20
#define READ_LEN 100

#define FRAGMENT_END_DIST 404 + READ_LEN // this is from the empirical e^6
#define FRAGMENT_START_DIST 53 + READ_LEN // this is from the empirical limit e^4
#define FRAGMENT_RANGE (FRAGMENT_END_DIST - FRAGMENT_START_DIST)



#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)


// Convert the the block size to minimal to limit the 
// memory usage 
struct MyTraits : public moodycamel::ConcurrentQueueDefaultTraits
{
	static const size_t BLOCK_SIZE = 2; // Use smaller blocks (32 default)
};


// Convention we assume closed interval in start and end [start, end] 
inline uint32_t normalizedStartPos(uint32_t fragmentStart, uint32_t fragmentEnd, uint32_t startPosOld){
    int newRange = static_cast<int>(fragmentEnd - fragmentStart) ; 
    return static_cast<uint32_t>( (((startPosOld - FRAGMENT_START_DIST ) * newRange) / FRAGMENT_RANGE ) + fragmentStart );
}

std::map<char , std::vector<char>> charMap = {
    {'A', {'T', 'G', 'C', 'N'}}, 
    {'T', {'A', 'G', 'C', 'N'}},
    {'G', {'A', 'T', 'C', 'N'}},
    {'C', {'A', 'G', 'T', 'N'}} 
};

std::map<char, size_t> nuclMap = {
    {'A', 0},
    {'T', 1},
    {'G', 2},
    {'C', 3},
    {'N', 4}
} ;

std::vector<char> nuclMapRev = {'A','T','G','C','N'} ;

std::map<char , std::vector<char>> charMapUMI = {
    {'A', {'T', 'G', 'C'}}, 
    {'T', {'A', 'G', 'C'}},
    {'G', {'A', 'T', 'C'}},
    {'C', {'A', 'G', 'T'}} 
};

using string_pair = std::pair<std::string, std::string> ;
using stream_pair_ptr = std::pair<std::stringstream*, std::stringstream*> ;

using ErrorMatrix = std::vector<std::vector<std::vector<double>>> ;

void indexDistribution(
    int UMIPoolSize,
    int totalNumberOfMolecules,
    std::vector<int>& indexList
){
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, UMIPoolSize-1);

    indexList.resize(totalNumberOfMolecules) ;
    for(int i = 0; i < totalNumberOfMolecules; ++i){
        indexList[i] = dist(rd) ;
    }
}


int binomialCoeff(int n, int k){ 
    int res = 1; 
    if ( k > n - k ) 
        k = n - k; 
    for (int i = 0; i < k; ++i) { 
        res *= (n - i); 
        res /= (i + 1); 
    } 
    return res; 
} 


void binomialExpansion(
    int n,
    double errorProbability,
    std::vector<double>& probVector
){
    int numOfElements = std::pow(2, n) ;
    int expandedNumberOfElements = n+1 ;

    int innerIndex{0} ;    
    probVector.resize(numOfElements) ;
    for(int i = 0; i < expandedNumberOfElements; ++i){
        int coeff = binomialCoeff(n, i);
        for(int j = 0 ; j < coeff; ++j){
            probVector[innerIndex++] = std::pow(errorProbability, i);
        }
    }
}



std::string imputeErrorInString(
    double errorProbability,
    std::string& seq
){
    
    //std::vector<double> prob(seq.size()) ;
    std::random_device r;
    std::mt19937 generator(r()) ;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    //auto generator = std::bind(distribution, engine) ;
    //std::generate_n(prob.begin(), seq.size(), generator) ;

    //bool isChanged{false} ;

    for(size_t i = 0; i < seq.size() ; ++i ){
        if(distribution(generator) < errorProbability){
            seq[i] = charMap[seq[i]][rand()%4] ;
        }
    }
    return seq ;
}

std::string imputeIlluminaModel(
    std::vector<std::vector<std::vector<double>>>& errorModel,
    std::string& seq 
){
    std::random_device rdg;
    std::mt19937 geng(rdg());

    for(size_t pos = 0 ; pos < seq.size() ; ++pos){
       int ind = nuclMap[seq[pos]] ;
       auto probVec = errorModel[pos][ind]  ;
	   
       std::discrete_distribution<> dg(probVec.begin(), probVec.end()) ;
       auto mutatedInd = dg(geng) ;
       if (mutatedInd != ind){

           seq[pos] = nuclMapRev[mutatedInd] ;
       }

    }
    return seq ;
}


std::string imputeErrorInStringUMI(
    double errorProbability,
    std::string& seq
){
    
    //std::vector<double> prob(seq.size()) ;
    std::random_device r;
    std::mt19937 generator(r()) ;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    //auto generator = std::bind(distribution, engine) ;
    //std::generate_n(prob.begin(), seq.size(), generator) ;

    //bool isChanged{false} ;

    for(size_t i = 0; i < seq.size() ; ++i ){
        if(distribution(generator) < errorProbability){
            seq[i] = charMapUMI[seq[i]][rand()%3] ;
        }
    }
    return seq ;
}

//Generate from Segment
uint32_t generateFragmentSegment(
    std::string& fragSeq,
    std::string& readSeq,
    bool ore 
){
    auto refLen = fragSeq.size() ;
    uint32_t startPos{0} ;
    //uint32_t endPos{0} ;
    if(refLen >= FRAGMENT_END_DIST){
        auto distanceFromEnd = util::generateStartPosition() ;
        if(distanceFromEnd < READ_LEN){
            distanceFromEnd += READ_LEN ;
            startPos = refLen - distanceFromEnd ;
        }
        //endPos = startPos + READ_LEN ;        
    }else{
        startPos = 0 ;
        //endPos = READ_LEN ;
    }

    if(startPos > refLen - READ_LEN){
        std::cerr << "RIASE RED FLAG\n" ;
    }

    if(!ore){
        //std::cerr << "hi\n" ;
        auto revStr = util::revcomp(fragSeq) ;
        readSeq = revStr.substr(startPos, READ_LEN) ;
    }else{
        readSeq = fragSeq.substr(startPos, READ_LEN) ;
    }

    return startPos ;
}


uint32_t generateFragmentSeq(
    Transcript& tr,
    std::string& fragSeq 
){
    
    auto& refLen = tr.RefLength ;

    if (refLen < READ_LEN){
        std::cerr << "Transcript is too short\n" ;
        fragSeq = std::string{""} ; 
    }
    if (refLen == READ_LEN){
        fragSeq = std::string(
            tr.Sequence(),
            READ_LEN
        );
        return 0 ;
    }

    uint32_t fragmentLength = MAX_FRAGLENGTH ;
    if (refLen < fragmentLength){
        fragmentLength = refLen ;
    }
    uint32_t fragmentStartPos = refLen - fragmentLength ;
    uint32_t fragmentEndPos = refLen - READ_LEN - 1 ;
    // Generate a random number from interval (4,6)
    auto distanceFromEnd = util::generateStartPosition() + READ_LEN;
    //auto fragmentRange = fragmentEndPos - fragmentStartPos ;
    int numPads{0} ;
    uint32_t startPos{0} ;
    uint32_t toSample{READ_LEN} ;

    if(distanceFromEnd > refLen){
        // Scale the distance to be in the range of fragment
        // Length
        startPos = normalizedStartPos(fragmentStartPos, fragmentEndPos, distanceFromEnd) ;
    }else{
        startPos = refLen - distanceFromEnd ;
        if(distanceFromEnd < READ_LEN){
            // We need to pad Ns
            numPads = READ_LEN - distanceFromEnd ;
            toSample = distanceFromEnd ;
        }
    }

    fragSeq = std::string(
            tr.Sequence() + startPos,
            toSample
    );
    // NOTE: rethink the logic 
    if(numPads > READ_LEN/2){
        std::cout << "fragmentLength " << fragmentLength << "\t"
                  << "RefLen " << refLen << "\t" 
                  << "startPos " << startPos << "\t" 
                  << "READ_LEN " << READ_LEN << "\n" ; 
         
    }
    
    fragSeq = fragSeq +  std::string( numPads, 'N') ; 
    return startPos ;
}

uint32_t truncatedNormalSampling(
    uint32_t mu ,
    uint32_t sigma,
    uint32_t max_length 
){
    std::random_device rd{};
    std::mt19937 gen{rd()};
    //uint32_t fragLen{0} ;

    size_t iterations{0} ;
    std::normal_distribution<> d{ static_cast<double>(mu), static_cast<double>(sigma)} ;
    while(true){
        auto fragLen = static_cast<uint32_t>(std::round(d(gen))); 
        if(fragLen >= READ_LEN && fragLen <= max_length){
            return fragLen ;
        }
        iterations++ ;
        if(iterations > 10000){
            std::cerr << "Try to set a smaller \\sigma, this random generation crosses "<< iterations 
            << " returning READ_LEN  \n"  ;
            return READ_LEN ;
        }
    }
}


uint32_t generateSegmentRSPD(
    std::string& trSeq,
    std::string& readSeq,
    std::vector<uint32_t>& rspdVec
){
    auto refLen = trSeq.size() ;


    if (refLen < READ_LEN){
        std::cerr << "Transcript is too short\n" ;
        std::exit(1) ;
    }

    if (refLen == READ_LEN){
        readSeq = std::string(trSeq) ;
        return 0 ;
    }
    // define a probVec, 
    // where first element is 
    // probVec[0] = prob(frag_length = READ_LEN) , 
    // probVec[1] = prob(frag_length = READ_LEN+1) ,
    // ...
    // probVec[segmentStartPos] = prob(frag_length = segmentStartPos) ,
    std::random_device rd;
    std::mt19937 gen(rd());

    std::discrete_distribution<> d(rspdVec.begin()+READ_LEN, rspdVec.begin()+(refLen-1)) ;

    uint32_t startPos = d(gen) ;
    if(startPos + READ_LEN > refLen){
        std::cerr<< "we are on the edge: d.size() " << refLen-READ_LEN << "\n" ; 
        std::cerr << startPos << "\t" << refLen << "\n" ;
        std::exit(1) ;
    }
    readSeq = trSeq.substr(startPos, READ_LEN) ;

    return startPos ;
}

uint32_t generateFragmentSeqFromSegment(
    std::string& trSeq,
    uint32_t mu,
    std::string& readSeq
){
    
    auto refLen = trSeq.size() ;


    if (refLen < READ_LEN){
        std::cerr << "Transcript is too short\n" ;
        std::exit(1) ;
    }

    if (refLen == READ_LEN){
        readSeq = std::string(trSeq) ;
        return 0 ;
    }
    
    uint32_t startPos{0} ;
    if(mu < READ_LEN){
        startPos = util::generateRandomNumber<uint32_t>(0, refLen - READ_LEN) ;
        readSeq = trSeq.substr(startPos, READ_LEN) ;
        if(readSeq == ""){
            std::cerr << "low mu : read seq empty ::: " << startPos << " \t" << refLen << "\n" ;
            std::exit(1) ;
        }
        return startPos ;
    }

    auto fragLen = truncatedNormalSampling(mu, 10, refLen) ;
    startPos = refLen - fragLen ;
    readSeq = trSeq.substr(startPos, READ_LEN) ;
        if(readSeq == ""){
            std::cerr << "high mu: read seq empty ::: " << startPos << " \t" << refLen << "\n" ;
            std::exit(1) ;
        }

    return startPos ;
}

// When equivaence classes are not used then 
// we samnple from an exponential distribution to 
// match that of the 10X dataset. 

uint32_t generateFragmentSeq(
    std::string& trSeq,
    std::string& fragSeq 
){
    
    auto refLen = trSeq.size() ;

    if (refLen < READ_LEN){
        std::cerr << "Transcript is too short\n" ;
        fragSeq = std::string{""} ; 
    }
    if (refLen == READ_LEN){
        fragSeq = std::string(trSeq) ;
        return 0 ;
    }

    uint32_t fragmentLength = MAX_FRAGLENGTH ;
    if (refLen < fragmentLength){
        fragmentLength = refLen ;
    }
    uint32_t fragmentStartPos = refLen - fragmentLength ;
    uint32_t fragmentEndPos = refLen - READ_LEN - 1 ;
    // Generate a random number from interval (4,6)
    auto distanceFromEnd = util::generateStartPosition() + READ_LEN;
    //auto fragmentRange = fragmentEndPos - fragmentStartPos ;
    int numPads{0} ;
    uint32_t startPos{0} ;
    uint32_t toSample{READ_LEN} ;

    if(distanceFromEnd > refLen){
        // Scale the distance to be in the range of fragment
        // Length
        startPos = normalizedStartPos(fragmentStartPos, fragmentEndPos, distanceFromEnd) ;
    }else{
        startPos = refLen - distanceFromEnd ;
        if(distanceFromEnd < READ_LEN){
            // We need to pad Ns
            numPads = READ_LEN - distanceFromEnd ;
            toSample = distanceFromEnd ;
        }
    }

    

    fragSeq = trSeq.substr(startPos, toSample) ;
    // NOTE: rethink the logic 
    if(numPads > READ_LEN/2){
        std::cout << "fragmentLength " << fragmentLength << "\t"
                  << "RefLen " << refLen << "\t" 
                  << "startPos " << startPos << "\t" 
                  << "READ_LEN " << READ_LEN << "\n" ; 
         
    }
    
    fragSeq = fragSeq +  std::string( numPads, 'N') ; 
    return startPos ;
}

uint32_t generateFragmentSeq(
    Transcript& tr,
    std::string& fragSeq,
    uint32_t& fragmentStart,
    uint32_t& fragmentEnd
){
    uint32_t toSample{READ_LEN} ;
    auto& refLen = tr.RefLength ;


    auto startPos = util::generateRandomNumber<uint32_t>(fragmentStart, (fragmentEnd + 1) - READ_LEN) ;
    int numPads{0} ;

    auto distanceFromEnd = refLen - startPos ;
    if (distanceFromEnd < READ_LEN){
        toSample = distanceFromEnd ;
        numPads = READ_LEN - distanceFromEnd ;
    }
    
    fragSeq = std::string(
        tr.Sequence() + startPos, 
        toSample
    ) ;
    fragSeq = fragSeq + std::string(numPads, 'N') ;
    return startPos ;

}

void clipSequence(std::string& parentSeq, uint32_t fStart, uint32_t fEnd, std::string& readSeq, uint32_t& startPos){
    std::string fragSeqPart = parentSeq.substr(CB_LENGTH + UMI_LENGTH) ;
    size_t absReadStartPos{0} ;
    size_t readStartPos{0} ;
    if (fragSeqPart.size() >= READ_LEN){
        absReadStartPos = util::generateRandomNumber<uint32_t>(fStart, (fEnd+1) - READ_LEN) ;
        readStartPos = absReadStartPos - fStart ;
        //absReadStartPos = util::generateRandomNumber<uint32_t>(fStart, (fEnd+1) - READ_LEN) ;
        //readStartPos = fStart ;                        
    }else{
        std::cerr << "DEBUG: ==== Check what's going on\n" ;
    }

    readSeq = fragSeqPart.substr(readStartPos, READ_LEN) ;
    startPos = readStartPos ;
}


//NOTE: deprecated
template <typename T>
void printVector(std::vector<T>& vec){
    for(auto& v: vec){
        std::cout << v << "\t" ;
    }
    std::cout << "\n" ;
}


void doPCRBulkDBG(
    std::string& cellName, 
    uint32_t dupCount,
    std::vector<util::CellBarcodeUMISegmentBlock>& uniqueMolecules,
    std::vector<uint32_t> rspdVec,
    std::vector<Transcript>& transcripts,
    GFAReader* dbgPtr,
    uint32_t& numOfPCRCycles,
    double& errorRate,
    bool switchOnEffModel,
    ErrorMatrix& errorModel,
    fmt::MemoryWriter& sstream_left,
    fmt::MemoryWriter& sstream_right
){
    PCRClass* pcrClassPtr{nullptr} ;
    bool rspdMode{false} ;

    if(rspdVec.size() > 0){
        rspdMode = true ;
    }

    //fragmentLength
    //size_t fragmentLength{MAX_FRAGLENGTH} ;
    double efficiencyProb{0.98} ;
    double mutationProb{0.01} ;
    double seqeneceErrorProb{0.005} ;

    auto totNum = std::pow(2, numOfPCRCycles) * uniqueMolecules.size() ;

    // Create a map of sequences
    // create parent sequences
    std::map<uint32_t, std::string> sequenceMap ;
    std::map<uint32_t, std::string> readSeqMap ;
    std::map<uint32_t, uint32_t> startPosMap ;

    std::map<uint32_t, uint32_t> absStartPosMap ;



    std::unordered_map<uint32_t, std::pair<size_t, size_t>> startPositionMap ;
    
    std::vector<uint32_t> muMap(uniqueMolecules.size()) ;
    size_t uniqReadId{0} ;

    auto unitigMap = dbgPtr->unitigMap ;  
    size_t smallerMu{0} ; 
    size_t biggerMu{0} ; 


    for(uint32_t i = 0 ; i < uniqueMolecules.size() ; ++i){
        
        auto& tid = uniqueMolecules[i].transcriptId ;
        //auto& sid = uniqueMolecules[i].segmentId ;

        auto& tr = transcripts[tid] ;

        auto segStartPos = uniqueMolecules[i].segmentStart ;
        auto segEndPos = uniqueMolecules[i].segmentEnd ;
        
        //auto oldSeLen = segEndPos - segStartPos + 1 ;

        // Sampling strategy 
        // Take a 0 truncated normal distribution where
        // mu = segStartPos + segmentLength - 1
        // sigma = 10 (take as input) 


        std::string readSeq ;
        uint32_t startPos{0} ;
        uint32_t absStartPos{0} ;
        std::string fragmentSeq{""} ;
        uint32_t fragmentStartPos{0} ;
        //uint32_t fragmentEnd{0} ;
        uint32_t fragLen = 0;
        
        if(tr.RefLength < MAX_FRAGLENGTH) {
                fragmentStartPos = 0 ;
               fragLen = tr.RefLength;
        }else {
            fragmentStartPos = tr.RefLength - MAX_FRAGLENGTH ;
            fragLen = MAX_FRAGLENGTH ;
        }
        // clip the segment start pos 
        auto segStartPosOld = segStartPos ;
        bool max_frag_on{false} ;

        if(segStartPos < fragmentStartPos){
            segStartPos = fragmentStartPos ;
            max_frag_on = true ;
        }
        
        uint32_t slack_val{50} ; // mohsen number

        if(rspdMode){
            //fragmentStartPos = std::min(segStartPos, static_cast<uint32_t>(MAX_FRAGLENGTH)) ;
            //fragLen = tr.RefLength - segStartPos ;

            auto segmentLength = segEndPos - segStartPos + 1 ;
            if(!max_frag_on){
                if (segStartPos > slack_val){
                    segStartPos = segStartPos - slack_val ;
                }
            }

            if(segmentLength < READ_LEN){
                fragLen = READ_LEN + slack_val ;
            }else{
                fragLen = segmentLength ;
            }

            if(segStartPos + fragLen > tr.RefLength){
                fragLen = tr.RefLength - segStartPos + 1; 
            }

            if(segStartPos + fragLen > tr.RefLength){
                std::cerr << fragLen << "\t" << tr.RefLength << "\t" << segStartPosOld << "\t" << segStartPos << "\t" << segEndPos << "\n" ;
            }

            //fragmentEnd = segStartPos + fragLen - 1 ;
            fragmentSeq = std::string{tr.Sequence() + segStartPos, fragLen} ;
            //store the abs start
            absStartPos = segStartPos ;

            startPos = generateSegmentRSPD(fragmentSeq, readSeq, rspdVec) ;
        }else{
            auto segmentLength = segEndPos - segStartPos + 1 ;
            // scale segment start pos
            segStartPos = segStartPos - fragmentStartPos ; 
            auto mu = static_cast<uint32_t>(fragLen - (segStartPos + std::round(segmentLength/2))) ;
            muMap[i] = mu ;
            // Sanity check
            if(mu < READ_LEN){
                smallerMu++ ;
            }else{
                biggerMu++ ;
            }
            //fragmentEnd = fragmentStartPos + fragLen - 1 ;
            // goes through PCR
            fragmentSeq = std::string{tr.Sequence() + fragmentStartPos, fragLen} ;
            startPos = generateFragmentSeqFromSegment(fragmentSeq, mu, readSeq) ;
            absStartPos = fragmentStartPos ;
        }

        

        sequenceMap[i] = uniqueMolecules[i].cellBarCode +
                         uniqueMolecules[i].UMICode + 
                         fragmentSeq ;
        
        readSeqMap[i] = readSeq ;
        startPosMap[i] = startPos ;
        absStartPosMap[i] = absStartPos ;
    }
    //std::cerr << "smaller mu " << smallerMu << "\t" << "bigger mu " << biggerMu << "\n" ;

    std::vector<util::EditBlock> editVector;
    editVector.resize(totNum) ;
    for(uint32_t i = 0 ; i < uniqueMolecules.size(); ++i){
        editVector[i] = {-1, -1, true} ;
    }

    bool doPCR{true} ;
    if(doPCR){
        pcrClassPtr = new PCRClass(
            uniqueMolecules.size(),
            numOfPCRCycles,
            mutationProb,
            efficiencyProb,
            errorRate,
            switchOnEffModel
        ) ;
        
        pcrClassPtr->runPCR(
            sequenceMap
        );
    }

    std::vector<uint32_t> indexArray ;
    indexArray.resize(totNum) ;
    for(uint32_t i = 0; i < totNum ; ++i){
        indexArray[i] = i ;
    }


    size_t numOfWritten{0} ;

    if(dupCount == 0){
      dupCount = MAXNUM ;
    }
    if (dupCount > totNum){
      dupCount = totNum ;
    }

    auto numMol = std::min(indexArray.size(), static_cast<size_t>(dupCount)) ;
    // Write down the normal stuff
    {
      for(uint32_t i = 0; i < uniqueMolecules.size(); ++i){
        uint32_t ind = i ;
                 numOfWritten++ ;

                 std::string modifiedCellName = sequenceMap[ind].substr(0, CB_LENGTH) ;
                //bool ore = uniqueMolecules[ind].count ; 

                //auto fragmentStr = sequenceMap[ind].substr(CB_LENGTH + UMI_LENGTH) ; 
                std::string readSeq = readSeqMap[ind] ;
                //auto startPos = generateFragmentSegment(fragmentStr, readSeq, ore) ;
                auto startPos = startPosMap[ind] ;
                auto absStartPos = absStartPosMap[ind] + startPos ;


                sstream_left << "@" << modifiedCellName 
                             << ":" << transcripts[uniqueMolecules[ind].transcriptId].RefName
                             << ":" << absStartPos
                             << ":" << ind
                             << ":" << uniqReadId
                             << "\n" ;

                std::string tagSeq{sequenceMap[ind].substr(0, CB_LENGTH + UMI_LENGTH)} ;
                sstream_left << imputeErrorInStringUMI(seqeneceErrorProb, tagSeq) << "\n" ;
                sstream_left << "+\n" ;
                sstream_left << std::string(CB_LENGTH + UMI_LENGTH, 'N') << "\n" ;                 


                sstream_right << "@" << modifiedCellName
                              << ":" << transcripts[uniqueMolecules[ind].transcriptId].RefName
                              << ":" << absStartPos
                  //<< ":" << (ore ? "+" : "-")
                              << ":" << ind
                              << ":" << uniqReadId 
                              << "\n" ;

                if(errorModel.size() != 0){
                    sstream_right << imputeIlluminaModel(errorModel, readSeq) << "\n" ;
                }else{
                    sstream_right << imputeErrorInString(seqeneceErrorProb, readSeq) << "\n" ;
                }
                sstream_right << "+\n" ;
                sstream_right << std::string(READ_LEN, 'N') << "\n" ; 
                uniqReadId++ ;              
            

                if(numOfWritten >= numMol)
                    break ;
      }
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(indexArray.begin(), indexArray.end(), g) ;

    //std::sort(indexArray.begin(), indexArray.begin() + numMol) ;
    size_t blockSize = uniqueMolecules.size() ;

    if(doPCR && numOfWritten < numMol){
        for(uint32_t i = 0 ; i < indexArray.size() ; ++i){
            //std::cerr << "iteration " << i << "\n" ;
            auto ind = indexArray[i] ;

            if(ind < blockSize){
            
                continue ;
            }

            bool validBlock = pcrClassPtr->constructOrCacheSequence(ind, sequenceMap) ;
            
            auto grandParentId = ind % blockSize ; 

            if(validBlock){

                auto it = sequenceMap.find(ind) ;
                if(it != sequenceMap.end()){
                    auto parentSeq = it->second ;

                    std::string modifiedCellName = sequenceMap[grandParentId].substr(0, CB_LENGTH) ;
                    auto fragmentStr = sequenceMap[grandParentId].substr(CB_LENGTH + UMI_LENGTH) ; 
                    auto mu = muMap[grandParentId] ;
                    uint32_t startPos{0} ;
                    
                    std::string readSeq ;
                    //auto startPos = generateFragmentSegment(fragmentStr, readSeq, ore) ;
                    if(rspdMode){
                        startPos = generateSegmentRSPD(fragmentStr, readSeq, rspdVec) ;
                    }else{
                        startPos = generateFragmentSeqFromSegment(fragmentStr, mu, readSeq) ;
                    }
                    
                    uint32_t absStartPos = absStartPosMap[grandParentId] + startPos ;

                    sstream_left << "@" << modifiedCellName
                             << ":" << transcripts[uniqueMolecules[grandParentId].transcriptId].RefName
                             << ":" << absStartPos 
                             << ":" << ind
                             << ":" << uniqReadId 
                             << "\n" ;

                    std::string tagSeq{parentSeq.substr(0, CB_LENGTH + UMI_LENGTH)} ;
                    sstream_left << imputeErrorInStringUMI(seqeneceErrorProb, tagSeq) << "\n" ;
                    sstream_left << "+\n" ;
                    sstream_left << std::string(CB_LENGTH + UMI_LENGTH, 'N') << "\n" ; 

                    //bool ore = uniqueMolecules[grandParentId].count ; 
                

                    sstream_right << "@" << modifiedCellName
                                << ":" << transcripts[uniqueMolecules[grandParentId].transcriptId].RefName
                                << ":" << absStartPos
                      //            << ":" << (ore ? "+" : "-")
                                << ":" << ind
                                << ":" << uniqReadId 
                                << "\n" ;


                    if(errorModel.size() != 0){
                        sstream_right << imputeIlluminaModel(errorModel, readSeq) << "\n" ;
                    }else{
                        sstream_right << imputeErrorInString(seqeneceErrorProb, readSeq) << "\n" ;
                    }
                    sstream_right << "+\n" ;
                    sstream_right << std::string(READ_LEN, 'N') << "\n" ; 
                    uniqReadId++ ;    
                    
                        
                }else{
                    std::cerr << "Should not happen\n" ;
                    std::exit(1) ;
                }

                numOfWritten++ ;
            }
        
            if(numOfWritten >= numMol)
                break ;
        }
        //pcrClassPtr->checkChainLength() ;
        if(numOfWritten > numMol){
            std::cerr << "WTF WE HAVE MORE THAN 100K mol " << totNum << "\n" ;
            std::exit(1) ; 
        }

        delete pcrClassPtr ;
    }



}

void doPCRBulk(
    std::string& cellName , 
    uint32_t& dupCount,
    std::vector<util::CellBarcodeUMIBasicBlock>& uniqueMolecules ,
    std::vector<Transcript>& transcripts ,
    uint32_t& numOfPCRCycles ,
    double& errorRate ,
    bool& simulateFromIntrons,
    bool& switchOnEffModel,
    bool& isNoisyCell,
    bool& isDoublet,
    ErrorMatrix errorModel,
    fmt::MemoryWriter& sstream_left ,
    fmt::MemoryWriter& sstream_right,
    bool debug = false 
){

    //size_t totalMismatches{0} ;
    std::string addedString = isNoisyCell ? ":Noise" : "" ;
    addedString += (isDoublet) ? ":Doublet" : "" ;

    PCRClass* pcrClassPtr{nullptr} ;


    // NOTE add a global debug flag 
    bool debug2 = false ;// (cellName == "AGCGCCTTCCCCTGAT") ;
    //bool debug3 = (cellName == "GTATTCTGTAGTGAAT") ;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double efficiencyProb{0.98} ;
    double mutationProb{0.01} ;
    double seqeneceErrorProb{0.005} ;

    //fragmentLength
    size_t fragmentLength{MAX_FRAGLENGTH} ;

    //size_t MAXNUM = 100000 ;
    //double intronSimulationRate{0.02} ;

    auto totNum = std::pow(2, numOfPCRCycles) * uniqueMolecules.size() ;


    if(debug && debug2){
        std::cerr << cellName << "\n" ;
        std::cerr << "=======================================\n" ;
        std::cerr << "cell expression: " << uniqueMolecules.size() << "\n" ;
        std::cerr << "cell barcode to be simulated " << uniqueMolecules[0].cellBarCode << "\n" ;
        std::cerr << "totNum: " << totNum << "\n" ;
        std::cerr << "MAXNUM: " << MAXNUM << "\n" ;
        std::cerr << "=======================================\n" ;
    }
    
    // Create a map of sequences
    // create parent sequences
    struct fragInfo{
        std::string fragmentSeq ;
        size_t absStartPos ;
    } ;


    // This map holds the uinqueMolecule ID to sequence string 
    std::map<uint32_t, std::string> sequenceMap ;
    std::map<uint32_t, std::string> readSeqMap ;
    // This map holds the uniqueMolecule ID to a pair of positions
    // aka start position and end position where the the fragment 
    // start position lies
    std::unordered_map<uint32_t, std::pair<size_t, size_t>> startPositionMap ;

    // Based on a random variable decide if we would like 
    // to include the intronic sequence or not 
    bool intronInclusion{false} ;
    double toss ;

    int polyAcount{0} ;
    int polyAcount2{0} ;
    size_t uniqReadId{0} ;

    // Generate original sequences
    for(uint32_t i = 0; i < uniqueMolecules.size(); ++i){
        auto& tid = uniqueMolecules[i].transcriptId ;
        auto& tr = transcripts[tid] ;

        bool polyASample{false} ;

        std::string polyASeq{""} ;
        // NOTE: Ideally this should be 
        // a bool variable that represents that 
        // to sample from provided polyA ends or not 
        if(!uniqueMolecules[i].count){
            
            // If this option is enabled then the polyA end  
            // must already be present as the part of the 
            // reference.    

            auto& polyAVec = tr.loadPolyASeq() ;
            if(polyAVec.size() > 0){
                polyASample = true ;
                polyASeq = polyAVec.back() ;
            }
            polyAcount2++ ;
        }

        std::string readSeq ;
        uint32_t startPos ;
        std::string fragmentSeq ;
        uint32_t fragmentStartPos ;
        uint32_t fragmentEnd ;
        uint32_t fragLen = 0;

        if(polyASample) {
            // Sample from PolyASeq. 

            startPos = generateFragmentSeq(polyASeq, readSeq) ;
            fragmentStartPos = 0 ;
            fragmentLength = polyASeq.size() ;
            fragmentEnd = fragmentLength - 1;
            fragmentSeq = polyASeq ;
            polyAcount++ ;
        } else {
            // Normal transncript 
            // Decide if velocity mode is on or not 
            toss = distribution(generator) ;
            std::string intronSeq{""} ;
            if (simulateFromIntrons and (toss < 0.02)){
                intronInclusion = tr.loadUnplicedIntronSeq() ;
                    
            }
            if(tr.RefLength < MAX_FRAGLENGTH) {
                fragmentStartPos = 0 ;
                fragLen = tr.RefLength;
            }else {
                fragmentStartPos = tr.RefLength - MAX_FRAGLENGTH ;
                fragLen = MAX_FRAGLENGTH ;
            }
            fragmentEnd = fragmentStartPos + fragLen - 1 ;
            fragmentSeq = std::string{tr.Sequence() + fragmentStartPos, fragLen} ;
            //uint32_t intronFragmentLen = 0 ;

            //Intron 
            if (intronInclusion) {
                uint32_t lastExonLength{tr.getLastExonLength()} ;
                std::string exonSeq(tr.Sequence() + (tr.RefLength - lastExonLength), lastExonLength) ;
                fragmentSeq = intronSeq + exonSeq ;
                if (fragmentSeq.size() == 0) {
                    std::cerr << "DEBUG : lastExonLength "  << lastExonLength 
                              << "For transcript "<< tr.RefName <<"\n" ;
                    std::exit(1) ;
                }
                if(exonSeq.size() < MAX_FRAGLENGTH) {
                    fragmentStartPos = 0 ;
                    fragLen = fragmentSeq.size();
                    fragmentEnd = fragLen - 1 ;
                } else {
                    fragmentStartPos = exonSeq.size() - MAX_FRAGLENGTH ;
                    fragLen = MAX_FRAGLENGTH ;
                    fragmentEnd = fragmentStartPos + fragLen - 1 ;
                }
                startPos = generateFragmentSeq(fragmentSeq, readSeq) ;
                if(readSeq.size() != READ_LEN) {
                    std::cerr << "The readSeq is initilized wrongly \t" ;
                    std::cerr << "startPos " << startPos << "\n" ;
                }
            } else {
                    startPos = generateFragmentSeq(tr, readSeq) ;
            }
        }

        sequenceMap[i] = uniqueMolecules[i].cellBarCode + 
                         uniqueMolecules[i].UMICode + 
                         std::string(fragmentSeq) ;

        startPositionMap[i] = { fragmentStartPos, fragmentEnd };
        readSeqMap[i] = readSeq ;
        std::string tagSeq{sequenceMap[i].substr(0, CB_LENGTH + UMI_LENGTH)} ;
    }

    std::vector<util::EditBlock> editVector;
    editVector.resize(totNum) ;
    for(uint32_t i = 0 ; i < uniqueMolecules.size(); ++i){
        editVector[i] = {-1, -1, true} ;
    }



    bool doPCR{true} ;
    if(doPCR){
        pcrClassPtr = new PCRClass(
            uniqueMolecules.size(),
            numOfPCRCycles,
            mutationProb,
            efficiencyProb,
            errorRate,
            switchOnEffModel
        ) ;
        
        pcrClassPtr->runPCR(
            sequenceMap
        );
    }

    
    std::vector<uint32_t> indexArray ;
    indexArray.resize(totNum) ;
    for(uint32_t i = 0; i < totNum ; ++i){
        indexArray[i] = i ;
    }

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(indexArray.begin(), indexArray.end(), g) ;
    
    if(dupCount == 0){
        dupCount = MAXNUM ;
    }
    if (dupCount > totNum){
        dupCount = totNum ;
    }


    auto numMol = std::min(indexArray.size(), static_cast<size_t>(dupCount)) ;

    size_t numOfWritten{0} ;
    size_t blockSize = uniqueMolecules.size() ;
    
    if(doPCR){
        for(uint32_t i = 0 ; i < indexArray.size() ; ++i){
            auto ind = indexArray[i] ;
            if(ind < blockSize){
                numOfWritten++ ;

                sstream_left << "@" << cellName
                             << ":" << transcripts[uniqueMolecules[ind].transcriptId].RefName + addedString 
                             << ":" << 0
                             << ":" << i 
                             << ":" << uniqReadId 
                             << "\n" ;
                std::string tagSeq{sequenceMap[ind].substr(0, CB_LENGTH + UMI_LENGTH)} ;
                sstream_left << imputeErrorInStringUMI(seqeneceErrorProb, tagSeq) << "\n" ;
                sstream_left << "+\n" ;
                sstream_left << std::string(CB_LENGTH + UMI_LENGTH, 'N') << "\n" ; 

                sstream_right << "@" << cellName 
                              << ":" << transcripts[uniqueMolecules[ind].transcriptId].RefName + addedString
                              << ":" << 0 
                              << ":" << i 
                              << ":" << uniqReadId 
                              << "\n" ;
                
                if(errorModel.size() != 0){
                    sstream_right << imputeIlluminaModel(errorModel, readSeqMap[ind]) << "\n" ;
                }else{
                    sstream_right << imputeErrorInString(seqeneceErrorProb, readSeqMap[ind]) << "\n" ;
                }

                //sstream_right << imputeErrorInString(seqeneceErrorProb, readSeqMap[ind]) << "\n" ;
                sstream_right << "+\n" ;
                sstream_right << std::string(READ_LEN, 'N') << "\n" ; 
                uniqReadId++ ;              
                
                if(numOfWritten >= numMol)
                    break ;
                
                continue ;
            }



            bool validBlock = pcrClassPtr->constructOrCacheSequence(ind, sequenceMap) ;
            
            auto grandParentId = ind % blockSize ; 

            if(validBlock){

                auto it = sequenceMap.find(ind) ;
                if(it != sequenceMap.end()){
                    auto parentSeq = it->second ;
                    std::string readSeq ;
                    uint32_t startPos ;
                    auto fStart = startPositionMap[grandParentId].first ;
                    auto fEnd = startPositionMap[grandParentId].second ;
                    clipSequence(parentSeq, fStart, fEnd, readSeq, startPos) ;
                    auto absStartPos = fStart + startPos ;                

                    // left end            
                    sstream_left << "@" << cellName 
                                    << ":" << transcripts[uniqueMolecules[grandParentId].transcriptId].RefName + addedString
                                    << ":" << absStartPos
                                    << ":" << grandParentId
                                    << ":" << uniqReadId 
                                    << "\n" ;
                    std::string tagSeq{parentSeq.substr(0, CB_LENGTH + UMI_LENGTH)} ;
                    //sstream_left << tagSeq << "\n" ;
                    sstream_left << imputeErrorInStringUMI(seqeneceErrorProb,tagSeq) << "\n" ;

                    
                    // wrtie '+' and quality
                    sstream_left << "+\n" ;
                    sstream_left << std::string(CB_LENGTH + UMI_LENGTH, 'N') << "\n" ; 

                    // right end
                    sstream_right << "@" << cellName
                                    << ":" << transcripts[uniqueMolecules[grandParentId].transcriptId].RefName + addedString
                                    << ":" << absStartPos
                                    << ":" << grandParentId
                                    << ":" << uniqReadId 
                                    << "\n" ;
                    uniqReadId++ ;
                    
                    if(errorModel.size() != 0){
                        sstream_right << imputeIlluminaModel(errorModel, readSeq) << "\n" ;
                    }else{
                        sstream_right << imputeErrorInString(seqeneceErrorProb, readSeq) << "\n" ;
                    }

                    //sstream_right << imputeErrorInString(seqeneceErrorProb, readSeq) << "\n" ;
                    // wrtie + and quality
                    sstream_right << "+\n" ;
                    sstream_right << std::string(READ_LEN, 'N') << "\n" ; 
                        
                }else{
                    std::cerr << "Should not happen\n" ;
                    std::exit(1) ;
                }

                numOfWritten++ ;
            }
            if(numOfWritten >= numMol){
                break ;
            }
       
        }
        //pcrClassPtr->checkChainLength() ;
        delete pcrClassPtr ;
    }
}


template <typename MutexT, typename T>
void generateSequencesForCellDBG(
    std::atomic<uint32_t>& threadsDone,
    std::atomic<uint32_t>& ccount,
    uint32_t& numCells,
    const DataMatrix<T>& dataMatrixObj,
    uint32_t& numOfPCRCycles,
    double& errorRate,
    bool& switchOnEffModel,
    std::vector<std::string>& CBList,
    std::vector<std::string>& UMIList,
    Reference& refInfo,
    std::vector<Transcript>& transcripts,
    ErrorMatrix& errorModel,
    moodycamel::ConcurrentQueue<stream_pair_ptr, MyTraits>& conQueue,
    MutexT* iomutex
){
   uint32_t cellId{0} ;
   while((cellId = ccount++) < numCells){
        auto cellName = dataMatrixObj.getCellNameConst(cellId) ;
        auto cellBarCode = CBList[cellId] ;

        uint32_t dupCounts{0} ;
        // get diplicated counts from alevin
        if(dataMatrixObj.cellNamesDupCount.size() != 0){
            auto it = dataMatrixObj.cellNamesDupCount.find(cellName) ;
            if (it != dataMatrixObj.cellNamesDupCount.end()){
                dupCounts = it->second ;
            }
        }

       std::vector<util::CellBarcodeUMISegmentBlock> uniqueMolecules ;    
       std::vector<int> trueCellExpression = dataMatrixObj.fetchTrueCellExp(cellId) ;

        int totalNumberOfMolecules = std::accumulate(trueCellExpression.begin(), trueCellExpression.end(), 0) ;
        std::vector<int> sampledList ;
        
        indexDistribution(
            UMIList.size(),
            totalNumberOfMolecules,
            sampledList
        ) ;

        int cumId{0} ;
        //int cumExpressionCount{0} ;
        //uint32_t skippedExpCounts{0} ;

        auto dbgPtr = dataMatrixObj.dbgPtr ;
        auto& segMap = 
            dataMatrixObj.cellSegCount[cellId] ;


        //size_t missedSegments{0} ;
        for(auto git : segMap){
            auto segIdMap = git.second ;
            for(auto sit: segIdMap){
                bool ore{false} ;
                auto success = dataMatrixObj.getSegOreMapVectorInfo(git.first, sit.first, ore) ;
                if(!success){
                    std::cerr << "The key for gene id and unitig id not found "
                              << "gene id " << git.first << " unitig id " << sit.first << "\n"
                              << " post an issue with the gfa file in https://github.com/COMBINE-lab/minnow/issues \n" ;
                              std::exit(1) ; 
                }
                for(int j = 0 ; j < sit.second ; ++j){
                    // choose a random transcript 
                    // that contains this segment
                    bool present ;
                    auto fixTidInfo = dataMatrixObj.getRandomTrInfo(git.first, sit.first, refInfo, present) ;
                    if(!present){
                        std::cerr << "DBG ::: GHOST ALERT ::: spooky tid " << fixTidInfo.tid << "\n" ;
                        std::exit(1) ;
                    }

                    if((fixTidInfo.end - fixTidInfo.start) < READ_LEN){
                        std::cerr << " In minnow simulate length of contig is smaller than read length \n"
                                  << fixTidInfo.start << "\t" << fixTidInfo.end << "\n" ;
                    }

                    uniqueMolecules.emplace_back(
                        fixTidInfo.tid,
                        sit.first,
                        fixTidInfo.start,
                        fixTidInfo.end,
                        cellBarCode,
                        UMIList[sampledList[cumId]],
                        ore
                    );
                    cumId++ ;

                }
            }
        }

        //std::cerr << "Missed segments " << missedSegments << "\n" ;

        // PCR
        auto logger = spdlog::get("stderrLog");
        fmt::MemoryWriter sstream_left;
        fmt::MemoryWriter sstream_right;

        doPCRBulkDBG(
            cellName,
            dupCounts,
            uniqueMolecules,
            dataMatrixObj.rspdVec,
            transcripts,
            dbgPtr,
            numOfPCRCycles,
            errorRate,
            switchOnEffModel,
            errorModel,
            sstream_left,
            sstream_right
        ) ;



        // PCR done 
        if(iomutex->try_lock()){
            if(ccount < numCells){
                fmt::print(
                    stderr, "\rProcessed {} cells",ccount
                );
            }

            iomutex->unlock();
        }

        std::stringstream* ssLeft = new std::stringstream ;

        std::stringstream* ssRight = new std::stringstream ;
        
        std::unique_ptr<std::ostream> comLeft = 
            std::unique_ptr<std::ostream>( new zstr::ostream(*ssLeft)) ;
        std::unique_ptr<std::ostream> comRight = 
            std::unique_ptr<std::ostream>( new zstr::ostream(*ssRight)) ;


        *comLeft << sstream_left.str();
        comLeft->flush() ;
        sstream_left.clear();
        *comRight << sstream_right.str();
        sstream_right.clear();
        comRight->flush() ;

        comLeft.reset() ;
        comRight.reset() ;
        //std::cerr << "Size of stream_pair: " << outStr_left.size() << " " << outStr_right.size() << "\n" ;
        while((!conQueue.try_enqueue({ssLeft,ssRight}))){ 
            //std::cerr << "In while waiting for enque " << conQueue.size_approx() << "\n" ;
        }  

   }
   threadsDone++ ;
}


template <typename MutexT, typename T>
void generateSequencesForCell(
    std::atomic<uint32_t>& threadsDone,
    std::atomic<uint32_t>& ccount,
    uint32_t& numCells,
    const DataMatrix<T>& dataMatrixObj,
    uint32_t& numOfPCRCycles,
    double& errorRate,
    bool& simulateFromIntrons,
    bool& switchOnEffModel,
    std::vector<std::string>& CBList,
    std::vector<std::string>& UMIList,
    std::unordered_set<uint32_t>& emptyCellVector,
    std::vector<Transcript>& transcripts,
    ErrorMatrix& errorModel,
    moodycamel::ConcurrentQueue<stream_pair_ptr, MyTraits>& conQueue,
    MutexT* iomutex
)
{
    uint32_t cellId{0} ;

    while((cellId = ccount++) < numCells){
        if(emptyCellVector.find(cellId) != emptyCellVector.end()){
            continue ;
        }
        //std::cerr << "-------------> start nonsense " << cellId << "\n" ;

        auto cellName = dataMatrixObj.getCellNameConst(cellId) ;
        auto cellBarCode = CBList[cellId] ;

        bool isNoisyCell = (cellId < dataMatrixObj.numOfWhiteListedCells) ? false : true ;
        bool isDoublet = false ;

        if(isNoisyCell){
            if(dataMatrixObj.cellDoubletMap.find(cellName) != dataMatrixObj.cellDoubletMap.end()){
                isNoisyCell = false ;
                isDoublet = true ;
            }
        
        }

        uint32_t dupCounts{0} ;
        // get diplicated counts from alevin
        if(dataMatrixObj.cellNamesDupCount.size() != 0){
            auto it = dataMatrixObj.cellNamesDupCount.find(cellName) ;
            if (it != dataMatrixObj.cellNamesDupCount.end()){
                dupCounts = it->second ;
            }
        }
    

        bool debug{false} ; // = (cellName == "GTATTCTGTAGTGAAT") ;

        //std::cerr << "Cell barcode for this cell " << cellBarCode << "\n" ;

        std::vector<util::CellBarcodeUMIBasicBlock> uniqueMolecules ;         
        std::vector<T> cellExpressionVector = dataMatrixObj.fetchCellByIdConst(cellId) ;
        int totalNumberOfMolecules = std::accumulate(cellExpressionVector.begin(), cellExpressionVector.end(), 0) ;
        
        int totalIntronCount{0} ;
        auto& intronCountMap = dataMatrixObj.intronCountMap ;
        std::vector<std::pair<uint32_t, int>> intronExpressionVector ;
        if (intronCountMap.size() > 0){
            auto it = intronCountMap.find(cellId) ;
            if(it != intronCountMap.end()){
                intronExpressionVector = it->second ;
                for (auto itpair : intronExpressionVector){
                    totalIntronCount += itpair.second ; 
                }
            }
        }
    


        if (debug) {std::cerr << "Total number of molecules " << totalNumberOfMolecules << "\n" ; }

        // Generate random set of UMIs to match with 
        // totalNumberOfMolecules
        std::vector<int> sampledList ;
        
        indexDistribution(
            UMIList.size(),
            totalNumberOfMolecules + totalIntronCount,
            sampledList
        ) ;

        int cumId{0} ;
        int cumExpressionCount{0} ;
        uint32_t skippedExpCounts{0} ;

        for(size_t transcriptId = 0; transcriptId < cellExpressionVector.size(); ++transcriptId){
            const auto refTid = dataMatrixObj.getrefId(transcriptId) ;
            if(refTid ==  std::numeric_limits<double>::max()){
                std::cerr<<"FATAL ERROR: transcript id "<< transcriptId <<" is not present in reference exiting !! "<<
                                  "\nput an issue in https://github.com/COMBINE-lab/minnow/issues\n" ;
            }
            //std::cerr << "remapped transcript id: " << refTid << "\t" << transcripts[refTid].RefLength << "\n" ;
            if((cellExpressionVector[transcriptId] == 0) || (transcripts[refTid].RefLength < READ_LEN)){
                //if  (transcripts[refTid].RefLength < 50){
                //    if(cellExpressionVector[transcriptId] != 0){
                //        std::cerr << "Short transcript has non zero value \n" ;
                //    }
                //}
                skippedExpCounts += cellExpressionVector[transcriptId] ;
                continue ;
            }
            auto expressionCount = cellExpressionVector[transcriptId] ;
            for(int moleculeId = 0; moleculeId < expressionCount; ++moleculeId){
                
                uniqueMolecules.emplace_back(
                    refTid, 
                    cellBarCode, 
                    UMIList[sampledList[cumId]],
                    true,
                    0,
                    0,
                    false 
                ) ;
                cumId++ ;
            }
            cumExpressionCount += expressionCount ; 
        }


        // If this cell is marked for polyA clipping then we need to take that sequence too 
        int droppedIntronCount{0} ;
        if (totalIntronCount > 0){
            for(auto& tcpair : intronExpressionVector){
                auto tid = tcpair.first ;
                auto& icount = tcpair.second ;
                //const auto refTid = dataMatrixObj.getrefId(tid) ;
                //auto refTid = dataMatrixObj.alevin2refTranscriptMap[tid] ;
                if((icount == 0) || (transcripts[tid].RefLength < READ_LEN)){
                    continue ;
                }
                for(int moleculeId = 0 ; moleculeId < icount ; ++moleculeId){
                    uniqueMolecules.emplace_back(
                        tid, 
                        cellBarCode,
                        UMIList[sampledList[cumId]],
                        false,// this signifies the molecule comes from intronic polyA
                        0,
                        0,
                        false
                    ) ;
                    cumId++ ;
                    droppedIntronCount++ ;
                }
                
            }
        }


        //std::cerr << "Total number of molecules " << totalNumberOfMolecules << " Exon counts " << uniqueMolecules.size() << "\n" ; 

        if(debug) {
            std::cerr << "====================================\n" ;
            std::cerr << "Total Intron count for " << cellName << " \t " << totalIntronCount 
                       << "After dropping " << droppedIntronCount << "\n";
        }

        //std::cerr << "\tcell name: " << cellName << " skipped " << skippedExpCounts << "\n" ;
       
        // PCR
        auto logger = spdlog::get("stderrLog");
        fmt::MemoryWriter sstream_left;
        fmt::MemoryWriter sstream_right;
        //fmt::MemoryWriter sstream_matrix;
        //std::cerr << "\n Before calling bulk \n" ;

        doPCRBulk(
            cellName, 
            dupCounts,
            uniqueMolecules, 
            transcripts, 
            numOfPCRCycles, 
            errorRate, 
            simulateFromIntrons,
            switchOnEffModel,
            isNoisyCell,
            isDoublet,
            errorModel,
            sstream_left, 
            sstream_right
           // sstream_matrix
        ) ;
        // PCR done 
        if(iomutex->try_lock()){
            if(ccount < numCells){
                fmt::print(
                    stderr, "\rProcessed {} cells",ccount
                );
            }

            iomutex->unlock();
        }
        // std::cerr << "\n okay until here \n" ;

        //std::unique_ptr<std::ostream> comLeft = 
        //    std::unique_ptr<std::ostream>( new zstr::ostream(*(streamVector[cellId].first))) ;
        //std::unique_ptr<std::ostream> comRight = 
        //    std::unique_ptr<std::ostream>( new zstr::ostream(*(streamVector[cellId].second))) ;


        std::stringstream* ssLeft = new std::stringstream ;

        std::stringstream* ssRight = new std::stringstream ;
        
        std::unique_ptr<std::ostream> comLeft = 
            std::unique_ptr<std::ostream>( new zstr::ostream(*ssLeft)) ;
        std::unique_ptr<std::ostream> comRight = 
            std::unique_ptr<std::ostream>( new zstr::ostream(*ssRight)) ;

        if(comLeft == nullptr){
            //std::cerr << " comLeft is nullPtr \n\n " ;
        }
        else{
            //std::cerr << "not null \n" ;
        }

        //comLeft->seekp(0, std::ios::end);
        //std::stringstream::pos_type offset = comLeft->tellp();
        //std::cerr << " TOKEN: " << token << " THREAD DEBUG: " << offset << std::endl ;

        *comLeft << sstream_left.str();
        comLeft->flush() ;
        sstream_left.clear();
        *comRight << sstream_right.str();
        sstream_right.clear();
        comRight->flush() ;

        comLeft.reset() ;
        comRight.reset() ;
        //std::cerr << "Size of stream_pair: " << outStr_left.size() << " " << outStr_right.size() << "\n" ;
        while((!conQueue.try_enqueue({ssLeft,ssRight}))){ 
            //std::cerr << "In while waiting for enque " << conQueue.size_approx() << "\n" ;
        } 
        //std::cerr << "waiting for this thread to finish " << conQueue.size_approx() << "\n" ;

        
    }
    threadsDone++;
    //std::cerr << "--------------> end nonsense \n" ;
}

void writeSequences(
    std::atomic<uint32_t>& threadsDone,
    std::ofstream& outLeft,
    std::ofstream& outRight,
    moodycamel::ConcurrentQueue<stream_pair_ptr, MyTraits>& conQueue,
    uint32_t& nthread 
){
        
    stream_pair_ptr pptr ;
    while((threadsDone < nthread-1)){
            while(conQueue.try_dequeue(pptr)){

                outLeft.write(pptr.first->str().data(), pptr.first->str().length()) ;
                pptr.first->clear() ;
                delete pptr.first ;
                outRight.write(pptr.second->str().data(), pptr.second->str().length()) ;
                pptr.second->clear() ;
                delete pptr.second ;
            }
    }
    
}


template<typename T>
bool spawnCellThreads(
    uint32_t& nthread,
    DataMatrix<T>& dataMatrixObj,
    std::string& outDir, 
    std::vector<std::string>& CBList,
    std::vector<std::string>& UMIList,
    std::unordered_set<uint32_t>& emptyCellVector,
    Reference& refInfo,
    std::vector<Transcript>& transcripts,
    SimulateOptions& simOpts,
    std::shared_ptr<spdlog::logger> consoleLog
){

    
    uint32_t numOfCells{dataMatrixObj.numCells} ;
    bool simulateFromIntrons = simOpts.velocityMode ;
    bool switchOnEffModel = simOpts.switchOnEffModel ;

    bool useDBG = simOpts.useDBG ;

    std::string illuminaModelFile = simOpts.illuminaModelFile ;

    ErrorMatrix errorModel ;

    if(illuminaModelFile != ""){
        util::readIlluminaErrorModel(illuminaModelFile, errorModel) ;
        consoleLog->info("Illumuna error model is loaded") ;
    }


    //concurrent queue to store the simulated reads
    moodycamel::ConcurrentQueue<stream_pair_ptr, MyTraits> conQueue(1, 0, nthread -1 ) ;  

    // create the output directory if it's not there 
    auto errorRate = simOpts.errorRate ;
    auto numOfPCRCycles = simOpts.numOfPCRCycles ;

    // set up the outout buffer
    std::ofstream outLeft(outDir + "/hg_100_S1_L001_R1_001.fastq.gz", std::ios::binary) ;
    std::ofstream outRight(outDir + "/hg_100_S1_L001_R2_001.fastq.gz", std::ios::binary) ;

    //std::ofstream outSkippedExpreesion(outDir + 'skipped_expression.txt') ;
    
    // set up the outout buffer

    // take care of thread safety 
    SpinLockT iomutex;
    consoleLog->info("Starting Minnow....") ;
    consoleLog->info("Number of cells to be processed {}", numOfCells) ;

    std::atomic<uint32_t> ccount{0} ;
    std::vector<std::thread> workerThreads ;

    std::atomic<uint32_t> threadsDone{0};
    for(size_t tn = 0; tn < nthread -1 ; ++tn){

        if(!useDBG){
            workerThreads.emplace_back(
                generateSequencesForCell<SpinLockT,T>,
                std::ref(threadsDone),
                std::ref(ccount),
                std::ref(numOfCells),
                std::ref(dataMatrixObj),
                std::ref(numOfPCRCycles),
                std::ref(errorRate),
                std::ref(simulateFromIntrons),
                std::ref(switchOnEffModel),
                std::ref(CBList),
                std::ref(UMIList),
                std::ref(emptyCellVector),
                std::ref(transcripts),
                std::ref(errorModel),
                std::ref(conQueue),
                &iomutex
            );

        }else{
            workerThreads.emplace_back(
                generateSequencesForCellDBG<SpinLockT,T>,
                std::ref(threadsDone),
                std::ref(ccount),
                std::ref(numOfCells),
                std::ref(dataMatrixObj),
                std::ref(numOfPCRCycles),
                std::ref(errorRate),
                std::ref(switchOnEffModel),
                std::ref(CBList),
                std::ref(UMIList),
                std::ref(refInfo),
                std::ref(transcripts),
                std::ref(errorModel),
                std::ref(conQueue),
                &iomutex
            );

        }
    }

    workerThreads.emplace_back(
        writeSequences,
        std::ref(threadsDone),
        std::ref(outLeft),
        std::ref(outRight),
        std::ref(conQueue),
        std::ref(nthread)
    );

    string_pair str_pair ;
    
    

    for(auto& t : workerThreads){
        t.join() ;
    }
    
    // flush remaining stuff 
    stream_pair_ptr pptr;
    while(conQueue.try_dequeue(pptr)){

                outLeft.write(pptr.first->str().data(), pptr.first->str().length()) ;
                pptr.first->clear() ;
                delete pptr.first ;
                outRight.write(pptr.second->str().data(), pptr.second->str().length()) ;
                pptr.second->clear() ;
                delete pptr.second ;
    }
    std::cout << "Size of conQueue: " << conQueue.size_approx() << "\n" ;
    
    consoleLog->info("\n\n") ;
    consoleLog->info("Done simulating reads ");
    //outFile_left.close();
    //outFile_right.close();
    
    return true ;
}


void printOptions(SimulateOptions& simOpts){
    std::cerr << "Matrix File Name "  << simOpts.matrixFile << "\n" 
              << "Reference Fasta " << simOpts.refFile << "\n" 
              << "Number of PCR cycles " << simOpts.numOfPCRCycles << "\n" 
              << "Erorr rate " << simOpts.errorRate << "\n" 
              << "Numeber of threads " << simOpts.numThreads << "\n" ; 
}


void minnowSimulate(SimulateOptions& simOpts){

    printOptions(simOpts) ;

    //bool alevinMode = simOpts.alevinMode ;
    //bool splatterMode = simOpts.splatterMode ;
    //auto& matrixFileName = simOpts.matrixFile ;
    auto& refFileName = simOpts.refFile ;
    auto& gene2txpFile = simOpts.gene2txpFile ; 

    auto& numOfCells = simOpts.numOfCells ;
    auto& numOfGenes = simOpts.numOfTranscripts ;

    //auto& numOfSampleCells = simOpts.sampleCells ;
    //auto& numOfPCRCycles = simOpts.numOfPCRCycles ;
    //auto& errorRate = simOpts.errorRate ;
    auto& numThreads = simOpts.numThreads ;

    std::string outDir = simOpts.outDir ;
    std::string intronFile = simOpts.intronFile ;

    bool useDBG = simOpts.useDBG ;

    //auto& sampleCells = simOpts.sampleCells ;

     
    //numOfPCRCycles = 10 ;

    // Set up the logger 
    auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>() ;
    auto consoleLog = spdlog::create("minnow-Log", {consoleSink});


    // Collect reference information
    // this will populate the transcript 
    // sequences too 

    
    consoleLog->info("Reading reference sequences ...") ;
    
    Reference refInfo(
        refFileName,
        gene2txpFile
    ) ;
    consoleLog->info("Reference sequence is loaded ...") ;

    // If the velocity mode is on we need to initialize 
    // One additional file with intronic reads from the 
    // fasta file 
    if(simOpts.velocityMode){
        consoleLog->info("RNA velocity mode is on; reading the intronic reads ...") ;
        refInfo.updateIntronSequence(simOpts.intronFile) ;
        consoleLog->info("Done reading intron files...") ;
        consoleLog->info("Reading the exon length sizes...") ;

        refInfo.updateLastExonLevelInfo(simOpts.exonLengthFile) ;
        
    }

    // Read the matrix
    // create explicit pointers 
    // remember to delete them after operation 
    // to avoid memory leaks 
     
    if (outDir.back() == '/'){
        outDir.pop_back();
    }

    // needed for vpolo
    std::string alevinLikeDir = outDir + "/alevin" ;
    util::fs::MakeDir(outDir.c_str()) ;
    if(!util::fs::DirExists(outDir.c_str())){
      consoleLog->error("we can not create nested direcory {} if that is what is intended", outDir) ;
      std::exit(2) ;
    }


    util::fs::MakeDir(alevinLikeDir.c_str()) ;

    
    DataMatrix<double> dataMatrixObj(
        simOpts,
        refInfo,
        consoleLog
    ) ;

    consoleLog->info("Done parsing matrix") ;


    uint32_t numOfWhiteListedCells ;
    uint32_t numOfNoisyCells ;
    uint32_t numOfDoublets ;
    //if (alevinMode ||  splatterMode || normalMode){
        numOfCells = dataMatrixObj.numCells ;
        numOfGenes = dataMatrixObj.numOfTranscripts ;
        numOfWhiteListedCells = dataMatrixObj.numOfWhiteListedCells ;
        numOfNoisyCells = dataMatrixObj.numOfNoisyCells ;
        numOfDoublets = dataMatrixObj.numOfDoublets ;
    //}


    //check the parsed matrix 
    //std::cerr << "Number of rows " <<  dataMatrixObj.data.size() << "\n" ; 
    //std::cerr << "Number of columns " <<  dataMatrixObj.data[0].size() << "\n" ; 

    // check for empty droplets according to Lun et al'17
    // algorithm 
    // std::vector<bool> emptyDropLabels ;
    // calculateEmptyDrops(dataMatrixObj, emptyDropLabels)


    // generate cell barcode list 
    

    // If we use usedbg this is useless 
    
    consoleLog->info("Checking for 0 expressed cells \n") ;
    uint32_t numSkipped{0} ;
    std::unordered_set<uint32_t> emptyCellVector ;


    if(!useDBG){
        for(size_t i = 0; i < dataMatrixObj.geneCounts.size() ; ++i){
            
            auto& v = dataMatrixObj.geneCounts[i] ;
            auto& t = dataMatrixObj.data[i] ;
            auto cellSum = std::accumulate(v.begin(), v.end(), 0) ;
            auto cellSumT = std::accumulate(t.begin(), t.end(), 0) ;

            
            if (cellSum == 0 or cellSumT == 0){
                numSkipped++ ;
                emptyCellVector.insert(i) ;
                //std::cout << dataMatrixObj.cellNames[i] << " .. skipping .." << "\n" ;
                
                //std::exit(1) ;
            }
            _verbose("\rNumber of cells checked : %lu", i);

        }
    }



    //consoleLog->info("\n{} cells will be skipped for having 0 expression in all genes",numSkipped) ;
    consoleLog->info("\n{} emptyCell size",emptyCellVector.size()) ;



    std::vector<std::string> CBList = util::generate10XCBList(numOfCells) ;

    
    consoleLog->info("CBList.size(): {}",CBList.size()) ;
    consoleLog->info("Number of cells {}",numOfCells) ;
    consoleLog->info("Number of whiteListedCells {}",numOfWhiteListedCells) ;
    consoleLog->info("Number of noisyCells {}",numOfNoisyCells) ;
    consoleLog->info("Number of Doublets {}",numOfDoublets) ;


    // don't write cells with all 0 count, it is unnecessary
    // we can log them if needed, now dropping them on floor
    // Leeroy Jenkins!!!!

    bool noDump{simOpts.noDump} ;
    
    std::ofstream whiteListStream(outDir + "/minnow_whitelist.txt") ;
    std::ofstream quantMatRowStream(alevinLikeDir + "/quants_mat_rows.txt") ;
    
    
    
    // read everything in quants_mat_rows first this is instrumental 
    // for reading the matrix bad, don't skip this 
    for(uint32_t cellId = 0; cellId < numOfCells ; ++cellId){
        if(emptyCellVector.find(cellId) == emptyCellVector.end()){
            quantMatRowStream << CBList[cellId] << "\n" ;
        }else{
            std::cerr << "this is empty\n" ; 
        }
    }



    uint32_t numOfCellsWritten{0} ; 
    for(uint32_t cellId = 0; cellId < numOfWhiteListedCells ; ++cellId){
        if (CBList[cellId] == ""){
            std::cerr << "whitelist is empty\n" ;
			std::exit(1) ; 
        }
        if(emptyCellVector.find(cellId) == emptyCellVector.end())
            whiteListStream << CBList[cellId] << "\n" ;
        if(numOfCells < numOfWhiteListedCells){
            if(cellId == numOfCells-1)
                break ;
        }
        numOfCellsWritten++ ;
    }

    std::ofstream noisyListStream(outDir + "/noisy_cells.txt") ;

    for(uint32_t cellId = numOfCellsWritten; (cellId < (numOfCellsWritten + numOfNoisyCells)) && (cellId < numOfCells) ; ++cellId){
        if (CBList[cellId] == ""){
            std::cerr << "whitelist is empty\n" ;
			std::exit(1) ; 
        }
        if (emptyCellVector.find(cellId) == emptyCellVector.end()){
            noisyListStream << CBList[cellId] << "\n" ;
        }
        numOfCellsWritten++ ;
    }
    std::ofstream doubletStream(outDir + "/doublet_cells.txt") ;
    for(uint32_t cellId = numOfCellsWritten; (cellId < (numOfCellsWritten + numOfDoublets)) && (cellId < numOfCells) ; ++cellId){
        if (CBList[cellId] == ""){
            std::cerr << "whitelist is empty\n" ;
			std::exit(1) ; 
        }
        if (emptyCellVector.find(cellId) == emptyCellVector.end()){
            doubletStream << CBList[cellId] << "\n" ;
        }
        numOfCellsWritten++ ;
    }
    


    consoleLog->info("Dumping the truth related files to the disk, this can take considerable time with respect to size, you can avoid it by --noDump") ;
    std::string geneCountFileName = alevinLikeDir + "/quants_mat.csv" ;
    std::string geneCountFileNameOriginal = alevinLikeDir + "/gene_count_original.txt" ;
    std::string cellNamesFile = alevinLikeDir + "/true_cell_names.txt" ;

    //
    if(!noDump){
        dataMatrixObj.dumpGeneCountMatrix(geneCountFileName, emptyCellVector) ;
        dataMatrixObj.dumpTrueCellNames(cellNamesFile) ;

        // DEBUG dump intron counts 
        //dataMatrixObj.dumpIntronCount() ;
        consoleLog->info("Truth files dumped") ;

    }



    // generate UMI list 
    std::vector<std::string> UMIList = util::generateUMIList() ; 

    // If we usedbg none of other stuff needed

    bool success = spawnCellThreads(
        numThreads,
        dataMatrixObj,
        outDir,
        CBList,
        UMIList,
        emptyCellVector,
        refInfo,
        refInfo.transcripts,
        simOpts,
        consoleLog
    ) ;

    dataMatrixObj.stop() ;
    if(success) {consoleLog->info("Successfully written simulated reads");} 

} 
