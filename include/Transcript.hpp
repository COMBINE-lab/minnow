#ifndef TRANSCRIPT_HPP
#define TRANSCRIPT_HPP

#include <iostream>
#include <numeric>
#include <memory>
#include <vector>
#include "macros.hpp"

class Transcript{
    public:
        Transcript(size_t idIn, std::string name, uint32_t len) : 
        RefName(name),
        RefLength(len),
        id(idIn)
        {} 

        // Will *not* delete seq on destruction
        /*
        void setSequenceBorrowed(
            const char* seq, 
            bool needGC = false
        ){
            Sequence_ = std::unique_ptr<const char, void (*)(const char*)>(
                seq,                 // store seq
                [](const char* p) {} // do nothing deleter
            );
        }*/

        // Will delete seq on destruction
        void setSequenceOwned(
            const char* seq
        ){
            Sequence_ = std::unique_ptr<const char, void (*)(const char*)>(
                seq,                              // store seq
                [](const char* p) { delete[] p; } // do nothing deleter
            );
    
        }
        
        // Will delete seq on destruction
        void setIntronSequenceOwned(
            const char* seq
        ){
            intronSequence_ = std::unique_ptr<const char, void (*)(const char*)>(
                seq,                              // store seq
                [](const char* p) { delete[] p; } // do nothing deleter
            );
            hasIntron = true;
        }

        // store coordinate and intron sequence
        void insertIntronSeq(uint32_t chrEndPos, std::string& seq){
            intronSeq.push_back( std::make_pair(chrEndPos, seq) ) ;
        }

        // load a random Intron sequence back 
        bool loadUnplicedIntronSeq(){
            if(lastExonLength > 404){
                return false;
            }else{
                if(intronSeq.size() > 0){
                    intronSeq.back().second ;
                }
                return true  ;
            }

        }


        // set the length of last exon sequence 
        void setLastExonLength(uint32_t l){
            lastExonLength = l ;
        } 
        
        // get the length of last exon sequence 
        uint32_t getLastExonLength(){
            return lastExonLength ;
        } 

        void insertPolyASeq(std::string& seq){
            polyASeq.push_back(seq) ;
        }

        void setGeneId(uint32_t& geneIdIn){
            geneId = geneIdIn;
        }

        uint32_t getGeneId(){
            return geneId;
        }

        // load intron sequence
        std::vector<std::pair<uint32_t,std::string>>& loadIntronSeq(){
            return intronSeq ;
        }
        std::vector<std::string>& loadPolyASeq(){
            return polyASeq ;
        }


        // get Sequence 
        const char* Sequence() const { return Sequence_.get(); }
        
        // get Sequence 
        const char* intronSequence() const { return intronSequence_.get(); }

        std::string RefName;
        uint32_t RefLength;
        uint32_t CompleteLength;
        uint32_t id;
        uint32_t geneId{std::numeric_limits<uint32_t>::max()};
        bool hasIntron{false};
    private:
        std::unique_ptr<const char, void (*)(const char*)> Sequence_ =
            std::unique_ptr<const char, void (*)(const char*)>(nullptr,
                                                         [](const char*) {});

        std::unique_ptr<const char, void (*)(const char*)> intronSequence_ =
            std::unique_ptr<const char, void (*)(const char*)>(nullptr,
                                                         [](const char*) {});

        std::vector<std::pair<uint32_t,std::string>> intronSeq ; // A vector of intron position and the sequences 
        std::vector<std::string> polyASeq ;
        
        uint32_t lastExonLength ;
};

#endif // TRANSCRIPT_HPP