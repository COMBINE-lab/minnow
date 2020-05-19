#ifndef REFERENCE_INFO_HPP
#define REFERENCE_INFO_HPP

#include <iostream>
#include <numeric>
#include <random>
#include <limits>
#include <vector>

#include "Transcript.hpp"
#include "MinnowUtil.hpp"
#include "MinnowFS.hpp"
#include "FASTAParser.hpp"

class Reference{
    public:
        Reference(
            std::string& fastaFileIn,
            std::string& gene2txpFileIn,
            std::shared_ptr<spdlog::logger>& consoleLogIn

        ){
            gene2txpFile = gene2txpFileIn ;
            fastaFile = fastaFileIn ;
            FASTAParser fastaParser(fastaFile) ;
            fastaParser.populateTargets(transcripts) ;
            size_t trId{0} ;
            for(auto& tr : transcripts){
                transcriptNameMap[tr.RefName] = trId++ ;
            }

            numOfTranscripts = transcripts.size();
            consoleLog = consoleLogIn ;

            consoleLog->info("Transcript file {} is read", fastaFileIn);

        }    

        void inline updateIntronSequence(std::string& intronfname_){
            FASTAParser fastaParser(intronfname_) ;
            if (transcripts.size() == 0){
                std::cerr << "This function should be called after populateTargets\n" ;
                std::exit(1) ;
            }
            fastaParser.updateTranscriptLevelIntron(transcripts, transcriptNameMap) ;
        }

        void inline updateLastExonLevelInfo(std::string& exonLengthFileName_){

            if (transcripts.size() == 0){
                std::cerr << "This function should be called after populateTargets\n" ;
                std::exit(1) ;
            }
            std::ifstream exonLengthFileHandle(exonLengthFileName_.c_str()) ;
            std::string line ;

            while(std::getline(exonLengthFileHandle,line)){
                std::stringstream lineStream(line) ;
                std::vector<std::string> tokens ;

                util::split(line, tokens, "\t") ;
                
                auto trName = tokens[0] ;
                auto lastExonLength = std::stoul(tokens[1]) ;

                auto it = transcriptNameMap.find(trName) ; 
                if (it != transcriptNameMap.end()){
                    auto trId = it->second ;
                    transcripts[trId].setLastExonLength(lastExonLength) ;
                }

            }
        }

        void inline updatePolyASequence(
            std::string& polyAFileName_,
            std::unordered_map<std::string, uint32_t>& geneMap,
            std::unordered_map<uint32_t, uint32_t>& gene2LastTrMap,
            std::unordered_set<uint32_t>& geneIntronMap
        ){
            FASTAParser fastaParser(polyAFileName_) ;
            if (transcripts.size() == 0){
                std::cerr << "This function should be called after populateTargets\n" ;
                std::exit(1) ;
            }
            fastaParser.updateGeneLevelIntron(
                transcripts,
                geneMap,
                gene2LastTrMap, 
                geneIntronMap
            ) ;

            int totT{0} ;
            for(auto& tr : transcripts){
                auto& sv = tr.loadPolyASeq() ; 
                if (sv.size() > 0){
                    totT++ ;
                    //std::cerr << tr.RefName << "\t" << sv.size() << "\n" ;
                }
            }
            consoleLog->info("Total nonzero polyA transcrtipts {}",totT) ;

        }
        
        void inline updateGene2TxpMap(std::string& gene2txpFile){
                
                std::ifstream gene2txpStream(gene2txpFile.c_str()) ;
                std::string line ;

                uint32_t geneId{0} ;
                uint32_t skippedTranscripts{0} ;

                while(std::getline(gene2txpStream, line)){
                    std::stringstream lineStream(line) ;
                    std::vector<std::string> valueOfCells ; 
                    //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    util::split(line, valueOfCells, "\t") ;

                    
                    // get the transcript name 
                    auto transcriptName = valueOfCells[0] ;
                    auto geneName = valueOfCells[1] ;

                    auto trIt = transcriptNameMap.find(transcriptName) ;
                    if( trIt != transcriptNameMap.end()){
                        auto geneIt = geneMap.find(geneName) ;
                        if (geneIt != geneMap.end()){
                            gene2transcriptMap[geneIt->second].push_back(trIt->second) ;
                            transcript2geneMap[trIt->second] = geneIt->second ;
                        }else{
                            geneMap[geneName] = geneId ;
                            gene2transcriptMap[geneId] = {trIt->second} ;
                            transcript2geneMap[trIt->second] = geneId ;
                            geneId++ ;
                        }
                    }else{
                        skippedTranscripts++ ;
                        //std::cerr << "The transcript " << transcriptName << " is not present in the reference \n" ;
                        //std::cerr << "skipping this row\n" ;
                        //std::exit(1) ;
                    }
                }
                //std::cerr << "Skipped " << skippedTranscripts << " transcripts because either short or not present in reference \n" ;
                consoleLog->info("Skipped {} transcripts" 
                                  "because either short or not present in reference ",
                                  skippedTranscripts) ;
                //std::cerr << "Done with parsing " << gene2txpFile << "\n" ;
                //std::cerr << "#of genes " << gene2transcriptMap.size() << "\n" ;
        }

        inline uint32_t getGeneId(std::string geneName){
            auto it = geneMap.find(geneName) ;
            if(it != geneMap.end()){
                return it->second ;
            }else{
                return std::numeric_limits<uint32_t>::max() ;
            }
        }

        inline void updateDuplicateMap(std::string filename){
            std::ifstream dupStream(filename.c_str()) ;
            std::string line ;

            while(std::getline(dupStream, line)){
                std::stringstream lineStream(line) ;
                std::vector<std::string> valueOfCells ; 
                //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                util::split(line, valueOfCells, "\t") ;

                
                // get the transcript name 
                auto transcriptKept = valueOfCells[0] ;
                auto transcriptRemoved = valueOfCells[1] ;

                auto trIt = transcriptNameMap.find(transcriptKept) ;
                auto trIt2 = transcriptNameMap.find(transcriptRemoved) ;
                if( trIt != transcriptNameMap.end() && trIt2 != transcriptNameMap.end() ){
                    duplicateList[trIt2->second] = trIt->second ;
                }
            }

        }

    // private:
        std::string fastaFile ;
        std::string gene2txpFile ;

        
        std::vector<Transcript> transcripts; 
        std::unordered_map<uint32_t, uint32_t> duplicateList ;
        std::unordered_map<std::string, uint32_t> geneMap ;
        std::map<uint32_t, std::vector<uint32_t>> gene2transcriptMap ;
        std::unordered_map<uint32_t, uint32_t> transcript2geneMap ;
        std::unordered_map<std::string, uint32_t> transcriptNameMap ;
        size_t numOfTranscripts ;
        std::shared_ptr<spdlog::logger> consoleLog;

} ;

#endif