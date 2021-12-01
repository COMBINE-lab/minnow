#ifndef MATRIX_PARSER
#define MATRIX_PARSER

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream> 
#include <vector>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "MinnowUtil.hpp"
#include "ReferenceInfo.hpp"
#include "ProgOpts.hpp"
#include "BFHClass.hpp"
#include "GFAReader.hpp"
#include "macros.hpp"


template <class T> 
class DataMatrix{
public :
	// read from file
	DataMatrix(
		SimulateOptions& simOpts,
		Reference& refInfo,
		std::shared_ptr<spdlog::logger>& consoleLogIn
	);

	inline void stop(){
		delete dbgPtr ;
		delete eqClassPtr ;
	}

	struct trInfo{
		uint32_t tid ;
		uint32_t start ;
		uint32_t end ; 

		trInfo() {}

		trInfo(
			uint32_t tidIn,
			uint32_t startIn,
			uint32_t endIn
		) :
		tid(tidIn),
		start(startIn),
		end(endIn)
		{} 

	} ;

	void loadAlevinData(	
		SimulateOptions& simOpts,
		Reference& refInfo 
	);

	void loadSplatterData(	
		SimulateOptions& simOpts,
		Reference& refInfo 
	);

	void dumpGeneCountMatrix(
		std::string& fileName,
		std::unordered_set<uint32_t>& emptyCells
	);
    
	void dumpTrueCellNames(
		std::string& fileName
	);

	void dumpGeneCountMatrixUnmolested(
		std::string& fileName
	);

	void dumpIntronCount() ;

	inline std::vector<std::vector<T>> fetchData(){
		return data ;
	}

	inline std::string getCellName(int cellId){
		return cellNames[cellId] ;
	}

	inline std::vector<T> fetchCellById(int cellId){
		return data[cellId] ;
	}

	inline std::string getCellNameConst (int cellId) const{
		return cellNames[cellId] ;
	}

	inline std::vector<T> fetchCellByIdConst (int cellId) const {
		return data[cellId] ;
	}
	
	inline std::vector<int> fetchTrueCellExp(int cellId) const {
		return trueGeneCounts[cellId] ;
	}

	inline bool getSegOreMapVectorInfo(const size_t gid, size_t sid, bool& ore) const{
		auto it = preSegOreMapVector[gid].find(sid);
		if(it != preSegOreMapVector[gid].end()){
			ore = it->second ;
			return true ;
		}
		return false ; 

	}

	inline trInfo getRandomTrInfo(const size_t gid, size_t sid,  Reference& refInfo, bool& present) const{
		auto it = geneSpecificTrVector[gid].find(sid) ;
		if(it != geneSpecificTrVector[gid].end()){
			auto trVec = it->second ;
			size_t checked{0} ;
			while(checked < trVec.size()){
				int ind = rand() % trVec.size() ;
				auto& fixTrInfo = trVec[ind] ; 

				if((fixTrInfo.end - fixTrInfo.start) < READ_LEN){
					std::cerr << "Should not happend REPORT IT \n" ;
					std::exit(1) ;
				}

				int mid = fixTrInfo.start + std::round((fixTrInfo.end - fixTrInfo.start)/2) ;
				if((refInfo.transcripts[fixTrInfo.tid].RefLength - mid) < READ_LEN){
					checked++;
					continue ;
				}
				present = true ;
				return trVec[ind] ;
			}
			int ind = rand() % trVec.size() ;
			present = true ;
			return trVec[ind] ;
		}else{
			present = false ;
			//std::cerr << "this pair does not exist in equivalence class !!!" ;
			return trInfo() ;
		}
	}

	inline uint32_t getrefId(uint32_t tid) const{
		
		auto it = alevin2refTranscriptMap.find(tid) ;
		if (it != alevin2refTranscriptMap.end())
			return it->second ;
		
		return  std::numeric_limits<uint32_t>::max() ;
	}


	// data structure related to  
	std::vector<std::vector<T>> data ; // Matrix containing the Cell x Transcriptome Matrix 
	std::vector<std::vector<T>> geneCounts ; // Matrix containing the Cell x Gene Matrix 
	std::vector<std::vector<int>> trueGeneCounts ; // True Matrix containing the Cell x Gene Matrix 
	std::vector<size_t> geneIdMap;
	
	// cell specific
	std::vector<std::string> cellNames ; // Vector of Cell Names 
	std::map<std::string, uint32_t> cellNamesMap ; // Cell Name -> Cell Id 
	std::map<std::string, uint32_t> cellNamesDupCount ; // Cell Name -> dedup count
	std::map<std::string, uint32_t> cellWhiteListMap ;
	std::unordered_map<std::string, uint32_t> cellNoisyMap ;
	std::unordered_map<std::string, uint32_t> cellDoubletMap ;

    std::unordered_map<uint32_t, uint32_t> cell2ClusterMap ;

	std::shared_ptr<spdlog::logger> consoleLog ; // Logger for outputting errors 
 
	std::map<uint32_t, uint32_t> alevin2refMap ; // Map Col 0f Input Matrix -> Gene ID from t2g tsv
	std::map<uint32_t, uint32_t> alevin2RevRefMap ; // Map Col 0f Input Matrix -> Gene ID from t2g tsv
	std::map<uint32_t, uint32_t> alevin2refTranscriptMap ; // Col of cell x transcript matrix -> Tr id from t2g tsv 

	uint32_t numCells{0} ;
	uint32_t numOfTranscripts{0} ;
	uint32_t numOfWhiteListedCells{0} ;
	uint32_t numOfNoisyCells{0} ;
	uint32_t numOfDoublets{0} ;

	// polyA related data structure 
	std::vector<double> fractionVector ; // Make the transcript id fractionVector[transcriptId] = fraction
	//std::unordered_map<uint32_t, std::vector<std::string>> transcriptIntronMap ; // Sequence of introns per gene
	std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, int>>> intronCountMap; // Map <Cellid, [<transcriptId , count>]>
	// Cached stuff
	// Vector of gene level probabiliies
	
	BFHClass* eqClassPtr{nullptr} ;
	GFAReader* dbgPtr{nullptr} ;
	// exon to transcript map
	// same for all clusters
	
	
	std::vector<uint32_t> rspdVec ;
	std::vector<std::unordered_map<size_t, uint32_t>> preCalculatedSegProb ;	
	std::vector<std::unordered_map<size_t, bool>> preSegOreMapVector ; // if this seg ever appears in rev 
	std::vector<std::unordered_map<uint32_t, std::unordered_map<size_t, int>>> cellSegCount ;
	//std::vector<std::unordered_map<size_t, std::vector<uint32_t>>> geneSpecificTrVector ; // Cache

	
	std::vector<std::unordered_map<size_t, std::vector<trInfo>>> geneSpecificTrVector ;// Ceched

};
#endif // MATRIX_PARSER
