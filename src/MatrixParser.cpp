#include <unordered_map>
#include <memory>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "MatrixParser.hpp"
#include "MinnowFS.hpp"
#include "BFHClass.hpp"
#include "GFAReader.hpp"

#include <unordered_set>
#include <numeric>

#include "zstr.hpp"


#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)

template <typename T>
void pupulateGeneCountFromTxt(
	std::string& countFile, 
	std::vector<std::vector<T>>& geneCount,
	uint32_t& numCells, 
	size_t& numOfGenes
){
	if(! util::fs::FileExists(countFile.c_str())){
		std::cerr << "quants_mat.csv file does not exist\n" ;
		std::exit(1) ; 
	}
	std::ifstream dataStream(countFile.c_str()) ;
	std::string line ; 

	size_t cellId{0} ;
	while(std::getline(dataStream, line)){
		cellId += 1 ;

		std::vector<std::string> cellGeneCountsString ; 
		util::split(line, cellGeneCountsString, ",") ;

		// should not have missing count 
		assert(cellGeneCountsString.size() == numOfGenes) ;

		std::transform(cellGeneCountsString.begin(), cellGeneCountsString.end(), geneCount[cellId].begin(), [](const std::string& val)
		{
			return std::stod(val);
		});
		if (cellId == numCells){
			break ;
		}
		_verbose("\rNumber of cells read genes : %lu", cellId);
	}

}

void loadTSVFile(
	std::string& tsvFile, 
	std::unordered_map<std::string, uint32_t>& geneMap,
	std::vector<double>& fractionVector
){
	if(! util::fs::FileExists(tsvFile.c_str())){
		std::cerr << "the tsv file does not exist\n" ;
		std::exit(1) ; 
	}

	fractionVector.resize(geneMap.size(), 0.0) ;
	std::ifstream dataStream(tsvFile.c_str()) ;
	std::string line ; 

	// fill the fraction vector
	while(std::getline(dataStream, line)){
		line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
		std::vector<std::string> tokens ;
		util::split(line, tokens, "\t") ;
		if(tokens.size() == 2){
			auto it = geneMap.find(tokens[0]) ;
			if(it != geneMap.end()){
				fractionVector[it->second] = std::stod(tokens[1]) ;
			}
		}
	}
	
}

template<typename T>
void populateGeneCountFromBinary(
	std::string& countFileBinary, 
	std::vector<std::vector<T>>& geneCount,
	uint32_t& numCells, 
	size_t& numOfGenes,
	std::unordered_map<uint32_t, uint32_t>& original2whitelistMap,
	std::unordered_map<uint32_t, uint32_t>& original2NoisyMap,
	uint32_t& numOfOriginalCells
)
{
	if(! util::fs::FileExists(countFileBinary.c_str())){
		std::cerr << "quants_mat.gz file does not exist\n" ;
		std::exit(1) ; 
	}
	bool predefinedCells = (original2whitelistMap.size() != 0) ;
	std::cerr << "Num cells ::: " << numCells << "\n" ;
	std::cerr << "Gene count dim " << geneCount.size() << "\t" <<  geneCount[geneCount.size()-1].size() << "\n" ;
	std::cerr << "noisyMap size " << original2NoisyMap.size() << std::endl  ;
	//std::cerr << "\n DEBUG ==> reading form binary with numCells: " << numCells <<  " numGenes: " << numOfGenes << "\n"  ;
	//std::cerr << "\n DEBUG ==> size of geneCount "<< geneCount.size() << " " << geneCount[geneCount.size() - 1].size() << "\n" ; 
	std::unique_ptr<std::istream> in =
    	std::unique_ptr<std::istream> (new zstr::ifstream(countFileBinary.c_str(), std::ios::in | std::ios::binary)) ;

	size_t elSize = sizeof(typename std::vector<double>::value_type);
	size_t cellCount{0} ;
	if(predefinedCells){
		for(size_t cellId = 0; cellId < numOfOriginalCells; ++cellId){
			std::vector<T> cell(numOfGenes) ;
			in->read(reinterpret_cast<char*>(cell.data()), elSize * numOfGenes) ;
			if (original2whitelistMap.find(cellId) != original2whitelistMap.end()){
				uint32_t whitelistCellId = original2whitelistMap[cellId] ;
				// std::cerr <<" cellid: " << cellId << "\tWhitelistCellId: " << whitelistCellId << std::endl ;
				if(whitelistCellId >= numCells)
					continue ;
				geneCount[whitelistCellId] = cell ;
				cellCount++ ;
				_verbose("\rNumber of cells read  : %lu", cellCount);
				if (cellCount == numCells){
					break ;
				}
			}
			else if(original2NoisyMap.find(cellId) != original2NoisyMap.end()){
				uint32_t noisyCellId = original2NoisyMap[cellId] ;
				geneCount[noisyCellId] = cell ;
				cellCount++ ;
				_verbose("\rNumber of cells read noisy  : %lu", cellCount);
				if (cellCount == numCells){
					break ;
				}
			}
			
		}
	}else{
		for(auto& cell : geneCount){
			cellCount += 1 ;
			//std::cerr << cellCount << "\n" ;
			_verbose("\rNumber of cells read  : %lu", cellCount);
			in->read(reinterpret_cast<char*>(cell.data()), elSize * numOfGenes) ;
			if (cellCount == numCells){
				break ;
			}
		}
		size_t numOfExpressedGenes{0} ;
		for(auto v : geneCount){
			for(auto c: v){
				if(c > 0){
					numOfExpressedGenes += 1 ;
				}
			}
		}

	}


	//std::cerr << "DEBUG: ==== right after reading binary expressedGeneCount " << numOfExpressedGenes << "\n" ;
	assert(cellCount == numCells) ;

}

/*
**	This module loads data from alevin folder, it reads 
**	the matrix and other related information from the same 
**	directory in order to simulate reads 
*/

template<typename T>
void DataMatrix<T>::loadAlevinData(	
	SimulateOptions& simOpts,
    Reference& refInfo
){

	// Load the values from simOpts
	auto& alevinDir = simOpts.matrixFile ;
	auto& binary = simOpts.binary ;
	auto& sampleCells = simOpts.sampleCells ;
	bool samplePolyA = simOpts.samplePolyA ;
	bool dupCounts = simOpts.dupCounts ; 
	bool generateNoisyCells = simOpts.generateNoisyCells ;
	auto cellClusterFile = simOpts.clusterFile ;
	auto outDir = simOpts.outDir ;
	bool createDoublet = (simOpts.numOfDoublets == 0) ? false : true ;

	bool useDBG = simOpts.useDBG ;
	auto gfaFile = simOpts.gfaFile ; 
	auto bfhFile = simOpts.bfhFile ;

	numOfDoublets = simOpts.numOfDoublets ;





	if (binary){
		consoleLog->info("The alevin matrix file is in binary mode") ;
	}else{
		std::cerr << "Should be in binary !!\n" ;
		std::exit(1) ; 
	}

	if(!simOpts.useWhiteList && generateNoisyCells){
		consoleLog->info("--generateNoisyCells needs to be invoked in conjunction with --useWhiteList") ;
	}	
	

	std::ifstream indata ;
	// load alevin related files 
	if(! util::fs::DirExists(alevinDir.c_str())){
		std::cerr << "Alevin directory does not exists\n" ;
		std::exit(1) ; 
	}

	// Here we model the exact same duplicated counts for each 
	// cell, this is proportional to the actually mapped the 
	// reference. This is explicitely obtained fron alevin run 
	 

	if (dupCounts){
		std::string dupCountFile = alevinDir + "/MappedUmi.txt" ;
		if(! util::fs::FileExists(dupCountFile.c_str())){
			std::cerr <<"filtered_cb_frequency.txt file does not exist\n" ;
			std::exit(1) ; 
		}
		std::ifstream dupCountStream(dupCountFile) ;
		std::string line ; 
		while(std::getline(dupCountStream, line)){
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			std::vector<std::string> tokens ; 
			util::split(line, tokens, "\t") ;
			if (tokens.size() == 2){
				cellNamesDupCount[tokens[0]] = std::stoul(tokens[1]) ;
			}
		}
		
	}

	// The map contains gid to tid map, where 
	// the gid is created first time here as, we 
	// get to know about them first here. The tids 
	// are taken from the fasta file. Other transcripts 
	// will be of no use as we don't know there 
	// sequence. 
	
	refInfo.updateGene2TxpMap(refInfo.gene2txpFile) ;
	auto& geneMap = refInfo.geneMap ;	
	consoleLog->info("Number of genes in the txp2gene file: {}", geneMap.size()) ;


	// This is an intermediate map that keeps a map 
 	// map between gene index and the gene names from 
	// the alevin produced file. This will be used to 
	// map the index of the matrix to gene index in 
	// the reference. 

	std::unordered_map<uint32_t, std::string> alevinGeneIndex2NameMap ;	

	// we might want to skip some genes
	std::unordered_set<uint32_t> skippedGenes ;
	uint32_t numOfSkippedGenes{0} ;
	size_t numOfOriginalGenes{0} ;


    consoleLog->info("Start parsing Alevin") ;
	consoleLog->info(
		"Parsing {}/quants_mat_cols.txt",alevinDir) ;

	
	{
		std::string geneListFile = alevinDir + "/quants_mat_cols.txt" ;
		if(! util::fs::FileExists(geneListFile.c_str())){
			consoleLog->error("quants_mat_cols.txt file does not exist") ;
			std::exit(1) ; 
		}
		std::ifstream geneNameStream(alevinDir + "/quants_mat_cols.txt") ;
		std::string line ;
		uint32_t geneIndex{0} ;
		while(std::getline(geneNameStream, line)){
			// strip new line
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			// search this gene name 
			auto it = geneMap.find(line) ;
			if(it != geneMap.end()){
				alevin2refMap[geneIndex] = it->second ;
				alevinGeneIndex2NameMap[geneIndex] = line ;
				geneIndex += 1 ;
			}else{
				//consoleLog->error("quant_mat_cols.txt contains gene {} that is not present in the reference ",line)  ;
				//std::exit(1) ;
				skippedGenes.insert(numOfOriginalGenes) ;
				numOfSkippedGenes++ ;
				
			}
			numOfOriginalGenes++ ;
		} 
	}

	std::string cellColFile = simOpts.outDir + "/alevin/quants_mat_cols.txt" ;
	std::ofstream cellColStream(cellColFile.c_str()) ;
	for(uint32_t i= 0; i < alevinGeneIndex2NameMap.size() ; ++i){
		cellColStream << alevinGeneIndex2NameMap[i] << "\n" ;
	}

	consoleLog->info("Original number of genes: {}\tNumber of genes skipped: {}", numOfOriginalGenes, numOfSkippedGenes) ;
	
	consoleLog->info("Number of genes in the alevin produced files: {}",alevin2refMap.size()) ;
	
	auto& gene2transcriptMap = refInfo.gene2transcriptMap ;
	
	// alevin2refTranscriptMap, a map from columns of the 
	// cell x transcrip count matrix to be fromed to the 
	// transcript ids of the reference. This map is *very*
	// important since we need the reference id to get the
	// real sequences. 

	std::map<uint32_t, uint32_t> alevinReverseMap ;
	{
		// create a transcript map
		uint32_t trId{0} ;
		for(auto it : alevin2refMap){
			auto gIt = gene2transcriptMap.find(it.second) ;
			if(gIt != gene2transcriptMap.end()){
				auto& trVec = gIt->second ;
				//std::cerr << "# of transcripts in this gene "<< gIt->first 
				//	      <<  ":\t" << trVec.size() << "\n" ;	
				for(auto tr : trVec){
					alevin2refTranscriptMap[trId] = tr ;
					alevinReverseMap[tr] = trId ;
					trId += 1 ;
				}
			}else{
				consoleLog->error("This should not happen: gene {} not found",gIt->first) ;
			}
		}
		 
	}
    
    consoleLog->info("Parsing {}/quants_mat_rows.txt",alevinDir) ;
	
	// all cell names irrespective of whitelist
	std::map<std::string, uint32_t> allCellListMap ; // map cell-name -> id
	std::vector<std::string> allCellNames ; // {cell names}
	
	
	// track a cell by barcode
	std::string cellToTrack  = ""; //"GTATTCTGTAGTGAAT"; // "TGCGGGTCAGACGCAA" ;
	//uint32_t cellIdToTrack{0} ;


	// read cell names/barcodes that alevin produce
	{
		std::string cellListFile = alevinDir + "/quants_mat_rows.txt" ;
		if(! util::fs::FileExists(cellListFile.c_str())){
			consoleLog->error("quants_mat_rows.txt file does not exist") ;
			std::exit(1) ; 
		}
		std::ifstream cellNameStream(alevinDir + "/quants_mat_rows.txt") ;
		std::string line ;
		while(std::getline(cellNameStream, line)){
			// strip new line
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			allCellNames.push_back(line) ;
			allCellListMap[line] = allCellNames.size() - 1 ;
		} 
	}


	uint32_t numOfOriginalCells = allCellNames.size() ;

	consoleLog->info("Number of cells in the alevin produced files: {}",numOfOriginalCells) ;
	
	
	// Read the whitelist to know which cells are going to be finally part of the matrix.
	// If the whitelist is not present then treat the rows as whitelist 


	std::string cellListFile = alevinDir + "/whitelist.txt" ;
	
	// When noisy cells are present and considered, then  
	// allCells are the combination of whitelist cells and 
	// noisy cells. These two maps would keep track of those 
	// two groups and map original id of the cell list to 
	// the whitelist and noisylist.

	// FIXME: If we sample from white-list, it might not be contiguous, 
	// it is not contiguous create a second order mapping. in that case
	// we should first
	  
	 
	std::unordered_map<uint32_t, uint32_t> original2whitelistMap ;  
	std::unordered_map<uint32_t, uint32_t> original2NoisyMap ;
	
	// Minnow would also keep a list of doublet cells, which will be
	// created later 
	std::unordered_map<uint32_t, uint32_t> original2DoubletMap ;
	std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> doubletContainer ;
	


	// Simulate noisy cell if asked by user. 
	
	// If whitelist by alevin is not present we will treat all cells as whitelist
	// cells. 
	// TODO: Given an in-built noise detector we can use to threshold 
	// and call it noise that would be great. Not now.   


	if((! util::fs::FileExists(cellListFile.c_str())) or !(simOpts.useWhiteList)){
		std::cerr << "whitelist.txt file does not exist/ or will NOT be used\n" 
				  << "we need to assume that the rows are \n" 
				  << "the whitelisted barcodes \n" ;
		
		// Copy the vector and map for now, and try something better
		// later

		cellNames = allCellNames ;
		cellWhiteListMap = allCellListMap ;

	}else{
		
		consoleLog->info("Parsing {}/whitelist.txt",alevinDir) ;
		
		{

			std::ifstream cellNameStream(alevinDir + "/whitelist.txt") ;
			std::string line ;
			while(std::getline(cellNameStream, line)){
				// strip new line
				line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
				cellNames.push_back(line) ;
				cellWhiteListMap[line] = cellNames.size() - 1 ;
				// cellsToBeRead.insert(allCellListMap[line]) ;
				//original2whitelistMap[allCellListMap[line]] = cellWhiteListMap[line] ;
			} 
		}
		// fill up the original2whitelistMap by going 
		// over the cellWhiteListMap
		for(size_t i = 0 ; i < cellNames.size(); ++i){
			original2whitelistMap[allCellListMap[cellNames[i]]] = i ; 
				
		} 

		

		// NOTE: Noisy cells if neccessary will be appended 
		// in the same vector of cell names. This makes 
		// cellWhiteListMap deprecated for the case where 
		// whitelist is provided. as original2whitelistMap and 
		// allCellListMap together provide the same information 
		
		// TODO: remove cellWhiteListMap in future release 

		consoleLog->info("Number of cells in whitelist file: {}", cellNames.size()) ;
		if(generateNoisyCells){
			consoleLog->info("Additionally reads from noisy cells will be generated too,"
							 "keeping track of noisy cells\n");

			if(numOfNoisyCells == 0){
				numOfNoisyCells = allCellNames.size() - cellNames.size() ;
			}
			uint32_t noisyCellId = cellNames.size() ;
			for(auto it : allCellListMap){
				if(cellWhiteListMap.find(it.first) == cellWhiteListMap.end()){
					cellNames.push_back(it.first) ;
					cellNoisyMap[it.first] = cellNames.size() - 1 ;
					original2NoisyMap[it.second] =  cellNoisyMap[it.first];
					noisyCellId++ ;
					// TODO: insert a condition to break

				}
			}			
			consoleLog->info("Additionally generate data from {} noisy experiments", noisyCellId) ;	
			
		}
	}

		
	if ((sampleCells > 0) and (sampleCells < cellNames.size())){
		std::cerr << "\n We will sample "<< sampleCells  << " from " << cellNames.size() << "\n" ;
	}else{
		sampleCells = cellNames.size() ;
	}

	numOfWhiteListedCells = cellWhiteListMap.size() ;
	numCells = sampleCells ;
	numOfTranscripts = alevin2refTranscriptMap.size() ;
	size_t numOfGenes = alevin2refMap.size() ;

	// TODO: add ablity to create 
	// more doublets 	
	//if(createDoublet){
	//	numOfDoublets =  1;
	//}

	// We create a temporary gene count matrix which can serve the purpose of before 
	// writing the counts to the data file 
	std::vector<std::vector<T>> originalGeneCountMatrix ;

	originalGeneCountMatrix.resize(sampleCells, std::vector<T>(numOfOriginalGenes)) ;

	data.resize(sampleCells + numOfDoublets, std::vector<T>(numOfTranscripts)) ;
	geneCounts.resize(sampleCells + numOfDoublets, std::vector<T>(numOfGenes)) ;
	trueGeneCounts.resize(sampleCells + numOfDoublets, std::vector<int>(numOfGenes)) ;

	// this matrix might have more genes 
	for(auto& v : originalGeneCountMatrix)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());
	
	for(auto& v : data)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());
	
	for(auto& v : geneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : trueGeneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	std::string countFile = alevinDir + "/quants_mat.csv" ;
	std::string countFileBinary = alevinDir + "/quants_mat.gz" ;


	std::unordered_map<uint32_t, uint32_t> gene2LastTrMap ;
	{
		// Make a map from gene id to last last tid
		for(auto& gmpair : geneMap){
			auto trIds = refInfo.gene2transcriptMap[gmpair.second] ;
			auto lastTid = trIds[trIds.size() - 1] ;
			auto alevinLastTid = alevinReverseMap[lastTid] ;
			gene2LastTrMap[gmpair.second] = alevinLastTid ; 
		}

	}

	// Decide whether to sample from the introns or not 
	// depending on that use a flag for that cell and gene id 
	// store a flag
	std::unordered_set<uint32_t> geneIntronMap ;

	if(samplePolyA){
		consoleLog->info("Sampling from polyA sites in between gene") ;
		auto& tsvFile = simOpts.polyAsiteFractionFile ; // fraction file 
		auto& polyAFastaFile = simOpts.polyAsiteFile ; // fasta file

		// load polyA fraction 
		loadTSVFile(
			tsvFile, 
			geneMap,
			fractionVector
		) ;	
		// load the fasta file with the sequences 
		refInfo.updatePolyASequence(
			polyAFastaFile,
			geneMap, 
			gene2LastTrMap,
			geneIntronMap	
		) ;
	} 

	bool useClusters{false} ;


	// TODO: 
	if(useDBG){
		// This will read the dbg now 
		if(!util::fs::FileExists(gfaFile.c_str())){
			std::cerr << gfaFile << " should exist EXITING !!"; 
			std::exit(2) ;
		}

		auto rspdFile = simOpts.rspdFile ;
		if(rspdFile != ""){
			if(util::fs::FileExists(rspdFile.c_str())){
				rspdVec.resize(MAX_FRAGLENGTH) ;
				std::ifstream rspdStream(rspdFile.c_str()) ;
				std::string line ;
				while(std::getline(rspdStream, line)){
					std::vector<std::string> tokens ;
					util::split(line, tokens, ",") ;
					if(tokens.size() != 2)
						continue ;
					uint32_t pos = std::stoul(tokens[0]) ;
					uint32_t prob = std::stoul(tokens[1]) ;
					if(pos < rspdVec.size()){
						rspdVec[pos] = prob ;
					}
				}
			}else{
				std::cerr << "RSPD file is not empty and doesn't exist, going with truncated sampling\n" ;
			}
		} 


		dbgPtr = new GFAReader(gfaFile) ;
		dbgPtr->parseFile(refInfo) ;

		// stand alone call to bfh
		// ignore other eqclass stuff 
		std::string eqFileDir = "/mnt/scratch1/hirak/minnow/metadata/hg/" ;

		eqClassPtr = new BFHClass(eqFileDir) ;
		eqClassPtr->loadBFH(
			bfhFile, 
			cellClusterFile, 
			refInfo, 
			cellWhiteListMap, 
			false, 
			cellNoisyMap
		) ;

		preCalculatedSegProb.resize(numOfGenes) ; // per gene
		preSegOreMapVector.resize(numOfGenes) ;
		geneSpecificTrVector.resize(numOfGenes) ;

		// fill segment based probability
		// gene to tr to prob mapping 

		//uint32_t trackGid{54819} ;
		for(uint32_t i = 0; i < numOfGenes; ++i){
			auto it = alevin2refMap.find(i) ;
			if(it != alevin2refMap.end()){
				auto originalGeneId = it->second ;
				auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

				std::unordered_map<size_t, uint32_t> localGeneProb ;
				std::unordered_map<size_t, bool> localSegOreMap ;
				std::unordered_map<size_t, std::vector<trInfo>> localTrVector ;

				for(auto tid : transcriptIds){
					if(dbgPtr->trSegmentMap.find(tid) != dbgPtr->trSegmentMap.end()){
						auto segVec = dbgPtr->trSegmentMap[tid] ;
						bool everRC{false} ;
						for(auto seg : segVec){
							auto segCount = dbgPtr->eqClassMapCounts[seg] ;
							auto bfhCount = eqClassPtr->getGeneLevelProbCount(
								originalGeneId,
								segCount 
							) ;
					
							if(bfhCount != 0){
								// old one
								// localTrVector[seg].emplace_back(tid) ;

								auto tcInfoVec = dbgPtr->eqClassMap[seg][tid] ;
								for(auto tInfo : tcInfoVec){
									if(refInfo.transcripts[tid].RefLength - tInfo.eposInContig <= MAX_FRAGLENGTH){
										localGeneProb[seg] = bfhCount ;
										localTrVector[seg].emplace_back(
											tid,
											tInfo.sposInContig,
											tInfo.eposInContig
										) ;
									}

									if(!tInfo.ore && !everRC){
										everRC = true ;
									}
								}
								if(everRC){
									localSegOreMap[seg] = false ;
								}else{
									localSegOreMap[seg] = true ;
								}
							}
						}
					}
				}

				preCalculatedSegProb[i] = localGeneProb ;
				preSegOreMapVector[i] = localSegOreMap ;
				geneSpecificTrVector[i] = localTrVector ;
				
			}
		} 

		std::cerr << " After loading bfh prob size " << preCalculatedSegProb.size() << "\n" ;
		cellSegCount.resize(numCells + numOfDoublets) ;

		std::ofstream probStream("/mnt/scratch1/hirak/minnow/geneLevelProb.txt") ;
		probStream << eqClassPtr->geneCountHistogram.size() << "\n" ;
		for(auto it : eqClassPtr->geneCountHistogram){
			probStream << it.second.size() << "\n" ;
			for(auto it2 : it.second){
				probStream << it2.first << "\t" << it2.second << "\n" ;
			}
		}

	}
	

	{

		size_t nonWhiteLisBarcodesSkipped{0} ;

		// The alevin data can be binary 
		// or can be non-binary, based on which
		// we decide it to read and fill the geneCount
		// matrix if useWhitelist is on we have to read
		// the entire matrix. 
		// Given N whitelistes and M noisy barcodes, if
		// asked minnow will generate (N+M) x numOfTranscripts
		// matrix in this order. 

		// if given in txt format.
		// the gene count matrix produced
		// by Alevin has the structure of 
		// cell to gene counts with a
		// comma as delimeter 
		// each line also ends with comma

		// we would simultaneously create the
		// cell x transcript matrix, it will have more 
		// columns presumably

		// Using a two level selection process for 
		// converting the the gene counts to transcript 
		// level counts. It is a multinomial distribution 
		// to begin with where gene counts are of type double	


		if(!binary){
			pupulateGeneCountFromTxt(countFile, geneCounts, numCells, numOfGenes) ;
		}else{
			std::cout << "\n" ;
			consoleLog->info("Loading the binary Matrix") ;
			populateGeneCountFromBinary(countFileBinary, originalGeneCountMatrix, numCells, numOfOriginalGenes, original2whitelistMap, original2NoisyMap, numOfOriginalCells) ;
		}

		consoleLog->info("Loaded the Matrix") ;

		// Make a truncated geneCounts file
		//T skippedCount{0} ;
		if(numOfSkippedGenes > 0){
			for(size_t cell_id = 0 ; cell_id < originalGeneCountMatrix.size(); ++cell_id){
				size_t geneCountsGeneIdx{0} ;
				for(size_t gene_id = 0 ; gene_id < numOfOriginalGenes ; ++ gene_id){
					if(skippedGenes.find(gene_id) == skippedGenes.end()){
						geneCounts[cell_id][geneCountsGeneIdx] = originalGeneCountMatrix[cell_id][gene_id] ;
						geneCountsGeneIdx++ ;
					}
				}
			}	
		}else{
			geneCounts = originalGeneCountMatrix ;
		}
		
		consoleLog->info("Truncated the matrix ") ;

		// CREATE Doublets 	
		if(createDoublet){
			// Treat doublets as normal cells and put them in the list of all cells 
			auto CB10XList = util::generate10XCBList(numOfDoublets) ;
			for(size_t i = 0; i < numOfDoublets; ++i){
				allCellNames.push_back(CB10XList[i]) ;
				allCellListMap[CB10XList[i]] = allCellNames.size() - 1 ;

				//std::cerr << "numOfWhiteListedCells "<< numOfWhiteListedCells << "\n" ;
				// generate doublet
				auto cell_id_1 = util::generateRandomNumber(0, static_cast<int>(numOfWhiteListedCells-1)) ;
				auto cell_id_2 = util::generateRandomNumber(0, static_cast<int>(numOfWhiteListedCells-1)) ;

				std::vector<T> doubletCount(numOfGenes) ;
				for(size_t j = 0; j < numOfGenes ; ++j){
					geneCounts[sampleCells+i][j] = 
						static_cast<T>((geneCounts[cell_id_1][j] + geneCounts[cell_id_2][j])/2) ;
				}

					
				// increase the cell numbers and update the maps 
				cellNames.push_back(CB10XList[i]);
				cellDoubletMap[CB10XList[i]] = cellNames.size()-1 ; 

				doubletContainer[cellDoubletMap[CB10XList[i]]] = {cell_id_1, cell_id_2} ;

			}

		}


		// By this stage the geneCount matrix should be 
		// filled up and ready to be converted to the transcript level 
		// counts 

		size_t cellId{0} ;
		size_t numOfExpressedGenesInput{0}; 
		//size_t numOfExpressedGenesOutput{0} ;
			

		for(auto& cellGeneCounts : geneCounts){
			// check if this cell Id is in cellWhiteListMap
			// std::cerr << "\nreading cell line\n" ;
			int droppedGeneExpression{0} ;
			uint32_t thisCellValidateGeneCount{0} ;
			//uint32_t thisCellValidateExonCount{0} ;
			uint32_t thisCellValidateTrCount{0} ;

			auto barcode = cellNames[cellId] ;
			uint32_t backupActualCellId ;

			// This step is unncessary when 
			// whitelist is present, well sort of 
			auto itw = cellWhiteListMap.find(barcode) ;
			auto itn = cellNoisyMap.find(barcode) ;
			auto itd = cellDoubletMap.find(barcode) ;
			if (
				(itw != cellWhiteListMap.end()) ||
				(itn != cellNoisyMap.end()) ||
				(itd != cellDoubletMap.end())
			){
				
				bool trackTheCell{false} ;
				bool isNoisyCell{false} ;
				bool isDoublet{false} ;

				uint32_t actualCellId ;				
				if(itn != cellNoisyMap.end()){
					isNoisyCell = true ;
					actualCellId = itn->second ;
				}else if(itw != cellWhiteListMap.end()){
					actualCellId = itw->second ;
				}else{

					actualCellId = itd->second ;
				}

				if(itw->first == cellToTrack){
					trackTheCell = true ;
				} 				
				//auto actualCellId = it->second ;
				backupActualCellId = actualCellId ;


				// If clusterMap is present then
				uint32_t cluster_id{0} ;
				if(useClusters && !isNoisyCell){
					// if it's doublet use one of the clusters
					if(isDoublet){
						auto parent_id = doubletContainer[actualCellId].first ;
						cluster_id = eqClassPtr->cell2ClusterMap[parent_id] ;
					}else{
						cluster_id = eqClassPtr->cell2ClusterMap[actualCellId] ;
					}
				
					if(trackTheCell){ std::cerr << "the cluster id for " << barcode << " is " << cluster_id << "\n" ;}
				}
				
				int validateGeneCount{0} ;
				//int validateExonCount{0} ;
				//int validateTrCount{0} ;


				// treat the vector as a distribution
				std::vector<int> cellGeneCountSampled(numOfGenes,0) ;

				//double totGeneCountDouble = std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0) ;

				int totGeneCount = static_cast<int>(std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0)) ;

				//std::cerr << "DEBUG: Double and Int " << totGeneCountDouble << "\t" << totGeneCount << "\n" ;

				int toSubTract{0} ; 
				for(uint32_t i = 0; i < cellGeneCounts.size() ; ++i){
					if(cellGeneCounts[i] >= 1) {
						--cellGeneCounts[i] ;
						++cellGeneCountSampled[i] ;
						++toSubTract ;
					}
				}
				
				totGeneCount -= toSubTract ;
				// Sample form a multinomial dist with prob dist cellGeneCounts
				

				std::random_device rdg;
    			std::mt19937 geng(rdg());
				std::discrete_distribution<> dg(cellGeneCounts.begin(), cellGeneCounts.end()) ;
				for(int i = 0; i < totGeneCount ; ++i){
					++cellGeneCountSampled[dg(geng)] ; 
				}
				//std::cerr << "\nDEBUG --> After multinomial sampling\n" ;
				//std::cerr << "\nDEBUG --> cellGeneCountSampled size "<< cellGeneCountSampled.size() <<"\n" ;
				trueGeneCounts[cellId].assign(cellGeneCountSampled.begin(), cellGeneCountSampled.end()) ;	
				size_t numOfExpressedGenes2{0} ;
				for(auto v : cellGeneCountSampled){
						if(v > 0){
							numOfExpressedGenes2 += 1 ;
					}
				}
				int numOfDroppedGenes{0} ;
				//totGeneCountDouble = std::accumulate(cellGeneCounts.begin(), cellGeneCounts.end(), 0.0) ;
				//std::cerr << "DEBUG: After multinomial sampling number of expressed genes " << numOfExpressedGenes2 
				//          << "  Before number of expressed genes " << toSubTract <<  "\n" ;
				for(size_t i = 0 ; i < cellGeneCountSampled.size(); ++i){
					// geneId 
					std::string geneName = alevinGeneIndex2NameMap[i] ;	

					int geneCount =  cellGeneCountSampled[i] ;
					//int geneLevelTrCount{0} ;
					//int geneLevelExCount{0} ;

					validateGeneCount += geneCount ;
					if (geneCount > 0){

						bool dropThisGene{false} ;
						numOfExpressedGenesInput++ ;

						auto it = alevin2refMap.find(i) ;
						

						if (it != alevin2refMap.end()){
							auto originalGeneId = it->second ;
							// fill geneCount
							// Gene count is greater than 0 and gene count 
							// Check if this gene is marked for polyA
							int intronicCount{0} ;

							if(samplePolyA){
								auto it2 = geneIntronMap.find(originalGeneId) ;
								if(it2 != geneIntronMap.end()){
									if(fractionVector[originalGeneId]){
										intronicCount = static_cast<int>(geneCount * fractionVector[originalGeneId]) ;
										geneCount = geneCount - intronicCount ;
									}		
								} 
							}

							auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

							// multinomial distribution with probabilities 
							// from the weibull distribution above
							std::vector<double> probVec(transcriptIds.size()) ;
							std::unordered_map<uint32_t, size_t> transcriptIdMap ;
							
							for(size_t ttid = 0 ; ttid < transcriptIds.size() ; ++ttid){
								transcriptIdMap[transcriptIds[ttid]] = ttid ; 
							}

							// DBG based allocation
							if(useDBG){
								dropThisGene = false ;
								auto segCountHist = preCalculatedSegProb[i] ;
								
								if(segCountHist.size() == 0){
									droppedGeneExpression++ ;
									dropThisGene = true ;
									geneCount = 0 ;
								}

								if(!dropThisGene){
									std::vector<size_t> segIndex(segCountHist.size()) ;
									std::vector<uint32_t> segProbVec(segCountHist.size()) ;
									std::vector<int> segCounts(segCountHist.size(),0) ;
									
									size_t ind{0} ;
									for(auto it : segCountHist){
										segIndex[ind] = it.first ;
										segProbVec[ind] = it.second ;
										ind++ ;
									}

									std::random_device rdt4;
									std::mt19937 gent4(rdt4());

									std::unordered_map<size_t, int> segCountMap ;
									
									std::discrete_distribution<> dmt(segProbVec.begin(), segProbVec.end()) ;
									for(int j = 0 ; j < geneCount ; ++j){
										++segCounts[dmt(gent4)] ;
									}
									for(size_t j = 0 ; j < segCounts.size() ; ++j){
										segCountMap[segIndex[j]] = segCounts[j] ;

									}

									cellSegCount[actualCellId][i] = segCountMap ;

								}
							}
							else{

								// These are the learned distribution weibull paeameters from the paper 
								// https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005761
								// Assign probability to individual transcripts
								
								// do it only if gene count is more than 0 		
								if(trackTheCell) {std::cerr << "Number of transcripts " << transcriptIds.size() << "\n";}

								std::random_device rdt;
								std::mt19937 gent(rdt());

								std::weibull_distribution<> dwt(0.44,0.6);

								std::random_device rd2;
								std::mt19937 g2(rd2()) ;
								std::vector<int> tIndex(transcriptIds.size()) ;
								std::iota(tIndex.begin(), tIndex.end(), 0);
								std::shuffle(tIndex.begin(), tIndex.end(), g2) ;

								for(size_t j = 0; j < transcriptIds.size() ; ++j ){
									probVec[tIndex[j]] = dwt(gent) ;
								}
								

							}


							// If dbg is used we don't work with transcripts 
							trueGeneCounts[cellId][i] = geneCount ;
							
							if(useDBG){
								continue ;
							}
								

							auto sum_prob = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
							if(sum_prob == 0){
								dropThisGene = true ;
							}

							if (dropThisGene){
								//cellExCount[actualCellId][i] = {} ;
								numOfDroppedGenes++ ;
								continue ;
							}
							
							std::random_device rdt3;
							std::mt19937 gent3(rdt3());


							std::discrete_distribution<> dmt(probVec.begin(), probVec.end()) ;
							// saple from this distribution geneCount times
							std::vector<int> transcriptCounts(transcriptIds.size(), 0) ; 
							for(int j=0 ; j < geneCount ; ++j){
								++transcriptCounts[dmt(gent3)] ;
							}

							//auto beforeFuckingUp = std::accumulate(
							//	transcriptCounts.begin(),
							//	transcriptCounts.end(),
							//	0
							//);
							// total gene count 
							T totCount{0} ;
							uint32_t lastAlevinTid{0} ;
							for(size_t j = 0; j < transcriptIds.size(); ++j){
								auto tid = transcriptIds[j] ;
								auto alevinTid = alevinReverseMap[tid] ;
								data[actualCellId][alevinTid] = transcriptCounts[j] ;
								totCount += data[actualCellId][alevinTid] ;
								lastAlevinTid = alevinTid ;
							}

							thisCellValidateTrCount += totCount ;
							thisCellValidateGeneCount += geneCount ;

							// Add intronic count if there is a 
							// non-zero quantity to add 
							if (intronicCount > 0){
								auto it3 = intronCountMap.find(actualCellId) ;
								if (it3 != intronCountMap.end()){
									intronCountMap[actualCellId].push_back(std::make_pair(lastAlevinTid, intronicCount)) ;
								}else{
									intronCountMap[actualCellId] = {{lastAlevinTid, intronicCount}} ;
								}
							}

						}else{
							consoleLog->error("This gene id: {} not found\n", i) ;
						}
					}

				}

				auto iVec = intronCountMap[backupActualCellId] ;
				int intronicSum = 0 ;
				for(auto i : iVec){
					intronicSum += i.second ;
				}
				//std::cerr << barcode << "\t" << backupActualCellId << "\t" << cellId <<  "\t" <<  validateGeneCount << "\t" << intronicSum  << "\n" ;

			}else{
				//std::cerr << "barcode not found " << barcode << "\n" ;
				nonWhiteLisBarcodesSkipped++ ;
			}

			if(useDBG){
				//consoleLog->info("DBG::: total Gene Expression skipped {}", droppedGeneExpression) ;
			}


			cellId += 1 ;
			_verbose("\rNumber of cells processed : %lu", cellId);

		}

		//std::cerr << "DEBUG: ==== Before dumping input expressedGeneCount " << numOfExpressedGenesInput << "\n" ;
		std::cout << "\n" ;
		consoleLog->info("{} cell barcodes are skipped as they are not in whitelist", nonWhiteLisBarcodesSkipped) ;
	}

	// FIXME: this should not be ad-hoc
	numCells = geneCounts.size() ;


	consoleLog->info("The transcript matrix is constructed, with dimension {} x {} \n" 
	                 "\t\t\t\t\tfrom gene count with dimention {} x {}, geneCounts.size(): {}",
					 data.size(), data[0].size(), cellNames.size(), alevin2refMap.size(), geneCounts.size()) ;

}


// NOTE: This matrix is numOfGenes x numCells 

template <typename T>
void readSplatterMatrix(
	std::string& countFile, 
	std::vector<std::vector<T>>& geneCount,
	size_t numCells, 
	size_t numOfGenes
){
	if(! util::fs::FileExists(countFile.c_str())){
		std::cerr << "quants_mat.csv file does not exist\n" ;
		std::exit(1) ; 
	}
	std::ifstream dataStream(countFile.c_str()) ;
	std::string line ; 

	size_t geneId{0} ;
	while(std::getline(dataStream, line)){

		std::vector<std::string> cellCountsString ; 
		util::split(line, cellCountsString, ",") ;

		for(size_t cellId = 0 ; cellId < numCells ; ++cellId){
			geneCount[cellId][geneId] = std::stoi(cellCountsString[cellId]) ; 
		}
		geneId += 1 ;
		_verbose("\rIn Splatter: Number of genes processed : %lu", geneId);
		if(geneId == numOfGenes){
			break ;
		}
	}

}


void readUniqueness(
     std::vector<std::pair<std::string, double>>& uniquenessInfo,
	 std::unordered_map<std::string, uint32_t>& geneMap, 
     bool useReverse
  
){

	std::string uniquenessListFile{"../data/hg_stranded_gene_uniqueness.txt"} ;
	std::ifstream uniqueFileStream(uniquenessListFile) ;

	std::string line ;
	// throw the fist line 
	std::getline(uniqueFileStream, line) ;

	while(std::getline(uniqueFileStream, line)){
		// strip new line
		line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
		std::vector<std::string> tokens ;
		util::split(line, tokens, "\t") ;
		auto geneName = tokens[0] ;
		double score{0.0} ;
		if(std::stod(tokens[1]) != 0){
			score = std::stod(tokens[1]) / std::stod(tokens[2]) ;
		}
		if(geneMap.find(geneName) != geneMap.end()){
			uniquenessInfo.push_back({tokens[0], score}) ;
		}
	}

	using ptype = std::pair<std::string, double> ;

  // This should keep uniq genes at top 
  if(useReverse){
		std::sort(uniquenessInfo.begin(), uniquenessInfo.end(), [](ptype &left, ptype &right) {
			return left.second > right.second;
		});
	}else{
    // This should keep least uniq genes at top
		std::sort(uniquenessInfo.begin(), uniquenessInfo.end(), [](ptype &left, ptype &right) {
			return left.second < right.second;
		});
	}
}

void makeUniquenessBins(
	std::vector<std::unordered_set<std::string>>& uniqBins,
	std::unordered_map<std::string, uint32_t>& geneMap
){


	std::string uniquenessListFile{"../data/hg_stranded_gene_uniqueness.txt"} ;
	std::ifstream uniqueFileStream(uniquenessListFile) ;
	std::vector<std::pair<std::string, double>> uniquenessInfo ;	
	uniqBins.resize(10) ;
	std::vector<double> scoreBins(11) ;

	int ind{0} ;
	for(double i = 0 ; i <= 1.0 ; i += 0.1){
		scoreBins[ind] = i ;
		ind++ ;
		//std::cerr << "score ::: " << ind << "\t" << scoreBins[ind-1] << "\n" ;
	}

	std::string line ;
	// throw the fist line 
	std::getline(uniqueFileStream, line) ;


	while(std::getline(uniqueFileStream, line)){
		// strip new line
		line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
		std::vector<std::string> tokens ;
		util::split(line, tokens, "\t") ;
		auto geneName = tokens[0] ;
		double score{0.0} ;
		if(std::stod(tokens[1]) != 0){
			score = std::stod(tokens[1]) / std::stod(tokens[2]) ;
		}
		bool gotin{false} ;
		if(geneMap.find(geneName) != geneMap.end()){
			uniquenessInfo.push_back({tokens[0], score}) ;
			for(int i = 0; i < 10 ; ++i){
				//if(lookin) std::cout << scoreBins[i] << "\t" << scoreBins[i+1] << "\t" << score << "\n";
				if((score >= scoreBins[i]) && (score <= scoreBins[i+1])){
					gotin = true ;
					uniqBins[i].insert(tokens[0]) ;
					break ;
				}
				if(score == 1){
					uniqBins[9].insert(tokens[0]) ;	
					gotin = true ;
				}
			}

			if(!gotin){
				std::cout <<"Score::: " << score << "\n" ;
				std::exit(1) ;
			}
		
		}
	}

	for(size_t i = 0 ; i < uniqBins.size() ; ++i){
		//std::cerr << scoreBins[i] << "\t" << uniqBins[i].size() << "\n" ;
	}


}


template<typename T>
void DataMatrix<T>::loadSplatterData(
	SimulateOptions& simOpts,
    Reference& refInfo
){

	// Load data from the simOpts 
	auto splatterDir = simOpts.matrixFile ;
	size_t sampleCells = simOpts.sampleCells ;
	std::string outDir = simOpts.outDir ;
	std::string bfhFile = simOpts.bfhFile ;
	std::string gfaFile = simOpts.gfaFile ;
	//bool useEqClass = simOpts.useEqClass ;
	bool useWeibull = simOpts.useWeibull ;
	bool useDBG = simOpts.useDBG ;


	std::string geneListFile = splatterDir + "/quants_mat_rows.txt" ;
	std::string cellListFile = splatterDir + "/quants_mat_cols.txt" ;

	// this is gene to cell matrix 
	std::string countFile = splatterDir + "/quants_mat.csv" ;

	if(! util::fs::FileExists(geneListFile.c_str())){
		std::cerr << geneListFile <<" file does not exist\n" ;
		std::exit(1) ; 
	} 
	if(! util::fs::FileExists(cellListFile.c_str())){
		std::cerr << cellListFile <<" file does not exist\n" ;
		std::exit(1) ; 
	} 


	if(! util::fs::FileExists(countFile.c_str())){
		std::cerr << countFile <<" file does not exist\n" ;
		std::exit(1) ; 
	} 
	
	// Splatter does not specify gene names 
	// We would randomly sample genes from from 
	// the reference provided. Everything will be treated 
	// as whitelisted cells. 

	// update txp to gene
	refInfo.updateGene2TxpMap(refInfo.gene2txpFile) ;

	auto& geneMap = refInfo.geneMap ;
	


	consoleLog->info("Number of genes in the txp2gene file: {}", geneMap.size()) ;
	
	consoleLog->info(
		"Parsing {}/quants_mat_rows.txt",splatterDir) ;	
	
	std::map<std::string, uint32_t> allCellListMap ;
	std::vector<std::string> allCellNames ;
	std::string cellToTrack  = ""; 
	//uint32_t cellIdToTrack{0} ;
	{
		std::ifstream cellNameStream(cellListFile) ;
		std::string line ;
		while(std::getline(cellNameStream, line)){
			// strip new line
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			allCellNames.push_back(line) ;
			allCellListMap[line] = allCellNames.size() - 1 ;
		} 
	}
	// In splatter there is no whitelist but there can be clusters 
	consoleLog->info("{} cells are present ", allCellNames.size()) ;
	consoleLog->info("Start parsing Splatter output") ;
	consoleLog->info(
		"Parsing {}/quants_mat_cols.txt",splatterDir) ;
	
	std::unordered_map<std::string, std::string> splatterGeneMap ; //splatter gene name -> real gene name
	std::vector<std::string> splatterGeneNames ;
	std::unordered_map<uint32_t, std::string> aleviGeneIndex2NameMap ;

	{
		std::ifstream geneNameStream(geneListFile) ;
		std::string line ;
		while(std::getline(geneNameStream, line)){
			// strip new line
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
			splatterGeneNames.push_back(line) ;
		} 
	}

	size_t numOfGenes = splatterGeneNames.size() ;


	// Let's assume no whitelist 
	{
		cellNames = allCellNames ;
		cellWhiteListMap = allCellListMap ;
	}

	if ((sampleCells > 0) and (sampleCells < cellNames.size())){
		std::cerr << "\n We will sample "<< sampleCells  << " from " << cellNames.size() << "\n" ;
	}else{
		sampleCells = cellNames.size() ;
	}

	numCells = sampleCells ;

	geneCounts.resize(sampleCells, std::vector<T>(numOfGenes)) ;

	trueGeneCounts.resize(sampleCells, std::vector<int>(numOfGenes)) ;

	
	for(auto& v : geneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());

	for(auto& v : trueGeneCounts)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());
	
	readSplatterMatrix(
		countFile,
		geneCounts, 
		allCellNames.size(),
		splatterGeneNames.size() 
	) ;


	std::cout << "\n" ;
	consoleLog->info("Splatter matrix is read, with dimension {} x {}",geneCounts.size(), geneCounts[0].size()) ;
	std::unordered_map<uint32_t, std::string> alevinGeneIndex2NameMap ;
	{
		// go over gene to transcript map, include a 
		// gene if contains a transcript that is present 
		// in the reference fasta file

		bool testUniqueness{simOpts.testUniqness} ;
    bool revUniqness{simOpts.reverseUniqness} ;

		if(testUniqueness){
			// read uniqueness file 

			std::vector<std::pair<std::string, double>> uniquenessInfo ;
			std::vector<std::pair<size_t, int>> countInfo ;

			readUniqueness(uniquenessInfo, refInfo.geneMap, revUniqness) ;
			size_t mostUniqGeneId{0} ;
			size_t leastUniqGeneId{uniquenessInfo.size() - 1} ;

			for(size_t geneId = 0 ; geneId < splatterGeneNames.size() ; geneId++){
				int gcount{0} ;
				for(size_t cellId = 0 ; cellId < allCellNames.size() ; cellId++){
					if(geneCounts[cellId][geneId] > 0){
						gcount++ ;
					}
				}
				countInfo.push_back({geneId, gcount}) ;
				_verbose("\rIn Splatter: Number of gene sum processed : %lu", geneId);
			}
			using ptype = std::pair<size_t, int> ;

			std::sort(countInfo.begin(), countInfo.end(), [](ptype &left, ptype &right) {
	    		return left.second > right.second;
			});	

			std::cerr << "\nIn Splattermode countInfo.size() " << countInfo.size() 
					  << "\t uniquenessInfo.size() " << uniquenessInfo.size() 
					  << "\n" ;

			for(auto it : countInfo){
				//search until you find a proper gene
				while( geneMap.find(uniquenessInfo[mostUniqGeneId].first) == geneMap.end()){
					mostUniqGeneId++ ;
				}

				//std::cerr << "found " << uniquenessInfo[mostUniqGeneId].first << "\n" ; 
				if (mostUniqGeneId <= leastUniqGeneId){
					auto splatterGeneId = it.first ;
					auto spgName = splatterGeneNames[splatterGeneId] ;

					auto realGeneName = uniquenessInfo[mostUniqGeneId].first ;
					auto refGeneId = geneMap[realGeneName] ;

					splatterGeneMap[realGeneName] = spgName ; // keep the real name

					alevin2refMap[splatterGeneId] = refGeneId ; // gene id to real gene id 
					// this one is extra
					alevinGeneIndex2NameMap[splatterGeneId] = realGeneName ; // keep mapping from real gene name  
					
				}
				mostUniqGeneId++ ;
			}
		}else{
			size_t splatterGeneId{0} ;
			// Assign a random gene to the columns of the matrix  
			// Go over geneMap assign a gene that exists
			// and has at least one transcript that has more than 
			// READ_LEN 
			std::vector<std::unordered_set<std::string>> uniqBins ;
			std::vector<int> uniqHist ;

			makeUniquenessBins(uniqBins, geneMap) ;
			uniqHist.resize(uniqBins.size()) ;

			//std::cout << "Made uniqueness bins " << uniqBins.size() << "\t" << uniqHist.size() << "\n"; 
			size_t totalGenes{0} ;

			for(size_t i = 0; i < uniqBins.size(); ++i){
				uniqHist[i] = uniqBins[i].size() ;
				//std::cout << "Bin " << i << " " << uniqHist[i] << "\n" ;
				totalGenes +=uniqHist[i] ; 
			}

			//std::cout << "Total genes " << totalGenes << "\t" << geneMap.size() << "\n"; 
					
			while(splatterGeneId < splatterGeneNames.size()){
				std::random_device rd;
    			std::mt19937 gen(rd());	

				std::discrete_distribution<> d(uniqHist.begin(), uniqHist.end());
				int ind = d(gen) ;
				while(!uniqHist[ind]){
					for(size_t i = 0; i < uniqBins.size(); ++i){
						uniqHist[i] = uniqBins[i].size() ;
						//std::cout << "Bin " << i << " " << uniqHist[i] << "\n" ;
					}	
					ind = d(gen) ;
					//std::exit(1) ;
				}

				auto randIt = uniqBins[ind].begin() ;
				auto selectedGeneName = *randIt ;
				auto selectedGeneId = geneMap[selectedGeneName] ;

				uniqBins[ind].erase(randIt) ;	
				uniqHist[ind] = uniqHist[ind] - 1 ;

				auto spgName = splatterGeneNames[splatterGeneId] ;
				splatterGeneMap[selectedGeneName] = spgName ;
				alevin2refMap[splatterGeneId] = selectedGeneId ;
				alevinGeneIndex2NameMap[splatterGeneId] = selectedGeneName ;

				splatterGeneId++ ;	

				_verbose("\rIn Splatter: Number of genes processed : %lu", splatterGeneId);
			}

		}
		if(splatterGeneMap.size() != splatterGeneNames.size()){
			consoleLog->error("The splatter matrix contains more genes than provided in reference, truncating") ;
			consoleLog->error("splatterGeneMap.size() {} splatterGeneNames.size() {} geneMap.size() {}", 
								splatterGeneMap.size(), 
								splatterGeneNames.size(),
								geneMap.size()
							) ;
			std::exit(2) ;
		}
	}


	auto& gene2transcriptMap = refInfo.gene2transcriptMap ;
	std::map<uint32_t, uint32_t> alevinReverseMap ;
	{
		// create a transcript map
		uint32_t trId{0} ;
		for(auto it : alevin2refMap){
			auto gIt = gene2transcriptMap.find(it.second) ;
			if(gIt != gene2transcriptMap.end()){
				auto& trVec = gIt->second ;
				//std::cerr << "# of transcripts in this gene "<< gIt->first 
				//	      <<  ":\t" << trVec.size() << "\n" ;	
				for(auto tr : trVec){
					alevin2refTranscriptMap[trId] = tr ;
					alevinReverseMap[tr] = trId ;
					trId += 1 ;
				}
			}else{
				consoleLog->error("This should not happen: gene {} not found",gIt->first) ;
			}
		}
		 
	}


	numOfTranscripts = alevin2refTranscriptMap.size() ;
	data.resize(sampleCells, std::vector<T>(numOfTranscripts)) ;
	for(auto& v : data)
    	std::memset(&v[0], 0, sizeof(v[0]) * v.size());
	

	// If we use dbg then the transcripts used above won't be needed  
	if(useDBG){

		if(!util::fs::FileExists(gfaFile.c_str())){
			std::cerr << gfaFile << " should exist EXITING !!"; 
			std::exit(2) ;
		}

		auto rspdFile = simOpts.rspdFile ;
		if(rspdFile != ""){
			if(util::fs::FileExists(rspdFile.c_str())){
				std::cerr << "RSPD file " << rspdFile << "\n" ;
				rspdVec.resize(MAX_FRAGLENGTH) ;
				std::ifstream rspdStream(rspdFile.c_str()) ;
				std::string line ;
				while(std::getline(rspdStream, line)){
					std::vector<std::string> tokens ;
					util::split(line, tokens, ",") ;
					if(tokens.size() != 2)
						continue ;
					uint32_t pos = std::stoul(tokens[0]) ;
					uint32_t prob = std::stoul(tokens[1]) ;
					if(pos < rspdVec.size()){
						rspdVec[pos] = prob ;
					}
				}
			}else{
				std::cerr << "RSPD file is not empty and doesn't exist, going with truncated sampling\n" ;
			}
		}else{
			std::cerr << "RSPD ::: \n" ;
		} 

		dbgPtr = new GFAReader(gfaFile) ;
		dbgPtr->parseFile(refInfo) ;

		// stand alone call to bfh
		// ignore other eqclass stuff 
		std::string eqFileDir = "/mnt/scratch1/hirak/minnow/metadata/hg/" ;

		eqClassPtr = new BFHClass(eqFileDir) ;
		eqClassPtr->loadBFH(
			bfhFile, 
			simOpts.clusterFile, 
			refInfo, 
			cellWhiteListMap, 
			false, 
			cellNoisyMap
		) ;

		consoleLog->info("Loaded the bfh.txt file") ;
		consoleLog->info("The size of probability Vector {}",eqClassPtr->countProbability.size()) ;

		

		preCalculatedSegProb.resize(numOfGenes) ; // per gene
		preSegOreMapVector.resize(numOfGenes) ;
		geneSpecificTrVector.resize(numOfGenes) ;

		// fill segment based probability
		// gene to tr to prob mapping 

		for(uint32_t i = 0; i < numOfGenes; ++i){
			auto it = alevin2refMap.find(i) ;
			if(it != alevin2refMap.end()){
				auto originalGeneId = it->second ;
				//std::cerr << "original gene id "  << originalGeneId << "\n" ;

				auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

				std::unordered_map<size_t, uint32_t> localGeneProb ;
				std::unordered_map<size_t, bool> localSegOreMap ;
				std::unordered_map<size_t, std::vector<trInfo>> localTrVector ;

				//std::cerr << "tr size " << transcriptIds.size() << "\n" ;

				for(auto tid : transcriptIds){
					if(dbgPtr->trSegmentMap.find(tid) != dbgPtr->trSegmentMap.end()){
						auto segVec = dbgPtr->trSegmentMap[tid] ;
						bool everRC{false} ;
						for(auto seg : segVec){
							auto segCount = dbgPtr->eqClassMapCounts[seg] ;
							//auto bfhCount = eqClassPtr->getGeneLevelProbCount(
							//	originalGeneId,
							//	segCount 
							//) ;
							auto bfhCount = eqClassPtr->countProbability[segCount] ;
					
							if(bfhCount != 0){
								// old one
								// localTrVector[seg].emplace_back(tid) ;

								auto tcInfoVec = dbgPtr->eqClassMap[seg][tid] ;
								for(auto tInfo : tcInfoVec){
									if(refInfo.transcripts[tid].RefLength - tInfo.eposInContig <= MAX_FRAGLENGTH){
										localGeneProb[seg] = bfhCount ;
										localTrVector[seg].emplace_back(
											tid,
											tInfo.sposInContig,
											tInfo.eposInContig
										) ;
									}

									if(!tInfo.ore && !everRC){
										everRC = true ;
									}
								}
								if(everRC){
									localSegOreMap[seg] = false ;
								}else{
									localSegOreMap[seg] = true ;
								}
							}else{
								//std::cerr << "segCount " << segCount << "\n" ;
							}
						}
					}
				}

				preCalculatedSegProb[i] = localGeneProb ;
				preSegOreMapVector[i] = localSegOreMap ;
				geneSpecificTrVector[i] = localTrVector ;
				
			}
		} 

		std::cerr << "SPLATTER MODE: After loading bfh prob size " << preCalculatedSegProb.size() << "\n" ;

		cellSegCount.resize(numCells + numOfDoublets) ;



	}	



	size_t numOfDroppedGenes{0} ;
	size_t numOfBadCells{0} ;
	// create transcript level matrix 
	{
		size_t cellId{0} ;
		for(auto&  cellGeneCounts : geneCounts){
			size_t droppedGeneExpression{0} ;
			size_t cellSum{0} ;
			auto barcode = cellNames[cellId] ;
			for(size_t i = 0 ; i < cellGeneCounts.size() ; ++i){
				
				std::string geneName = aleviGeneIndex2NameMap[i] ;

				int geneCount = cellGeneCounts[i] ;
				//int geneLevelTrCount{0} ;
				//int geneLevelExCount{0} ;

				if(geneCount > 0){

					cellSum += geneCount ;

					bool dropThisGene{false} ;
					
					auto originalGeneId = alevin2refMap[i] ;
					auto transcriptIds = refInfo.gene2transcriptMap[originalGeneId] ;

					std::vector<double> probVec(transcriptIds.size()) ;
					std::unordered_map<uint32_t, size_t> transcriptIdMap ;

					for(size_t ttid = 0 ; ttid < transcriptIds.size() ; ++ttid){
						transcriptIdMap[transcriptIds[ttid]] = ttid ; 
					}

					
					// follow weibull 
					if(useWeibull){
						std::random_device rdt;
						std::mt19937 gent(rdt());

						// this is the learned distribution paeameters from the paper 
						// https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005761
						std::weibull_distribution<> dwt(0.44,0.6);

						std::random_device rd2;
						std::mt19937 g2(rd2()) ;
						std::vector<int> tIndex(transcriptIds.size()) ;
						std::iota(tIndex.begin(), tIndex.end(), 0);
						std::shuffle(tIndex.begin(), tIndex.end(), g2) ;

						for(size_t j = 0; j < transcriptIds.size() ; ++j ){
							probVec[tIndex[j]] = dwt(gent) ;
						}
					}else if(useDBG){

						dropThisGene = false ;

						auto segCountHist = preCalculatedSegProb[i] ;

						//if(i == 53593){
					    //std::cerr << "segCountHist Size " << segCountHist.size() << "\n" ;
						//}	
						
						if(segCountHist.size() == 0){
							//std::cerr << "No seg hist for this gene " << i << "\n" ; 
							droppedGeneExpression += geneCount ;
							dropThisGene = true ;
							geneCount = 0 ;
						}

						if(!dropThisGene){
							std::vector<size_t> segIndex(segCountHist.size()) ;
							std::vector<uint32_t> segProbVec(segCountHist.size()) ;
							std::vector<int> segCounts(segCountHist.size(),0) ;
							
							size_t ind{0} ;
							for(auto it : segCountHist){
								segIndex[ind] = it.first ;
								segProbVec[ind] = it.second ;
								ind++ ;
							}

							std::random_device rdt4;
							std::mt19937 gent4(rdt4());

							std::unordered_map<size_t, int> segCountMap ;
							
							std::discrete_distribution<> dmt(segProbVec.begin(), segProbVec.end()) ;
							for(int j = 0 ; j < geneCount ; ++j){
								++segCounts[dmt(gent4)] ;
							}
							for(size_t j = 0 ; j < segCounts.size() ; ++j){
								segCountMap[segIndex[j]] = segCounts[j] ;
							}

							cellSegCount[cellId][i] = segCountMap ;

						}
					}
					else{
						consoleLog->error("This version takes either of --useEqClass/--weibull/--dbg") ;
						std::exit(1) ;
					}


					trueGeneCounts[cellId][i] = geneCount ;
					
					if(useDBG){
						continue ;	
					}

					auto sum_prob = std::accumulate(probVec.begin(), probVec.end(), 0.0) ;
					if(sum_prob == 0){
						dropThisGene = true ;
					}

					if (dropThisGene){
						numOfDroppedGenes++ ;
						continue ;
					}
					// fill up the cell to transcriptome matrix 

					std::random_device rdt3;
					std::mt19937 gent3(rdt3());

					std::discrete_distribution<> dmt(probVec.begin(), probVec.end()) ;
					// saple from this distribution geneCount times
					std::vector<int> transcriptCounts(transcriptIds.size(), 0) ; 
					for(int j=0 ; j < geneCount ; ++j){
						++transcriptCounts[dmt(gent3)] ;
					}
					// total gene count 
					// This should not be needed
					T totCount{0} ;
					//uint32_t lastAlevinTid{0} ;
					for(size_t j = 0; j < transcriptIds.size(); ++j){
						auto tid = transcriptIds[j] ;
						auto alevinTid = alevinReverseMap[tid] ;
						data[cellId][alevinTid] = transcriptCounts[j] ;
						totCount += data[cellId][alevinTid] ;
						//lastAlevinTid = alevinTid ;
					}
				}
			}
			cellId++ ;
			if(numOfDroppedGenes > numOfGenes - 100){
				numOfBadCells++ ;
			}
			numOfDroppedGenes = 0 ;
			_verbose("\rIn Splatter: Number of cells processed : %lu", cellId);
		
			if(useDBG){
				std::cerr << " Num of genes dropped due to dbg " << droppedGeneExpression 
				<< " out of " <<  cellSum << "\n" ;

			}	
		}

	}
	if(numOfBadCells > 0)
		consoleLog->warn("{} cells have less than 100 genes present",numOfBadCells);
	
	numOfWhiteListedCells = cellWhiteListMap.size() ;
	//dumpSplatterMaps(outDir) ;
	//dump the generated cell cols 
	//which are genes 
	std::string cellColFile = outDir + "/alevin/quants_mat_cols.txt" ;
	std::ofstream cellColStream(cellColFile.c_str()) ;
	
	for(size_t splatterGeneId = 0 ; splatterGeneId < alevinGeneIndex2NameMap.size() ; ++splatterGeneId){
		cellColStream << alevinGeneIndex2NameMap[splatterGeneId] << "\n" ; 
	}
}


template<typename T>
void DataMatrix<T>::dumpGeneCountMatrix(
	std::string& fileName,
	std::unordered_set<uint32_t>& emptyCells
){
	std::ofstream geneMatrixStream(fileName) ;

	size_t numOfExpressedGenes{0} ;

	if(trueGeneCounts.size() != 0){
		for(size_t cell_id = 0 ; cell_id < trueGeneCounts.size() ; ++cell_id){
			if(emptyCells.find(cell_id) != emptyCells.end()){
				continue ;
			}
			auto v = trueGeneCounts[cell_id] ;
			for(auto c : v){
				geneMatrixStream << c << "," ;
				if(c > 0){
					numOfExpressedGenes++ ;
				}
			}
			geneMatrixStream << "\n" ;
		}

	}else{
		for(size_t cell_id = 0 ; cell_id < geneCounts.size() ; ++cell_id){
			if(emptyCells.find(cell_id) != emptyCells.end()){
				continue ;
			}
			auto v = geneCounts[cell_id] ;
			for(auto c : v){
				geneMatrixStream << c << "," ;
				if(c > 0){
					numOfExpressedGenes++ ;
				}
			}
			geneMatrixStream << "\n" ;
		}
	}

	//std::cerr<< "DEBUG: ==== At the end Expressed Gene count "<< numOfExpressedGenes << "\n" ;
}


template<typename T>
void DataMatrix<T>::dumpTrueCellNames(
	std::string& fileName
){
	std::ofstream cellNameStream(fileName) ;

	for(auto cn : cellNames){
		cellNameStream << cn << "\n" ;
	}

	//std::cerr<< "DEBUG: ==== At the end Expressed Gene count "<< numOfExpressedGenes << "\n" ;
}

template<typename T>
void DataMatrix<T>::dumpGeneCountMatrixUnmolested(
	std::string& fileName
){
	std::ofstream geneMatrixStream(fileName) ;
	for(auto& v : geneCounts){
		for(auto c : v){
			geneMatrixStream << c << "," ;
		}
		geneMatrixStream << "\n" ;
	}
}

template<typename T>
void DataMatrix<T>::dumpIntronCount(
	//std::string& fileName
){
	//std::ofstream geneMatrixStream(fileName) ;
	uint32_t intronCountSum{0} ;
	for(auto& v : intronCountMap){
		
		if(cellWhiteListMap["GTATTCTGTAGTGAAT"] == v.first){
			for(auto& i : v.second)
				intronCountSum += i.second ;
		}
		
	}
	std::cout << "BEFORE INTRON COUNT FOR GTATTCTGTAGTGAAT:  " << intronCountSum << "\n" ;
}

template<typename T>
DataMatrix<T>::DataMatrix(
	SimulateOptions& simOpts,
    Reference& refInfo,
	std::shared_ptr<spdlog::logger>& consoleLogIn
){

	consoleLog = consoleLogIn ; 
	
    if (simOpts.alevinMode){
        loadAlevinData(simOpts, refInfo) ;
    }else if(simOpts.splatterMode){
        loadSplatterData(simOpts, refInfo) ;
    }else{
		consoleLog->info("Please specify either of --splatter-mode or --alevin-mode") ;
	}
    //std::cerr << "Parsed the matrix\n" ; 
}

template class DataMatrix<int> ;
template class DataMatrix<double> ;


//template 
//void pupulateGeneCountFromTxt<int>(
//	std::string& countFile, 
//	std::vector<std::vector<int>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes	
//); 
//template 
//void pupulateGeneCountFromTxt<double>(
//	std::string& countFile, 
//	std::vector<std::vector<double>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes	
//); 
//
//template
//void populateGeneCountFromBinary<int>(
//	std::string& countFileBinary, 
//	std::vector<std::vector<int>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes
//);
//
//template
//void populateGeneCountFromBinary<double>(
//	std::string& countFileBinary, 
//	std::vector<std::vector<double>>& geneCount,
//	uint32_t& numCells, 
//	size_t& numOfGenes
//);
