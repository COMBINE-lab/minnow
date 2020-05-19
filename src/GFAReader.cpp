#include "GFAReader.hpp"
#include "macros.hpp"

// Parts of the methods are taken from here 
// https://github.com/COMBINE-lab/pufferfish/blob/master/src/GFAConverter.cpp

bool is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

std::string getGencodeTranscript(std::string v){
    std::vector<std::string> tokens ; 
    util::split(v, tokens, "|") ;
    return tokens[0] ;
}


std::vector<std::pair<size_t, bool>> GFAReader::explode(
  const std::string str, 
  const char& ch
) {
  std::string next;
  std::vector<std::pair<size_t, bool>> result;
  // For each character in the string
  for (auto it = str.begin(); it != str.end(); it++) {
    // If we've hit the terminal character
    if (*it == '+' or *it == '-') {
      bool orientation = true;
      // If we have some characters accumulated
      // Add them to the result vector
      if (!next.empty()) {
        if (*it == '-') {
          orientation = false;
        }
        try {
          uint64_t nid = std::stoll(next);
          result.emplace_back(nid, orientation);
        } catch (std::exception& e) {
          // not a numeric contig id
          consoleLog->error("tried to convert {}"
                            " into a long long",
                            next
          );
          std::exit(1);
        }
        next.clear();
      }
    } else if (*it != ch) {
      // Accumulate the next character into the sequence
      next += *it;
    }
  }
  if (!next.empty()) {
    consoleLog->error("impossible is the opposite of possible {} "
                      "The line is {} ",
                      next,str
                     );
    result.emplace_back(std::stoll(next),
                        true); // this case shouldn't even happen
  }
  return result;
}


void GFAReader::updateEqClass(
  std::string& transcriptName, 
  std::vector<std::pair<size_t, bool>>& contigVec,
  Reference& refInfo,
  size_t overlap
  ){

    auto it = refInfo.transcriptNameMap.find(transcriptName) ; 
    if(it != refInfo.transcriptNameMap.end()){

      auto tid = it->second ;

      
      if(tid == 11393){
        //std::cout << "[DEBUG]-----I'm here" << transcriptName <<"\n";
      }

	  size_t numOfMultiMapped{0} ;
      uint32_t tStart = 0 ;
      for(auto contigInfo : contigVec){
          auto contigId = contigInfo.first ;
          auto ore = contigInfo.second ;
          if(unitigMap.find(contigId) != unitigMap.end()){
            auto sizeOfUnitig = unitigMap[contigInfo.first].size() ;

            if(sizeOfUnitig < READ_LEN){
              std::cerr << sizeOfUnitig <<  " --- possible twopaco bug !!!!\n" ;
              std::exit(1) ;
            }

            auto tEnd = tStart + sizeOfUnitig - 1 ;
            eqClassMap[contigId][tid].emplace_back(tStart, tEnd, ore) ;
            if(eqClassMap[contigId][tid].size() > 1){
              numOfMultiMapped++ ;
            }

            //if(transcriptName == "ENST00000621744.4"){
            //  std::cout << "\n[DEBUG] tr: " << transcriptName << "\t"
            //            << "tid " << tid << "\t" << tStart << "\t"
            //            << tEnd << " contigId " << contigId;
            //}

            //trSegmentMap[tid].push_back(contigId) ;
            tStart += (sizeOfUnitig - overlap) ;
          }
      }
  
    }

}

void GFAReader::readUnitigs(){
	std::cerr << "Filling unitigs... \n" ;
	size_t contigCt{0} ;
    std::string ln, tag, id, value ;

	file.reset(new std::ifstream(gfaFileName_)) ;


    while(std::getline(*file, ln)){
        char fastC = ln[0] ;
        if(fastC != 'S')
            continue ; 
        
        std::vector<std::string> tokens  ; 
        util::split(ln, tokens, "\t") ;

        id = tokens[1] ;
        value = tokens[2] ;
        if(is_number(id)){
            size_t contigId = std::stoll(id) ;
            unitigMap[contigId] = value ;
        }
        contigCt++ ;
    }

	std::cerr << "Saw " << contigCt << " contigs in total, unitigMap.size(): " << unitigMap.size() << "\n" ;

}

void GFAReader::parseFile(
	Reference& refInfo,
  std::string& outDir
){

  consoleLog->info("Parsing GFA file {}", gfaFileName_);
  consoleLog->info("Start loading segments...") ;
  size_t contigCt{0} ;
  std::string ln, tag, id, value ;


  //We need to calculate the overlap size on the fly
  //if and only if the new pufferize TwoPaCo is run
  //otherwise it is not important, we can go without
  
  //ThreePaCo 
  //overlapSize = READ_LEN - 1
  //TwoPaCo
  //overlapSize = READ_LEN + 1

  bool foundOverlapSize{false};
  // size_t overlapsize = READ_LEN-1;
  size_t overlapsize = READ_LEN + 1;
  consoleLog->info("Predicted overlap size: {}", overlapsize);

	file.reset(new std::ifstream(gfaFileName_)) ;

	size_t maxContigId{0} ;
  while(std::getline(*file, ln)){
      char fastC = ln[0] ;
      if(fastC != 'S')
          continue ; 
      
      std::vector<std::string> tokens  ; 
      util::split(ln, tokens, "\t") ;

      id = tokens[1] ;
      value = tokens[2] ;
      if(is_number(id)){
          size_t contigId = std::stoll(id) ;
          if(contigId > maxContigId)
            maxContigId = contigId ;
          if (value != "*"){
            unitigMap[contigId] = value ;
          }
      }
      contigCt++ ;
  }

  consoleLog->info("Saw {} segment lines, number of unitigs {}",contigCt, unitigMap.size()) ; 
	//std::cerr << "Max contig id " << maxContigId << "\n" ; 
  //std::cerr << "Starting to load paths \n" ;

  // reset the file
  file.reset(new std::ifstream(gfaFileName_)) ;
  // find overlap size
  while(std::getline(*file, ln)){
      char fastC = ln[0] ;
      if(fastC != 'P')
          continue ;
      std::vector<std::string> tokens ;
      util::split(ln, tokens, "\t") ;
      if(tokens.size() != 4){
          continue ;
      }
      // sparse this to get the transcript name
      id = getGencodeTranscript(tokens[1]) ;
      if(id == ""){
          continue ;
      }

      // A valid line
      auto pvalue = tokens[2] ;
      std::vector<std::pair<size_t, bool>> contigVec = explode(pvalue, ',') ;
      // calculate overlap size
      if(!foundOverlapSize and contigVec.size() > 2){
        auto selem1 = contigVec[0];
        auto selem2 = contigVec[1];
        std::string elem1, elem2;
        if (!selem1.second){
          elem1 = util::revcomp(unitigMap[selem1.first]);
        }else{
          elem1 = unitigMap[selem1.first];
        }
        if (!selem2.second){
          elem2 = util::revcomp(unitigMap[selem2.first]);
        }else{
          elem2 = unitigMap[selem2.first];
        }
        while(elem2.substr(0,overlapsize) != elem1.substr(elem1.size()-overlapsize)){
          overlapsize++;
          if(overlapsize == elem1.size() || overlapsize == elem2.size()){
            consoleLog->error("GFA is ill-constructed") ;
            //std::cout << elem1 << "\t" << elem2 << "\t" << overlapsize << "\n";
            std::exit(1);
          }
        }
        consoleLog->info("Calculated overlap size {}", overlapsize);

        foundOverlapSize = true;
      }
      if (foundOverlapSize){
          break;
      }
  }

  // now update the eqclasses
  file.reset(new std::ifstream(gfaFileName_)) ;
  consoleLog->info("Start loading paths...");
  while(std::getline(*file, ln)){
      char fastC = ln[0] ;
      if(fastC != 'P')
          continue ;
      std::vector<std::string> tokens ;
      util::split(ln, tokens, "\t") ;
      if(tokens.size() != 4){
          continue ;
      }
      // parse this to get the transcript name
      id = getGencodeTranscript(tokens[1]) ;

      if(id == ""){
          continue ;
      }

      // A valid line
      auto pvalue = tokens[2] ;
      std::vector<std::pair<size_t, bool>> contigVec = explode(pvalue, ',') ;

      updateEqClass(id, contigVec, refInfo, overlapsize) ;

  }

	// filter segements that is does not have 
	// any transcript where it is within 
	// MAX_FRAGLEN within 3'

	std::unordered_map<uint32_t, uint32_t> distanceFromEndMap ;
	for(auto tr : refInfo.transcriptNameMap){
		auto tid = tr.second ;
		auto refLen = refInfo.transcripts[tid].RefLength ;
		if(refLen > MAX_FRAGLENGTH){
			distanceFromEndMap[tid] = refLen - MAX_FRAGLENGTH ;
		}else{
			distanceFromEndMap[tid] = 0 ;
		}
    // if(tid == 11393){
    //   std::cout << "[DEBUG]-----" << distanceFromEndMap[tid]<<"\n";
    // }
	}

	std::unordered_set<size_t> removeKeys ;	

  // debug
  //size_t unitig_id_to_track{0};
	for(auto unitig : eqClassMap){
		auto id = unitig.first ;
		auto trInfoMap = unitig.second ;
		bool keepIt{false} ;

		for(auto tinfo : trInfoMap){
			auto tid = tinfo.first ;
			auto tidVec = tinfo.second ;
			for(auto info : tidVec){
        // [DEBUG]
        //if(id == 38557 and tid == 11393){
        //   std::cout << "[DEBUG]-----" << "\t"
        //             <<info.sposInContig << "\t"
        //             <<info.eposInContig
        //             << "\t" << distanceFromEndMap[tid] << "\t"
        //             << refInfo.transcripts[tid].RefLength <<  "\n" ;
        //}
				if(info.eposInContig >= distanceFromEndMap[tid]){
					keepIt = true ;
					break ;
				}
			}
			if(keepIt)
				break ;
		}

		if(!keepIt){
			removeKeys.insert(id) ;
		}
	}

	for(auto k : removeKeys){
		eqClassMap.erase(k) ;
	}


  std::vector<bool> includedTranscripts(refInfo.numOfTranscripts, false) ;
	for(auto unitig : eqClassMap){

		auto id = unitig.first ;

		auto trInfoMap = unitig.second ;

		eqClassMapCounts[id] = static_cast<uint32_t>(trInfoMap.size()) ;

		for(auto tinfo : trInfoMap){
			auto tid = tinfo.first ;
			trSegmentMap[tid].push_back(id) ;
      includedTranscripts[tid] = true;
		}
	}

  if(refInfo.duplicateList.size() and (trSegmentMap.size() < refInfo.transcriptNameMap.size())){
    consoleLog->info("Trying to rescue some transcripts that might be dropped due to duplicates");
    for(size_t i = 0 ; i < refInfo.numOfTranscripts ; ++i){
      if(!includedTranscripts[i]){
        auto trIt = refInfo.duplicateList.find(i);
        if(trIt != refInfo.duplicateList.end()){
          trSegmentMap[i] = trSegmentMap[trIt->second] ;
          includedTranscripts[i] = true ;
        }
      }
    }
  }
  
  consoleLog->info("Done with GFA "
                   "Equivalence class size {} "
                   "Segment map size after filtering {} "
                   "number of transcripts {}",
                   eqClassMap.size(),
                   trSegmentMap.size(),
                   refInfo.transcriptNameMap.size()
  );
  
  if(trSegmentMap.size() < refInfo.transcriptNameMap.size()){
    consoleLog->warn("{} transcripts will not be included as " 
                    "they don't have suitable segments",
                    refInfo.transcriptNameMap.size() - trSegmentMap.size()
    );
    std::string skippedGenesFile = outDir + "/skipped_transcripts.txt" ;
		std::ofstream skippedGeneStream(skippedGenesFile.c_str()) ;
		
    for(size_t i = 0; i < includedTranscripts.size(); ++i){
      if (!includedTranscripts[i]){
        skippedGeneStream << refInfo.transcripts[i].RefName << "\t" 
        << refInfo.transcripts[i].RefLength << "\n" ;
      }
    }
  }


}


