#include "GFAReader.hpp"
#include "macros.hpp"

// Taken from here 
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


std::vector<std::pair<size_t, bool>> explode(const std::string str, const char& ch) {
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
          std::cerr << "tried to convert " << next << " into a long long\n";
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
    std::cerr << "impossible is the opposite of possible " << next << "\n";
    std::cerr << "The line is " << str << "\n";
    result.emplace_back(std::stoll(next),
                        true); // this case shouldn't even happen
  }
  return result;
}


void GFAReader::updateEqClass(
  std::string& transcriptName, 
  std::vector<std::pair<size_t, bool>>& contigVec,
  Reference& refInfo 
  ){

	uint32_t overlap{READ_LEN + 1} ;
    auto it = refInfo.transcriptNameMap.find(transcriptName) ; 
    if(it != refInfo.transcriptNameMap.end()){

      auto tid = it->second ;
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
	Reference& refInfo 
){
    std::cerr << "Start loading segments... \n" ;
    size_t contigCt{0} ;
    std::string ln, tag, id, value ;

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

            unitigMap[contigId] = value ;
        }
        contigCt++ ;
    }

    std::cerr << "Saw " << contigCt << " contigs in total, unitigMap.size(): " << unitigMap.size() << "\n" ; 
	std::cerr << "Max contig id " << maxContigId << "\n" ; 
    std::cerr << "Starting to load paths \n" ;

    file.reset(new std::ifstream(gfaFileName_)) ;

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
        updateEqClass(id, contigVec, refInfo) ;

        
    }


	
    
    std::cerr << "Done with GFA \n" 
              << "Equivalece class size " << eqClassMap.size() 
			  << "\ttrSegmentMap size " << trSegmentMap.size() 
			  << "\ttranscript map size " << refInfo.transcriptNameMap.size() << "\n" ;

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
	}

	std::unordered_set<size_t> removeKeys ;	

	for(auto unitig : eqClassMap){
		auto id = unitig.first ;
		auto trInfoMap = unitig.second ;
		bool keepIt{false} ;

		for(auto tinfo : trInfoMap){
			auto tid = tinfo.first ;
			auto tidVec = tinfo.second ;
			for(auto info : tidVec){
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


	for(auto unitig : eqClassMap){

		auto id = unitig.first ;

		auto trInfoMap = unitig.second ;

		eqClassMapCounts[id] = static_cast<uint32_t>(trInfoMap.size()) ;

		for(auto tinfo : trInfoMap){
			auto tid = tinfo.first ;
			trSegmentMap[tid].push_back(id) ;
		}
	}
   
    std::cerr << "Done Filtering \n" 
              << "Equivalece class size " << eqClassMap.size() 
			  << "\ttrSegmentMap size " << trSegmentMap.size() 
			  << "\ttranscript map size " << refInfo.transcriptNameMap.size() << "\n" ;

    //std::cerr << "After filtering eqclassSize " << eqClassMap.size() << " \n" ;

}


