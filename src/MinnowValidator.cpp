#include "CLI/CLI.hpp"
#include "clipp.h"

#include "FastxParser.hpp"
#include "FASTAParser.hpp"
#include "GFAReader.hpp"
#include "MinnowUtil.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <unordered_map>
#include <map>
#include "ScopedTimer.hpp"
#include "string_view.hpp"
#include "edlib.h"
#include "macros.hpp"

#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)

#define READ_LEN 100

struct ValidateOpt{
  std::string gfaFile{""} ;
  std::string fastqFile{""} ;
  std::string referenceFile{""} ;
  std::string outFile{""} ;
  int numThreads{0} ;
  int edit_max_lim{0} ; 

} ;

void split2(const std::string& str, std::vector<std::string>& tokens, const std::string& delim){
    const std::string whiteSpace = " " ;
    size_t prev = 0, pos = 0;
    do
      {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty() and token != whiteSpace) tokens.push_back(token);
        prev = pos + delim.length();
      }
    while (pos < str.length() && prev < str.length());
}


stx::string_view sampleSequence(
    Transcript& tr,
    size_t start_pos 
){
    
    return stx::string_view(
            tr.Sequence() + start_pos,
            READ_LEN
    );
}


void refValidate(
	std::string& fastqFile,
	std::string& referenceFile,
  int& edit_max_lim,
	std::string& outFile
){
	std::vector<Transcript> transcripts; 
	std::cout << "Reading reference file\n" ; 	
	FASTAParser fastaParser(referenceFile) ;
    fastaParser.populateTargets(transcripts) ;

	std::unordered_map<std::string, size_t> trMap ;
	trMap.reserve(transcripts.size()) ;
	for(size_t i = 0 ; i < transcripts.size(); ++i){
		trMap[transcripts[i].RefName] = i ;
	}
	std::map<int, int> editDistanceMap ;

	{
		ScopedTimer st ;
		std::vector<std::string> read_file = {fastqFile} ;
		fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
		parser.start() ;

		// Get the read group by which this thread will
		// communicate with the parser (*once per-thread*)
		size_t rn{0};
		// size_t kmer_pos{0};
		auto rg = parser.getReadGroup();
		while (parser.refill(rg)) {
		// Here, rg will contain a chunk of read pairs
		// we can process.
			for (auto& rp : rg) {
				// kmer_pos = 0;
				++rn;
        if(rn % 1000 == 0)
          _verbose("\rNumber of reads processed : %lu", rn);
				auto& r1 = rp.seq;
				auto& header = rp.name ;
				std::vector<std::string> headerVec; 
				split2(header, headerVec, ":") ;
				std::string transcriptName = headerVec[1] ;
        if (transcriptName.find("polyA") != std::string::npos){
          //std::cout << "Here found polyA, this should not be present in the reference\n" ;
          //std::cout << "But getting this " << trMap[transcriptName] << "\n" ;
          //std::exit(1) ;
        }
        size_t position{0} ;

        if(headerVec.size() == 5){
				  position = std::stol(headerVec[3]) ;
        }else{

				  position = std::stol(headerVec[2]) ;
        }


				auto& index = trMap[transcriptName] ;
				stx::string_view refSeq = sampleSequence(transcripts[index], position) ;

				AlignerEngine ae_ ;
				ae_(
					r1.data(), 
					static_cast<int>(r1.size()), 
					refSeq.data(), static_cast<int>(refSeq.size()),
				 	edlibNewAlignConfig(edit_max_lim, EDLIB_MODE_NW, EDLIB_TASK_LOC)
				);
				auto& result = ae_.result() ;
				auto it = editDistanceMap.find(result.editDistance) ;
				if(it != editDistanceMap.end()){
					it->second += 1 ;
				}else{
					editDistanceMap.insert({result.editDistance, 1}) ;
				}	
			}
		}
    std::cout << "Done with " << rn << " reads \n" ;

    parser.stop() ;
	}
	std::ofstream dictFile{outFile.c_str()} ;
	for(auto it: editDistanceMap){
		dictFile << it.first << "\t" << it.second << "\n" ; 
	}
	dictFile.close() ;

}

void gfaValidate(
  std::string& gfaFile, 
  std::string& fastqFile,
  int& edit_max_lim,
	std::string& outFile

){
  GFAReader gfaObj(gfaFile) ;

  gfaObj.readUnitigs() ;

  std::map<int, int> editDistanceMap ;

	{
		ScopedTimer st ;
		std::vector<std::string> read_file = {fastqFile} ;
		fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
		parser.start() ;

		// Get the read group by which this thread will
		// communicate with the parser (*once per-thread*)
		size_t rn{0};
		// size_t kmer_pos{0};
		auto rg = parser.getReadGroup();
		while (parser.refill(rg)) {
		// Here, rg will contain a chunk of read pairs
		// we can process.
			for (auto& rp : rg) {
				// kmer_pos = 0;
				++rn;
        if(rn % 1000 == 0)
          _verbose("\rNumber of reads processed : %lu", rn);
				auto& r1 = rp.seq;
				auto& header = rp.name ;
				std::vector<std::string> headerVec; 
				split2(header, headerVec, ":") ;
				std::string transcriptName = headerVec[1] ;



        if(headerVec.size() < 6){
				  std::cerr << "something is wrong" ;
          std::exit(1) ;
        }

        size_t contigId = std::stoull(headerVec[1]) ;
        size_t position = std::stoull(headerVec[2]) ;
        bool ore ;
        if(headerVec[3] == "-"){
          ore = false ;
        }else{
          ore = true ;
        }
        std::string refSeq{""} ;
        std::string refSeqPart{""} ;
        if(!ore){
          refSeq = util::revcomp(gfaObj.unitigMap[contigId]) ;
        }else{
          refSeq = gfaObj.unitigMap[contigId] ; 
        }

        refSeqPart = refSeq.substr(position, READ_LEN) ;

				AlignerEngine ae_ ;
				ae_(
					r1.data(), 
					static_cast<int>(r1.size()), 
					refSeqPart.data(), static_cast<int>(refSeqPart.size()),
				 	edlibNewAlignConfig(edit_max_lim, EDLIB_MODE_NW, EDLIB_TASK_LOC)
				);
				auto& result = ae_.result() ;
				auto it = editDistanceMap.find(result.editDistance) ;
				if(it != editDistanceMap.end()){
					it->second += 1 ;
				}else{
					editDistanceMap.insert({result.editDistance, 1}) ;
				}	
			}
		}
    std::cout << "Done with " << rn << " reads \n" ;

    parser.stop() ;
	}
	std::ofstream dictFile{outFile.c_str()} ;
	for(auto it: editDistanceMap){
		dictFile << it.first << "\t" << it.second << "\n" ; 
	}
	dictFile.close() ;


}

void validate(ValidateOpt& valOpts){
  if(valOpts.gfaFile == ""){
    refValidate(
      valOpts.fastqFile,
      valOpts.referenceFile,
      valOpts.edit_max_lim,
      valOpts.outFile
    ) ;
  }else{
    gfaValidate(
      valOpts.gfaFile,
      valOpts.fastqFile,
      valOpts.edit_max_lim,
      valOpts.outFile
    ) ;
  }

}

int main(int argc, char* argv[]) {
  using namespace clipp ;
  using std::cout ;
  enum class mode {help, validate} ;

  mode selected = mode::help ;

  ValidateOpt valOpts ;
  

  auto simulateMode = (
    command("validate").set(selected, mode::validate),

    
    (option("--gfa") &
    value("alevin-dir", valOpts.gfaFile)) %
    "fastq file",
    
    (option("-f", "--fastq-file") &
    value("alevin-dir", valOpts.fastqFile)) %
    "fastq file",
    
    (option("-t", "--reference") &
    value("transcript", valOpts.referenceFile)) %
    "Directory to contain alevin related files",
    
    (option("-m", "--max_dist") &
    value("transcript", valOpts.edit_max_lim)) %
    "Directory to contain alevin related files",
    
    (option("-p", "--num-threads") &
    value("number-threads", valOpts.numThreads)) %
    "number of Threads needed for the parsing",
    
    (option("-o", "--output") &
    value("output", valOpts.outFile)) %
    "output file to dump the results"
  );

  auto cli = (
    (simulateMode | 
     command("--help").set(selected,mode::help) |
     command("-h").set(selected,mode::help) |
     command("help").set(selected,mode::help)
    ),
    option("-v", "--version").call([]{std::cout << "version 0.1.0\n\n";}).doc("show version")
  );

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "validate");
    return 1;
  }

  if(res){
    switch(selected){
      case mode::validate: validate(valOpts) ; break ; 
      case mode::help: std::cout << make_man_page(cli, "validate") ; 
    }
  }
  std::exit(0) ;

}