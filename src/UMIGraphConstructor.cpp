#include <unordered_map>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "MinnowFS.hpp"
#include "TranscriptGroup.hpp"
#include "MinnowUtil.hpp"
#include "Graph.hpp"

#include "CLI/CLI.hpp"
#include "clipp.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>

#include <sparsepp/spp.h>

#include "cuckoohash_map.hh"

using BFSType = std::unordered_map<std::string, spp::sparse_hash_map<TranscriptGroup, std::unordered_map<std::string,int>, TranscriptGroupHasher>> ;
using TxGrpType = spp::sparse_hash_map<TranscriptGroup, std::unordered_map<std::string,int>, TranscriptGroupHasher> ;

uint32_t stoui32(const std::string& s)
{
    std::istringstream reader(s);
    uint32_t val = 0;
    reader >> val;
    return val;
}

int cleanLineToInt(std::string& line){
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    return std::stoi(line) ;
}




void parseBFH(
    std::string& alevinDir,
    std::string& t2g_file,
    BFSType& bfh 
){
    std::cerr << "transcript to gene mapping can be found here " << t2g_file << "\n" ;
    std::cerr << "Start reading the BFH file\n" ; 
    // location to the bfh file
    if(! util::fs::DirExists(alevinDir.c_str())){
		std::cerr << "Alevin directory does not exists\n" ;
		std::exit(1) ; 
	}
    auto bfh_file = alevinDir + "/alevin/bfh.txt" ;
    {
        std::ifstream fileStream(bfh_file) ;
        // read first three lines
        std::string line ;
        // number of transcripts
        std::getline(fileStream, line) ;
        int T = cleanLineToInt(line) ;
        // number of cells 
        std::getline(fileStream, line) ;
        int C = cleanLineToInt(line) ;
        // number of experiments 
        std::getline(fileStream, line) ;
        int E = cleanLineToInt(line) ;
        // Now in loop read read rest of the lines 

        //read all the transcripts 
        std::vector<std::string> tname ;
        tname.resize(T) ;
        std::cerr << "Reading "<<T<< " transcripts \n" ;
        for(int i = 0; i < T; ++i){
            std::getline(fileStream, line) ;
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            tname[i] = line ;
        }

        // read the set of predicted cells 
        std::vector<std::string> cellbarcodeMap ;
        std::cerr << "Reading "<<C<< " cells \n" ;
        cellbarcodeMap.resize(C) ; 
        for(int i = 0; i < C; ++i){
            std::getline(fileStream, line) ;
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            cellbarcodeMap[i] = line ;
        }

        // size of bfh is known
        bfh.reserve(static_cast<uint32_t>(C)) ;

        // read line by line to construct three 
        // level hash 
        // cell barcode --> tgGroup --> UMI --> count


        std::cerr << "Will read "<<E<< "  equivalence classes\n" ;
        uint32_t lineNum{0} ;
        while(std::getline(fileStream, line)){
            std::vector<std::string> toks ;
            util::split(line, toks, "\t") ;
            int num_of_labels = std::stoi(toks[0]) ;
            std::vector<uint32_t> txps ;
            txps.resize(num_of_labels) ;
            for(int i = 0; i < num_of_labels; ++i){
                txps[i] = stoui32(toks[1 + i]) ;
            }
            std::sort(txps.begin(), txps.end()) ;
            TranscriptGroup txGrp(txps) ;
            int total_number_of_reads = std::stoi(toks[num_of_labels + 1]) ;

            int currIdx = num_of_labels + 2 ;
            int num_bc = std::stoi(toks[currIdx]) ;

            int read_validator = 0 ;

            for(int i = 0; i < num_bc; ++i){
                currIdx += 1 ;
                std::string bc_name = cellbarcodeMap[std::stoi(toks[currIdx])] ;
                currIdx += 1 ;
                int num_umi = std::stoi(toks[currIdx]) ;
                int num_reads = 0 ;

                for(int j = 0; j < num_umi ; ++j){
                    currIdx += 2 ;
                    num_reads += std::stoi(toks[currIdx]) ;
                    std::string umi = toks[currIdx-1] ;
                    auto findCellBC = bfh.find(bc_name) ;
                    if (findCellBC != bfh.end()){
                        // this cell barcode is seen
                        // reference of the object 
                        auto& thisCellBC = findCellBC->second ;
                        // find if tgGrp is there 
                        auto findTgGrp = thisCellBC.contains(txGrp) ;
                        if (findTgGrp){
                            // tgGrp is also present
                            // reference of the object 
                            auto& thisUMIMap = thisCellBC[txGrp] ;
                            // find if umi is present
                            auto findUMI = thisUMIMap.find(umi) ;
                            if(findUMI != thisUMIMap.end()){
                                // UMI is present 
                                // increase the count 
                                findUMI->second = findUMI->second + std::stoi(toks[currIdx]) ;
                            }else{
                                thisUMIMap[umi] = std::stoi(toks[currIdx])  ;
                            }
                        }else{
                            // tgGrp is absent
                            // we can create the umi object
                            // and add this group 
                            std::unordered_map<std::string, int> tmpMap ;
                            tmpMap[umi] = std::stoi(toks[currIdx]) ; 
                            thisCellBC[txGrp] = tmpMap ;
                        }
                    }else{
                        // this cell barcode is not present 
                        std::unordered_map<std::string, int> tmpMap ;
                        tmpMap[umi] = std::stoi(toks[currIdx]) ; 
                        TxGrpType tmpGrpMap ;
                        tmpGrpMap[txGrp] = tmpMap ;
                        bfh[bc_name] = tmpGrpMap ;
                    }
                }
                    read_validator += num_reads ;
            }
            assert(total_number_of_reads == read_validator) ;
            lineNum++ ;
            if (lineNum%100 == 0){
                std::cout << lineNum << "\r" ;
            }
        }
        std::cout << "\nDone reading equivalence classes \n\n" ;

    }

}


void constructGraph(
    std::string& alevinDir,
    std::string& t2g_file
){

    // read line by line to construct three 
    // level hash 
    // cell barcode --> tgGroup --> UMI --> count 
    BFSType bfh ;

    parseBFH(
        alevinDir,
        t2g_file,
        bfh  
    );
}


int main(int argc, char* argv[]) {
  using namespace clipp ;
  using std::cout ;
  enum class mode {help, simulate} ;

  mode selected = mode::help ;

  std::string alevinDir ;
  std::string t2g_file ;

  auto simulateMode = (
    command("construct").set(selected, mode::simulate),

    (option("-a", "--alevin-dir") &
    value("alevin-dir", alevinDir)) %
    "Directory to contain alevin related files",
    
    (option("-t", "--transcript-gene") &
    value("transcript-gene", t2g_file)) %
    "Directory to contain alevin related files"
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
    std::cout << make_man_page(cli, "constructgraph");
    return 1;
  }

  if(res){
    switch(selected){
      case mode::simulate: constructGraph(alevinDir, t2g_file) ; break ; 
      case mode::help: std::cout << make_man_page(cli, "constructgraph") ; 
    }
  }
  std::exit(0) ;

}