#include <unordered_map>
#include <memory>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"


#include <unordered_set>
#include <numeric>

#include "MinnowUtil.hpp"
#include "BFHClass.hpp"
#include "ProgOpts.hpp"
#include "ReferenceInfo.hpp"
#include "ghc/filesystem.hpp"

int minnowEstimate(EstimateOptions& eopts){

     // Set up the logger 
    auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>() ;
    auto consoleLog = spdlog::create("minnow-estimate-Log", {consoleSink});


    auto refFileName = eopts.refFile;
    auto gene2txpFile = eopts.gene2txpFile;
    auto outDir = eopts.outDir;
    auto bfhFile = eopts.bfhFile;
    auto readLength = eopts.ReadLength;


    consoleLog->info("Reading reference sequences ...") ;
    Reference refInfo(
        refFileName,
        gene2txpFile,
        readLength,
        consoleLog
    ) ;
    refInfo.updateGene2TxpMap();
    consoleLog->info("Reference sequence is loaded ...") ;

    // Create output directory
    if (outDir.back() == '/'){
        outDir.pop_back();
    }

    if (ghc::filesystem::exists(outDir.c_str())) {
      if (!ghc::filesystem::is_directory(outDir.c_str())) {
        consoleLog->error("{} exists as a file. Cannot create a directory of the same name.", outDir.c_str());
        std::exit(1);
      }
    } else {
      ghc::filesystem::create_directories(outDir.c_str());
    }
    
    consoleLog->info("Loaded estimated options with eq class folder: {}",eopts.eqClassFolder) ;
    
    BFHClass* eqClassPtr = new BFHClass(consoleLog) ;
    if(bfhFile != ""){
        // NOTE: This branch is currently
        // active
        eqClassPtr->loadBFHLite(
            bfhFile,
            refInfo,
            outDir
        ) ;    
    }


    consoleLog->info("Wrote the estimated files") ;

    return 0 ;
}

