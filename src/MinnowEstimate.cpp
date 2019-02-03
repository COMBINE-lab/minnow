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

int minnowEstimate(EstimateOptions& eopts){

     // Set up the logger 
    auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>() ;
    auto consoleLog = spdlog::create("minnow-Log", {consoleSink});


    consoleLog->info("Reading reference sequences ...") ;
    consoleLog->info("Loaded estimated options with eq class folder: {}",eopts.eqClassFolder) ;
    
    consoleLog->info("Reference sequence is loaded ...") ;

    return 0 ;
}

