#include <unordered_map>
#include <memory>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "ghc/filesystem.hpp"


#include <unordered_set>
#include <numeric>
#include "ProgOpts.hpp"

uint64_t getNumDistinctKmers(unsigned kmlen, const std::string& ifile);
int fixFastaMain(std::vector<std::string>& args,
                 std::vector<uint32_t>& refIdExtension,
                 std::vector<std::pair<std::string, uint16_t>>& shortRefsNameLen);
int buildGraphMain(std::vector<std::string>& args);
int dumpGraphMain(std::vector<std::string>& args);

int puffIndex(IndexOptions& indexOpts){

  // Set up the logger
  auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>() ;
  auto console = spdlog::create("minnow-index-Log", {consoleSink});

  uint32_t k = indexOpts.k ;
  std::vector<std::string> rfiles = indexOpts.rfile ;
  std::string rfile ;
  std::string outdir = indexOpts.outdir ;
  std::string gfafile = indexOpts.outdir ;

  std::vector<uint32_t> refIdExtensions;
  std::vector<std::pair<std::string, uint16_t>> shortRefsNameLen;

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};

  // Create the output directory
  if (outdir.back() == '/') {
    outdir.pop_back();
  }
  if (ghc::filesystem::exists(outdir.c_str())) {
    if (!ghc::filesystem::is_directory(outdir.c_str())) {
      console->error("{} exists as a file. Cannot create a directory of the same name.", outdir.c_str());
      std::exit(1);
    }
  } else {
    ghc::filesystem::create_directories(outdir.c_str());
  }


  {
    // Running fixfasta
    console->info("Running Fixfasta") ;
    std::vector<std::string> args ;
    args.push_back("--klen");
    args.push_back(std::to_string(k));
    args.push_back("--input");
    args.insert(args.end(), rfiles.begin(), rfiles.end());
    args.push_back("--output");
    args.push_back(outdir+"/ref_k"+std::to_string(k)+"_fixed.fa");

    int ffres = fixFastaMain(args, refIdExtensions, shortRefsNameLen);
    if (ffres != 0) {
      console->error("The fixFasta phase failed with exit code {}", ffres);
      std::exit(ffres);
    }
    // replacing rfile with the new fixed fasta file
    rfile = outdir+"/ref_k"+std::to_string(k)+"_fixed.fa";
  }

  // If the filter size isn't set by the user, estimate it with ntCard
  if (indexOpts.filt_size == -1){
    console->info("Filter size not provided; estimating from number of distinct k-mers");
    auto nk = getNumDistinctKmers(k, rfile);
    double p = 0.001;
    double k = 5.0;
    double logp_k = std::log(p) / k;
    double r = (-k) / std::log(1.0 - std::exp(logp_k));
    indexOpts.filt_size = static_cast<int32_t>(std::ceil(std::log2(std::ceil(nk * r))));
    console->info("ntHll estimated {} distinct k-mers, setting filter size to 2^{}", nk, indexOpts.filt_size);
    /*
      lgp_k = l($p) / $k
      r = (-$k) / l(1 - e(lgp_k))
      ceil(l(ceil($n * r))/l(2))
    */
  }

  {

    console->info("running twopaco") ;
    std::vector<std::string> args;
    args.push_back("twopaco");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-t");
    args.push_back(std::to_string(indexOpts.p));
    args.push_back("-f");
    args.push_back(std::to_string(indexOpts.filt_size));
    args.push_back("--outfile");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("--tmpdir");

    std::string twopaco_tmp_path = indexOpts.twopaco_tmp_dir;
    // if the tmp path wasn't set, then use a subdirectory 
    // of the index directory (that we will later remove).
    if (twopaco_tmp_path.empty()) {
      twopaco_tmp_path = outdir + "/twopaco_tmp";
    }

    // create the tmp directory if we need to (and can). Complain and exit
    // if the user passed an existing file as the target path. 
    if (ghc::filesystem::exists(twopaco_tmp_path.c_str())) {
        if (!ghc::filesystem::is_directory(twopaco_tmp_path.c_str())) {
            console->error("{} exists as a file. Cannot create a directory of the same name.", twopaco_tmp_path.c_str());
            console->flush();
            std::exit(1);
        }
    } else {
        ghc::filesystem::create_directories(twopaco_tmp_path.c_str());
    }
    args.push_back(twopaco_tmp_path);

    args.push_back(rfile);
    buildGraphMain(args);

    // cleanup tmp
    ghc::filesystem::remove_all(twopaco_tmp_path);
  }

  {
    console->info("Running graphdump") ;
    std::vector<std::string> args;
    args.push_back("graphdump");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-s");
    args.push_back(rfile);
    args.push_back("-f");
    args.push_back("pufferized");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("-p");
    args.push_back(outdir);
    dumpGraphMain(args);

    // cleanup what we no longer need
    ghc::filesystem::path outpath{outdir};
    ghc::filesystem::path tmpDBG = outdir / ghc::filesystem::path{"tmp_dbg.bin"};
    if (ghc::filesystem::exists(tmpDBG)) {
      ghc::filesystem::remove(tmpDBG);
    }
  }

  return 0 ;

}
