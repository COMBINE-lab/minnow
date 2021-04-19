// Minnow - An efficient simulator of single cell RNA-seq dataset
// This is the main module that supports all the command line parameter
// required for simulating UMI counts, and read sequences. 

#include "CLI/CLI.hpp"
#include "clipp.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <unordered_map>

#include <ghc/filesystem.hpp>
#include "ProgOpts.hpp"
#include "json.hpp"

using json = nlohmann::json;

int minnowSimulate(SimulateOptions& simOpts) ;
int minnowEstimate(EstimateOptions& eopts) ;
int puffIndex(IndexOptions& iopts) ;

void writeCmdParams(std::string& outDir, json& obj){
  if (outDir.back() == '/'){
    outDir.pop_back();
  }

  if (ghc::filesystem::exists(outDir.c_str())) {
    if (!ghc::filesystem::is_directory(outDir.c_str())) {
      std::cerr << outDir << " exists as a file. Cannot create a directory of the same name.";
      std::exit(1);
    }
  } else {
    ghc::filesystem::create_directories(outDir.c_str());
  }

  std::string filename = outDir + "/arg.json" ;
  std::ofstream stream(filename.c_str());
  stream << obj.dump(4) << std::endl ;
}


int main(int argc, char* argv[]) {
  using namespace clipp ;
  using std::cout ;

  IndexOptions indexOpt ;
  SimulateOptions simulateOpt ;
  EstimateOptions estimateOpt ;

  enum class mode {help, index, simulate, estimate} ;

  mode selected = mode::help ;

  auto ensure_file_exists = [](const std::string& s) -> bool {
                              bool exists = ghc::filesystem::exists(s);
                              if (!exists) {
                                std::string e = "The required input file " + s + " does not seem to exist.";
                                throw std::runtime_error{e};
                              }
                              return true;
                            };
  
  auto ensure_input_dir = [](const std::string& inputdir) -> bool {
                              if (ghc::filesystem::is_directory(inputdir)){
                                auto c = inputdir + "/quants_mat_cols.txt";
                                auto r = inputdir + "/quants_mat_rows.txt";
                                if(!ghc::filesystem::exists(c) || !ghc::filesystem::exists(r)){
                                  std::string e = "Either " + c + " or " + r + " does not exist" \
                                                  + "check the " + inputdir ;
                                  throw std::runtime_error{e};
                                }
                              }else{
                                std::string e = "The input directory" + inputdir + "does not exist";
                                throw std::runtime_error{e};
                              }
                              return true ;
  };


  std::vector<std::string> indexWrong;
  auto indexMode = (
    command("index").set(selected, mode::index),
    (required("-r", "--ref") & values(ensure_file_exists, "ref_file", indexOpt.rfile)) % "path to the reference fasta file",
    (required("-o", "--output") & value("output_dir", indexOpt.outdir)) % "directory where index is written",
    (option("-f", "--filt-size") & value("filt_size", indexOpt.filt_size)) % "filter size to pass to TwoPaCo when building the reference dBG",
    (option("--tmpdir") & value("twopaco_tmp_dir", indexOpt.twopaco_tmp_dir)) % "temporary work directory to pass to TwoPaCo when building the reference dBG",
    (option("-k", "--klen") & value("kmer_length", indexOpt.k))  % "length of the k-mer with which the dBG was built (default = 101)",
    (option("-p", "--threads") & value("threads", indexOpt.p))  % "total number of threads to use for building MPHF (default = 16)",
    any_other(indexWrong)
  );

  std::vector<std::string> estimateWrong;
  auto estimateMode = (
    command("estimate").set(selected, mode::estimate),

    // (required("-eq", "--eqClassDir") &
    // value("eqclass_dir", estimateOpt.eqClassFolder)) %
    // "directory containing relevent files produced by the python script",

    (required("-o", "--outdir") &
    value("outdir", estimateOpt.outDir)) %
    "the simulated models will be written",

    (required("-r") &
    value(ensure_file_exists,"reference", estimateOpt.refFile)) %
    "transcript fasta file",

    (option("--ReadLength") & 
    value("Read length", estimateOpt.ReadLength)) % "read length by default is 100",

    (required("--g2t") &
    value(ensure_file_exists,"gene_tr", estimateOpt.gene2txpFile)) %
    "tab separated list of Gene to Transcirpt mapping",

    (required("--bfh") &
    value("bfh_file", estimateOpt.bfhFile)) %
    "BFH file produced by alevin",

    (option("--cluster") &
    value("cluster", estimateOpt.clusterFile)) %
    "Optional cluster file to model cluster based histogram",
    any_other(estimateWrong)

  );


  std::vector<std::string> simulateWrong;
  auto simulateMode = (
    command("simulate").set(selected, mode::simulate),

    // it can either be --splatter-mode or --alevin mode
    // but not both
    required("--alevin-mode").set(simulateOpt.alevinMode, true) | 
    required("--splatter-mode").set(simulateOpt.splatterMode, true) 
                    .if_missing   ( []  { std::cerr << "\n\033[31mNone of the "
                                                    << "--alevin-mode "
                                                    << "--splatter-mode " 
                                                    << "are specified\033[0m \n\n"; } )
                    .if_conflicted( []  { std::cerr << "\n\033[31mUse either of "
                                                    << "two modes"
                                                    << " but not both\033[0m\n\n"; } ),

    // required options
    (required("-i", "--inputdir")
                          .if_missing ([] {std::cerr << "\033[31mrequired option -i/--inputdir missing\033[0m\n\n";})
                          & value(ensure_input_dir,"inputdir", simulateOpt.inputdir))
                          % "directory with matrix file/ if this is a file instead of a dir",
    (required("-r", "--reffile")
                          .if_missing ([] {std::cerr << "\033[31mrequired option -r/--reffile missing\033[0m\n\n";})
                          &  value(ensure_file_exists, "ref_file", simulateOpt.refFile)) % "transcriptome reference file (assumed from fasta file)",
    (required("-w", "--whitelistFile") 
              .if_missing ([] {std::cerr << "\033[31mrequired option -w/--whitelistFile missing\033[0m\n\n";})
              & value(ensure_file_exists, "whitelist_file", simulateOpt.whitelistFile)) % "10X provided cell barcodes, generally named as 737K-august-2016.txt",
    (required("-o", "--outdir") 
              .if_missing ([] {std::cerr << "\033[31m required option -o/--outdir missing \033[0m\n\n";})
              & value("output directory", simulateOpt.outDir)) % "the simulated reads will be written here",
    (required("--t2g") | required("--g2t") 
              .if_missing ([] {std::cerr << "\033[31m required option --t2g missing \033[0m \n\n";})
              & value(ensure_file_exists, "gene_tr", simulateOpt.gene2txpFile)) % "tab separated list of Gene to Transcirpt mapping",


    (option("--metadataDir") & value("meta data folder", simulateOpt.metadataDir)) % "A directory containing metadata files in case the user defined files don't exit",
    (option("--numMolFile") & value("num mol file", simulateOpt.numMolFile)) % "Number of molecules generated from each cell",
    
    (option("--CBLength") & value("Cell barcode length", simulateOpt.CBLength)) % "Cell barcode length by default is 16",
    (option("--UMILength") & value("UMI length length", simulateOpt.UMILength)) % "Cell barcode length by default is 10",
    (option("--ReadLength") & value("Read length", simulateOpt.ReadLength)) % "read length by default is 100",

    // (option("--alevin-mode").set(simulateOpt.alevinMode, true)) %
    // "The program would assume that the input matrix is obtained from Alevin",
    
    // (option("--splatter-mode").set(simulateOpt.splatterMode, true)) %
    // "matrix file is obtained from running splatter",

    (option("--custom").set(simulateOpt.customNames, true)) % "Read custom gene names instead of assigning genes creatively",
    
    // (option("--normal-mode").set(simulateOpt.normalMode, true)) %"user provided matrix",
    
    (option("--testUniqness").set(simulateOpt.testUniqness, true)) % "In splatter mode the gene names are not specified, therefore the genes should are selected according to uniqueness score",
    (option("--reverseUniqness").set(simulateOpt.reverseUniqness, true)) % "Genes are selected with reversed uniqueness scores",
    (option("--useWeibull").set(simulateOpt.useWeibull, true)) % "The transcripts within one gene get the gene counts following weibull distribution",
    
    
    (option("--numOfDoublets") & value("number of Doublets", simulateOpt.numOfDoublets)) % "Number of doublets to be generated",
    
    (option("--gencode").set(simulateOpt.gencode, true)) % "gencode reference has | separator",
    
    (option("--rspd") & value("rspd_dist", simulateOpt.rspdFile)) % "tab separated read start position distribution",
    
    (option("--bfh") & value(ensure_file_exists, "BFH file", simulateOpt.bfhFile)) % "BFH file",
    (option("--gfa") & value(ensure_file_exists, "gfa_file", simulateOpt.gfaFile)) % "gfa file for contigs",
    (option("--dbg").set(simulateOpt.useDBG, true)) % "Use the provided GFA file and BFH",
    (option("--duplicateList") & value(ensure_file_exists, "duplicate list of files", simulateOpt.duplicateFile)) % "list of transcripts that are actually duplicated",
    
    (option("--geneProb") & 
    value("gene level probability", simulateOpt.geneProbFile)) %
    "Gene level probability file (TSV)",
    
    (option("--countProb") & 
    value("global count probability", simulateOpt.countProbFile)) %
    "global scale count probability file",


    // Veclocity related commands 
    (option("--velocity").set(simulateOpt.velocityMode, true)) %
    "In velocity mode we generate reads from out side exons junction",
    
    (option("--spliceVec") & 
    value(ensure_file_exists,"Splicing rate vector", simulateOpt.spliceVecFile)) %
    "splicing rate per cell in a vetor (values should vary from 0 to 1)",
    
    (option("--intronfile") & 
    value(ensure_file_exists,"intron_file", simulateOpt.intronFile)) %
    "Intron bed file which contains the intron boundaries per transcript",
    // end  

    
    (option("--binary").set(simulateOpt.binary, true)) %
    "If the matrix file is written in binary",
    
    // This option has to be implemented properly
    (option("--noDump").set(simulateOpt.noDump, true)) %
    "The output fules will not be dumped",
    

    (option("--uniq") & value("sequence uniqueness file", simulateOpt.uniquenessFile)) % "sequence uniqueness file",
    
    (option("--illum") & 
    value("illumina model", simulateOpt.illuminaModelFile)) %
    "illumina error model file",
    
    //read dedup counts 
    (option("--dupCounts").set(simulateOpt.dupCounts, true)) %
    "Flag for making minnow read the dup counts TSV filtered_cb_frequency.txt in the same folder",
    
    //librarySize provided 
    (option("--librarySize") & 
    value("specify library size", simulateOpt.librarySize))  % 
    "Total library size, if not specified it will be set to the maximum of 100,000",
    
    //Use white list  
    (option("--useWhiteList").set(simulateOpt.useWhiteList, true)) %
    "A switch effective for in --alevin-mode, where the whitelisted cells will be used for simulation",

    //generate Noise, should come with useWhiteList  
    (option("--generateNoisyCells").set(simulateOpt.generateNoisyCells, true)) %
    "Add some noise to the generated count data [OPTION NOT IN USE]",
    
    //polyAsiteFile
    (option("--polyA").set(simulateOpt.samplePolyA, true)) %
    "Flag to sample with polyA sites this should accompany --polyAsite and --polyAfraction",
    
    //polyAsiteFile
    (option("--polyAsite") & 
    value("polyA_site", simulateOpt.polyAsiteFile)) %
    "Fasta file with polyA sites",

    //polyAsiteFile
    (option("--polyAfraction") & 
    value("polyA_site", simulateOpt.polyAsiteFractionFile)) %
    "File with polyA site fraction ",
    
    (option("-s", "--sampleCells") & 
    value("sample_cells", simulateOpt.sampleCells)) %
    "sample this many cells from the set of all cells",
    
    (option("--genome") & 
    value("genome", simulateOpt.genomefile)) %
    "genome FASTA file",
    
    (option("-c", "--numberOfCells") & 
    value("number of PCR cycles", simulateOpt.numOfCells))  % 
    "Number of cells required for simulation (default = 10)",
    
    (option("-g", "--numberOfTranscripts") & 
    value("number of transcripts", simulateOpt.numOfTranscripts))  % 
    "Number of transcripts for simulation (default = 100)",
    
    (option("--clusters") & 
    value("number of transcripts", simulateOpt.clusterFile))  % 
    "Gene cluster file (should be ported with --dbg)",
    
    (option("--PCR") & 
    value("number of PCR cycles", simulateOpt.numOfPCRCycles))  % 
    "Maximum cycles to repeat PCR (default = 5)",
    
    (option("--PCRmodel6").set(simulateOpt.switchOnEffModel, true)) %
    "model6 from Best, Katharine et al. (2015)",
    
    (option("-e", "--error-rate") & 
    value("error rate for sequences", simulateOpt.errorRate))  % 
    "error rate to be introduced while producing sequences",
    
    (option("-p", "--num-threads") & 
    value("number of threads", simulateOpt.numThreads))  % 
    "number of threads to parallelize the process",

    any_other(simulateWrong)
  );

  std::vector<std::string> wrong;
  auto cli = (
    (simulateMode |
     estimateMode |
     indexMode |
     command("--help").set(selected,mode::help) |
     command("-h").set(selected,mode::help) |
     command("help").set(selected,mode::help)
    ),
    option("-v", "--version").call([]{std::cout << "version 0.1.0\n\n";}).doc("show version"),
    any_other(wrong)
  );

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "minnow");
    return 1;
  }

  // debug::print(std::cout, res);
  if(selected != mode::help){
    if(!res.any_error()){
      std::string lastCommand{""};
      std::unordered_map<std::string, std::string> commandMap ;
      json commandMapJson;
      for(const auto& m : res) {
        if(!m.param()->flags().empty()){
          //std::cout << "boolean" << m.index() << ": " << m.arg() << " -> " << m.param()->flags().front();
          commandMap[m.arg()] = m.param()->flags().front(); 
          lastCommand = m.arg();
        }
        
        if(!m.param()->label().empty()){
          commandMap[lastCommand] = m.arg();
        }
      }  
      auto it1 = commandMap.find("-o");
      auto it2 = commandMap.find("--outdir");
      std::string outDir{""}; 
      if( it1 != commandMap.end()){
        outDir = it1->second;
      }else if(it2 != commandMap.end()){
        outDir = it2->second ;
      }else{
        std::cerr << "\033[31m required option -o/--outdir missing \033[0m\n\n" ;
        return 1;
      }

      for(auto it : commandMap){
        if(it.first == it.second){
          commandMapJson[it.first] = "true";
        }else{
          commandMapJson[it.first] = it.second ;
        }
      }
      
      writeCmdParams(outDir, commandMapJson) ;
    }
  }

  if(res){
    switch(selected){
    case mode::index: puffIndex(indexOpt) ; break ;
    case mode::simulate: minnowSimulate(simulateOpt) ; break ;
    case mode::estimate: minnowEstimate(estimateOpt) ; break ;
    case mode::help: std::cout << make_man_page(cli, "minnow") ;
    }
  }else{
    auto b = res.begin() ;
    auto e = res.end() ;
    if(std::distance(b,e) > 0){
      if(b->arg() == "index"){
        std::cout << make_man_page(indexMode, "minnow") ;
      }else if(b->arg() == "simulate"){
        std::cout << usage_lines(simulateMode, "minnow") << "\n" ;
      }else if(b->arg() == "estimate"){
        std::cout << usage_lines(estimateMode, "minnow") << "\n" ;
      }else{
        std::cout << "There is no command \"" << b->arg() << "\"\n" ;
        std::cout << usage_lines(cli, "minnow") << '\n';
      }
    } else{
      std::cout << usage_lines(cli, "minnow") << '\n';
      return 1;
    }
  }

  return 0 ;
}
