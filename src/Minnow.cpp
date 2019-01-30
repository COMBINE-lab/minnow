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

#include "ProgOpts.hpp"

int minnowSimulate(SimulateOptions& simOpts) ;
int minnowEstimate(EstimateOptions& eopts) ;

int main(int argc, char* argv[]) {
  using namespace clipp ;
  using std::cout ;

  SimulateOptions simulateOpt ;
  EstimateOptions estimateOpt ;

  enum class mode {help, simulate, estimate} ;

  mode selected = mode::help ;

  auto estimateMode = (
    command("estimate").set(selected, mode::estimate),

    (required("-eq", "--eqClassDir") & 
    value("eqclass_dir", estimateOpt.eqClassFolder)) %
    "directory containing relevent files produced by the python script",

    (required("-o", "--outdir") & 
    value("mat_file", estimateOpt.outDir)) %
    "the simulated models will be written", 

    (required("--g2t") & 
    value("gene_tr", estimateOpt.gene2txpFile)) %
    "tab separated list of Gene to Transcirpt mapping", 

    (required("--bfh") & 
    value("bfh_file", estimateOpt.bfhFile)) %
    "BFH file produced by alevin", 

    (option("--cluster") & 
    value("cluster", estimateOpt.clusterFile)) %
    "Optional cluster file to model cluster based histogram" 

  );


  auto simulateMode = (
    command("simulate").set(selected, mode::simulate),
    
    (option("--alevin-mode").set(simulateOpt.alevinMode, true)) %
    "The program would assume that the input matrix is obtained from Alevin",
    
    (option("--splatter-mode").set(simulateOpt.splatterMode, true)) %
    "matrix file is obtained from running splatter",
    
    (option("--testUniqness").set(simulateOpt.testUniqness, true)) %
    "matrix file is obtained from running splatter",


    (option("--reverseUniqness").set(simulateOpt.reverseUniqness, true)) %
    "matrix file is obtained from running splatter",

    (option("--useWeibull").set(simulateOpt.useWeibull, true)) %
    "matrix file is obtained from running splatter",
    
    
    (option("--numOfDoublets") & 
    value("number of Doublets", simulateOpt.numOfDoublets)) %
    "Number of doublets to be generated",
    
    (required("-m", "--matfile") & 
    value("mat_file", simulateOpt.matrixFile)) %
    "matrix file is currently obtained from other program",
    
    
    (required("-o", "--outdir") & 
    value("mat_file", simulateOpt.outDir)) %
    "the simulated reads will be written here",
    
    (option("--gencode").set(simulateOpt.gencode, true)) %
    "gencode reference has | separator",
    
    (required("-r", "--reffile") & 
    value("ref_file", simulateOpt.refFile)) %
    "transcriptome reference file (assumed from fasta file)",
    
    (option("--g2t") & 
    value("gene_tr", simulateOpt.gene2txpFile)) %
    "tab separated list of Gene to Transcirpt mapping",
    
    (option("--rspd") & 
    value("rspd_dist", simulateOpt.rspdFile)) %
    "tab separated read start position distribution",
    
    (option("--bfh") & 
    value("BFH gile", simulateOpt.bfhFile)) %
    "BFH file",
    
    (option("--velocity").set(simulateOpt.velocityMode, true)) %
    "In velocity mode we generate reads from exon-exon junction",
    
    (option("--binary").set(simulateOpt.binary, true)) %
    "If the matrix file is written in binary",
    
    (option("--useEqClass").set(simulateOpt.useEqClass, true)) %
    "will use the model file made",
    
    (option("--dbg").set(simulateOpt.useDBG, true)) %
    "Use the provided GFA file and BFH",
    
    (option("--noDump").set(simulateOpt.noDump, true)) %
    "will use the model file made",
    
    (option("--gfa") & 
    value("gfa_file", simulateOpt.gfaFile)) %
    "gfa file for contigs",
    
    (option("--exonLengthFile") & 
    value("exon_length", simulateOpt.exonLengthFile)) %
    "TSV file for transcripts last exon length ",

    //read dedup counts 
    (option("--dupCounts").set(simulateOpt.dupCounts, true)) %
    "Flag for making minnow read the dup counts TSV filtered_cb_frequency.txt in the same folder",
    
    //Use white list  
    (option("--useWhiteList").set(simulateOpt.useWhiteList, true)) %
    "Flag for making minnow read the dup counts TSV filtered_cb_frequency.txt in the same folder",

    //generate Noise, should come with useWhiteList  
    (option("--generateNoisyCells").set(simulateOpt.generateNoisyCells, true)) %
    "Flag for making minnow read the dup counts TSV filtered_cb_frequency.txt in the same folder",
    
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
    
    (option("--intronfile") & 
    value("intron_file", simulateOpt.intronFile)) %
    "Intron bed file which contains the intron boundaries per transcript",
    
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
    "Gene cluster file (should be ported with --useEqClass)",
    
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
    "number of threads to parallelize the process"
  );

  auto cli = (
    (simulateMode | estimateMode |  
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
    std::cout << make_man_page(cli, "minnow");
    return 1;
  }

  if(res){
    switch(selected){
      case mode::simulate: minnowSimulate(simulateOpt) ; break ; 
      case mode::estimate: minnowEstimate(estimateOpt) ; break ;
      case mode::help: std::cout << make_man_page(cli, "minnow") ; 
    }    
  }else{
      std::cout << usage_lines(cli, "minnow") << '\n';
      return 1;
  }
  
}
