link_directories(
	${GAT_SOURCE_DIR}/external/install/lib
)

set(
    GRAPH_LIB_SRC
    Graph.cpp
    MinnowFS.cpp 
    TranscriptGroup.cpp
    xxhash.c
)

set(
    VALIDATE_LIB_SRC
    FASTAParser.cpp
    FastxParser.cpp
    edlib.cpp
    GFAReader.cpp
)
set(
    FIXFASTA_LIB_SRC
    FASTAParser.cpp
    FastxParser.cpp
    xxhash.c
)

set(minnow_lib_srcs
    MinnowSimulate.cpp
    MinnowEstimate.cpp
    GFAReader.cpp
    MatrixParser.cpp
    xxhash.c
    FASTAParser.cpp
    FastxParser.cpp
    MinnowFS.cpp
    MinnowUtil.cpp
    BFHClass.cpp
    PCR.cpp
)


add_library(minnowlib STATIC ${minnow_lib_srcs})
add_executable(minnow Minnow.cpp)
target_link_libraries(minnow minnowlib z twopaco graphdump)

add_library(graphlib STATIC ${GRAPH_LIB_SRC})
add_executable(constructgraph UMIGraphConstructor.cpp)
target_link_libraries(constructgraph graphlib)

add_library(validatelib STATIC ${VALIDATE_LIB_SRC})
add_executable(validate MinnowValidator.cpp)
target_link_libraries(validate validatelib z)

add_library(fixfastalib STATIC ${FIXFASTA_LIB_SRC})
add_executable(fixfasta FixFasta.cpp)
target_link_libraries(fixfasta fixfastalib z)

#add_executable(gfftest TestGFFParser.cpp)
#target_link_libraries(gfftest gff)


#add_dependencies(minnow libgff)
#add_dependencies(gfftest libgff)
