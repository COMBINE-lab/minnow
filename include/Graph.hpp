#ifndef GRAPH_HPP
#define GRAPH_HPP

//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/functional/hash.hpp>

#include <sparsepp/spp.h>
#include <utility>
#include <cstdint>
#include <vector>

#include "xxhash.h"

//#include "AlevinUtils.hpp"
//#include "edlib.h"

//namespace alevin {

  
  namespace graph {

    using VertexT = std::pair<uint32_t, uint32_t>;    

    struct vertexHasher{
      std::size_t operator()(const VertexT& v) const {
        std::size_t seed{0} ;
        spp::hash_combine(seed, v.first);
        spp::hash_combine(seed, v.second);
        return seed ;
      }
    } ;

    enum EdgeType {
      NoEdge,
      BiDirected,
      XToY,
      YToX,
    };

    struct Graph {
      std::vector<VertexT> vertexNames;
      spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t>> edges;

      void add_edge(uint32_t source, uint32_t sink) {
        edges[source].insert(sink);
      }

      size_t num_vertices() {
        return vertexNames.size();
      }

      size_t num_edges() {
        size_t i = 0;
        for(auto& it: edges) {
          i += it.second.size();
        }
        return i;
      }

      uint32_t getEqclassId(uint32_t vertex) {
        return vertexNames[vertex].first;
      }

      spp::sparse_hash_set<uint32_t> getNeighbors(uint32_t vertex) {
        return edges[vertex];
      }

      /*
      // don't need connected component for now
      uint32_t connected_components(std::vector<uint32_t>& component){
        // resize the component based on the number of vertices
        component.resize( vertexNames.size() );

        typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > AdjList;
        AdjList adjList ( vertexNames.size() );

        // iterating over edges and filling the graph
        for (auto& it: edges) {
          uint32_t source = it.first;
          for(uint32_t target: it.second) {
            boost::add_edge(source, target, adjList);
          }
        }

        return boost::connected_components(adjList, component.data());
      }*/
    };

    uint32_t getVertexIndex(spp::sparse_hash_map<VertexT, uint32_t, vertexHasher>& vertMap,
                            VertexT& node);

    EdgeType hasEdge(std::pair<uint64_t, uint32_t> &x,
                     std::pair<uint64_t, uint32_t> &y);
  }
//}

#endif // GRAPH_HPP
