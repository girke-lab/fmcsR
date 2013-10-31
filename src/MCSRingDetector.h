

#ifndef _MCSRINGDETECTOR_H
#define _MCSRINGDETECTOR_H

#include "../config.h"

#ifndef HAVE_LIBOPENBABEL

#include <set>
#include <vector>
#include <cmath>
#include <map>
#include <utility>

#include "MCSCompound.h"

namespace FMCS {
    
    class MCSRingDetector {
        struct Vertex {
            std::vector<int> connectedEdges;
            Vertex() {}
            Vertex(std::vector<int> connectedEdges) : connectedEdges(connectedEdges) {}
            friend bool operator<(const Vertex& left, const Vertex& right) {
                return left.connectedEdges.size() < right.connectedEdges.size();
            }
            void removeConnectedEdge(int edgeId) {
                std::vector<int>::iterator iter;
                for (iter = connectedEdges.begin(); iter != connectedEdges.end(); ++iter) {
                    if (*iter == edgeId) {
                        break;
                    }
                }
                if (iter != connectedEdges.end()) {
                    connectedEdges.erase(iter);
                }
            }
        };
        
        struct Edge {
            std::vector<int> vertexPath;
            std::vector<int> edgePath;
            
            Edge() {}
            Edge(const std::vector<int>& vertexPath, const std::vector<int>& edgePath) 
            : vertexPath(vertexPath), edgePath(edgePath) {}
            
            const std::vector<int>& get() const {
                return vertexPath;
            }
            
            int front() const {
                return vertexPath.front();
            }
            
            int back() const {
                return vertexPath.back();
            }
            
            size_t size() const {
                return vertexPath.size();
            }
        };
        
        struct Ring {
            
            static bool electronMapInitialized;
            static std::map<std::string, int> electronMap;
            std::vector<int> vertexPath;
            std::vector<int> edgePath;
            
            std::map<int, int> vertexIndexMap;
            
            const MCSCompound* compoundPtr;
            
            Ring() : compoundPtr(NULL) {}
            Ring(const Edge& edge, const MCSCompound* compoundPtr) : compoundPtr(compoundPtr) {
                if (edge.front() != edge.back()) {
                    return;
                }
                this->edgePath = edge.edgePath;
                this->vertexPath = edge.vertexPath;
                this->vertexPath.pop_back();
                
                for (int i = 0; i < this->vertexPath.size(); ++i) {
                    vertexIndexMap[vertexPath[i]] = i;
                }
            }
            
            ~Ring() {
                compoundPtr = NULL;
            }
            
            static bool electronMapInit() {
                electronMap.insert(std::pair<std::string, int>("Al", 3));
                electronMap.insert(std::pair<std::string, int>("As", 5));
                electronMap.insert(std::pair<std::string, int>("At", 7));
                electronMap.insert(std::pair<std::string, int>("B", 3));
                electronMap.insert(std::pair<std::string, int>("Bi", 5));
                electronMap.insert(std::pair<std::string, int>("Br", 7));
                electronMap.insert(std::pair<std::string, int>("C", 4));
                electronMap.insert(std::pair<std::string, int>("Cl", 7));
                electronMap.insert(std::pair<std::string, int>("F", 7));
                electronMap.insert(std::pair<std::string, int>("Ga", 3));
                electronMap.insert(std::pair<std::string, int>("Ge", 4));
                electronMap.insert(std::pair<std::string, int>("I", 7));
                electronMap.insert(std::pair<std::string, int>("In", 3));
                electronMap.insert(std::pair<std::string, int>("N", 5));
                electronMap.insert(std::pair<std::string, int>("O", 6));
                electronMap.insert(std::pair<std::string, int>("P", 5));
                electronMap.insert(std::pair<std::string, int>("Pb", 4));
                electronMap.insert(std::pair<std::string, int>("Po", 6));
                electronMap.insert(std::pair<std::string, int>("S", 6));
                electronMap.insert(std::pair<std::string, int>("Sb", 5));
                electronMap.insert(std::pair<std::string, int>("Se", 6));
                electronMap.insert(std::pair<std::string, int>("Si", 4));
                electronMap.insert(std::pair<std::string, int>("Sn", 4));
                electronMap.insert(std::pair<std::string, int>("Te", 6));
                electronMap.insert(std::pair<std::string, int>("Tl", 3));
                
                return true;
            }
            
            
            bool isSp2Hybridized(size_t vertex, int level, bool& hasLonePair) const;
            bool isAromatic() const;
            
            int leftVertex(size_t vertex) const {
                int vertexIndex = vertexIndexMap.find(vertex)->second;
                
                if (vertexIndex > 0) {
                    return vertexPath[vertexIndex-1];
                } else {
                    return vertexPath[vertexPath.size()-1];
                }
            }
            
            int rightVertex(size_t vertex) const {
                int vertexIndex = vertexIndexMap.find(vertex)->second;
                
                if (vertexIndex < vertexPath.size()-1) {
                    return vertexPath[vertexIndex + 1];
                } else {
                    return vertexPath[0];
                }
            }
            
            int leftEdge(size_t vertex) const {
                int vertexIndex = vertexIndexMap.find(vertex)->second;
                if (vertexIndex > 0) {
                    return edgePath[vertexIndex-1];
                } else {
                    return edgePath[edgePath.size()-1];
                }
            }
            
            int rightEdge(size_t vertex) const {
                int vertexIndex = vertexIndexMap.find(vertex)->second;
                return edgePath[vertexIndex];
            }
            
        };
        
        int vertexId;
        int edgeId;
        
        MCSCompound& compound;
        std::map<int, Vertex> vertexMap;
        std::map<int, Edge> edgeMap;
        
        std::vector<int> vertexQueue;
        std::vector<Ring> rings;
        
        
        void addEdge(const Edge& newEdge);
        Edge catEdge(const Edge& one, const Edge& another);
        bool canCat(const Edge& one, const Edge& another);

        void sortVertexQueue();
        void remove(int vertex);
        void convert();
        
        int nextVertexId() { return ++vertexId; }
        int nextEdgeId() { return ++edgeId; }
        
        int currVertexId() { return vertexId; }
        int currEdgeId() { return edgeId; }
    public:

        MCSRingDetector(MCSCompound& compound) : compound(compound), vertexId(-1), edgeId(-1) {
            convert();
        }
        void detect();
        
    };
    
}

#endif 

#endif // _MCSRINGDETECTOR_H
