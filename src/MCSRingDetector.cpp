#include "../config.h"

#ifndef HAVE_LIBOPENBABEL

#include <set>
#include <map>
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>

#include "MCSRingDetector.h"
#include "MCSCompound.h"

using namespace std;

namespace FMCS {
    
    void MCSRingDetector::addEdge(const Edge& newEdge) {
        int edgeId = nextEdgeId();
        edgeMap[edgeId] = newEdge;
        vertexMap[newEdge.front()].connectedEdges.push_back(edgeId);
        vertexMap[newEdge.back()].connectedEdges.push_back(edgeId);
    }
    
    MCSRingDetector::Edge MCSRingDetector::catEdge(const Edge& one, const Edge& another) {
        vector<int> vertexPath;
        vector<int> edgePath;
        
        if (one.back() == another.front()) {
        
            vertexPath.insert(vertexPath.begin(), one.vertexPath.begin(), one.vertexPath.end());
            vertexPath.insert(vertexPath.end(), another.vertexPath.begin()+1, another.vertexPath.end());
            
            edgePath.insert(edgePath.begin(), one.edgePath.begin(), one.edgePath.end());
            edgePath.insert(edgePath.end(), another.edgePath.begin(), another.edgePath.end());
            
            
        } else if (one.front() == another.back()) {
        
            vertexPath.insert(vertexPath.begin(), another.vertexPath.begin(), another.vertexPath.end());
            vertexPath.insert(vertexPath.end(), one.vertexPath.begin()+1, one.vertexPath.end());
            
            edgePath.insert(edgePath.begin(), another.edgePath.begin(), another.edgePath.end());
            edgePath.insert(edgePath.end(), one.edgePath.begin(), one.edgePath.end());
            
        } else if (one.back() == another.back()) {
        
            vertexPath.insert(vertexPath.begin(), one.vertexPath.begin(), one.vertexPath.end());
            vertexPath.insert(vertexPath.end(), another.vertexPath.rbegin()+1, another.vertexPath.rend());
            
            edgePath.insert(edgePath.begin(), one.edgePath.begin(), one.edgePath.end());
            edgePath.insert(edgePath.end(), another.edgePath.rbegin(), another.edgePath.rend());
            
        } else if (one.front() == another.front()) {
        
            vertexPath.insert(vertexPath.begin(), another.vertexPath.rbegin(), another.vertexPath.rend());
            vertexPath.insert(vertexPath.end(), one.vertexPath.begin()+1, one.vertexPath.end());
            
            edgePath.insert(edgePath.begin(), another.edgePath.rbegin(), another.edgePath.rend());
            edgePath.insert(edgePath.end(), one.edgePath.begin(), one.edgePath.end());
            
        }
        return Edge(vertexPath, edgePath);
    }
    
    bool MCSRingDetector::canCat(const Edge& one, const Edge& another) {
        if (!(one.back() == another.back()) 
            && !(one.back() == another.front()) 
            && !(one.front() == another.front())
            && !(one.front() == another.back())) {
            return false;
        }

        if (!(one.vertexPath.size() > 2 && another.vertexPath.size() > 2)) {
            return true;
        }

        set<int> oneVertexSet(one.vertexPath.begin()+1, one.vertexPath.end()-1);
        set<int> anotherVertexSet(another.vertexPath.begin()+1, another.vertexPath.end()-1);
        
        vector<int> intersection;
        set_intersection (oneVertexSet.begin(), oneVertexSet.end(), 
                          anotherVertexSet.begin(), anotherVertexSet.end(), 
                          back_inserter(intersection));
        
        return intersection.size() == 0;
    }
    
    void MCSRingDetector::sortVertexQueue() {
        size_t queueSize = vertexQueue.size();
        for (int i = 0; i < vertexQueue.size(); ++i) {
            for (int j = 0; j < queueSize-1-i; ++j ) {
                if (vertexMap[vertexQueue[j]] < vertexMap[vertexQueue[j+1]]) {
                    int tmp = vertexQueue[j];
                    vertexQueue[j] = vertexQueue[j+1];
                    vertexQueue[j+1] = tmp;
                }
            }
        }
    }
    
    void MCSRingDetector::remove(int vertex) {
        
        int connectedEdgeCount = vertexMap[vertex].connectedEdges.size();
        
        for (int i = 0; i < connectedEdgeCount-1; ++i) {
            for (int j = i + 1; j < connectedEdgeCount; ++j) {
                const Edge& one = edgeMap[vertexMap[vertex].connectedEdges[i]];
                const Edge& another = edgeMap[vertexMap[vertex].connectedEdges[j]];
                if (!canCat(one, another)) {
                    continue;
                }
                
                Edge newEdge = catEdge(one, another);
                if (newEdge.front() == newEdge.back()) {
                    rings.push_back(Ring(newEdge, &compound));
                } else {
                    addEdge(newEdge);
                }
            }
        }
        
        for (int i = 0; i < connectedEdgeCount; ++i) {
            int edgeId = vertexMap[vertex].connectedEdges[i];
            const Edge& edge = edgeMap[edgeId];
            if (edge.front() == vertex) {
                Vertex& otherVertex = vertexMap[edge.back()];
                otherVertex.removeConnectedEdge(edgeId);
            } else {
                Vertex& otherVertex = vertexMap[edge.front()];
                otherVertex.removeConnectedEdge(edgeId);
            }
            
            edgeMap.erase(edgeId);
        }
    }
    
    void MCSRingDetector::convert() {
        const MCSCompound::Atom* atoms = compound.getAtoms();
        int atomCount = compound.getAtomCount();
        const MCSCompound::Bond* bonds = compound.getBonds();
        int bondCount = compound.getBondCount();
        
        for (int i = 0; i < atomCount; ++i) {
            vector<int> connectedEdges;
            size_t degree = atoms[i].neighborBonds.size();
            for (int j = 0; j < degree; ++j) {
                connectedEdges.push_back(atoms[i].neighborBonds[j]->bondId); 
            }
            vertexMap[nextVertexId()] = Vertex(connectedEdges);
            vertexQueue.push_back(currVertexId());
        }
        
        for (int i = 0; i < bondCount; ++i) {
            vector<int> vertexPath;
            vertexPath.push_back(bonds[i].firstAtom);
            vertexPath.push_back(bonds[i].secondAtom);
            
            vector<int> edgePath;
            edgePath.push_back(nextEdgeId());
            
            edgeMap[currEdgeId()] = Edge(vertexPath, edgePath);
        }
    }
    
    void  MCSRingDetector::detect() {
        while (!vertexQueue.empty()) {
            int vertex = vertexQueue.back();
            vertexQueue.pop_back();
            remove(vertex);
            sortVertexQueue();
        }
        int aromaticCount = 0;
        for (vector<Ring>::const_iterator ringIterator = rings.begin(); ringIterator != rings.end(); ++ringIterator) {
            const vector<int>& ringEdges = ringIterator->edgePath;
            for (vector<int>::const_iterator ringEdgeIter = ringEdges.begin(); ringEdgeIter != ringEdges.end(); ++ringEdgeIter) {
                compound.setRingBond(*ringEdgeIter);
            }
            if (ringIterator->isAromatic()) {
                for (vector<int>::const_iterator ringEdgeIter = ringEdges.begin(); ringEdgeIter != ringEdges.end(); ++ringEdgeIter) {
                    compound.setAromaticBond(*ringEdgeIter);
                }
            }
        }
    }
    
    map<string, int> MCSRingDetector::Ring::electronMap;
    bool MCSRingDetector::Ring::electronMapInitialized = electronMapInit();
    
    bool MCSRingDetector::Ring::isSp2Hybridized(size_t vertex, int level, bool& hasLonePair) const {
        
        if (level > vertexPath.size()) {
            return false;
        }
        
        const MCSCompound::Atom* atoms= compoundPtr->getAtoms();
        const MCSCompound::Atom& atom = atoms[vertex];
        
        if (electronMap[atom.atomSymbol] == 0) {
            return false;
        }
        
        int lonePairCount = 0;
        int neighborCount = atom.neighborBonds.size();
        int orbitalCount = 0;
        int electronCount = 0;
        
        const MCSCompound::Bond* const* neighborBonds = atom.neighborBonds.get();
        for (int i = 0; i < neighborCount; ++i) {
            if(neighborBonds[i]->isSingleBond()) {
                electronCount += 1;
            } else if (neighborBonds[i]->isDoubleBond()) {
                electronCount += 2;
            } else {
                electronCount +=3;
            }
        }
        
        //cout << electronCount << endl;
        
        orbitalCount = electronCount;
        
        int emptyOrbitalCount = 4 - orbitalCount;
        int freeElectornCount = electronMap[atom.atomSymbol] - electronCount;
        
        if (freeElectornCount > emptyOrbitalCount) {
            neighborCount += (emptyOrbitalCount*2 - freeElectornCount);
            lonePairCount = freeElectornCount / 2;
        } else {
            neighborCount += freeElectornCount;
        }
        
        if (lonePairCount > 0) {
            hasLonePair = true;
        }
        bool dumm;
        
        if (lonePairCount + neighborCount == 3) {
            return true;
        } else if (lonePairCount > 0) {
            return isSp2Hybridized(leftVertex(vertex), level+1, dumm) || isSp2Hybridized(rightVertex(vertex), level+1, dumm);
        } else {
            return false;
        }
    }
    
    bool MCSRingDetector::Ring::isAromatic() const {
        
        const MCSCompound::Bond* bonds= compoundPtr->getBonds();
        int piElectornCount = 0;
        for (vector<int>::const_iterator i = vertexPath.begin(); i != vertexPath.end(); ++i) {
            bool hasLonePair = false;
            if (! isSp2Hybridized(*i, 1, hasLonePair)) {
                return false;
            }
            
            size_t leftEdge = this->leftEdge(*i);
            size_t rightEdge = this->rightEdge(*i);
            
            if (bonds[leftEdge].isDoubleBond()) {
                ++piElectornCount;
            }
            if (bonds[rightEdge].isDoubleBond()) {
                ++piElectornCount;
            }
            
            if (!bonds[leftEdge].isDoubleBond() && !bonds[rightEdge].isDoubleBond() && hasLonePair) {
                piElectornCount += 2;
            }
        }
        
        return (piElectornCount - 2) % 4 == 0;
    }

}


#endif
