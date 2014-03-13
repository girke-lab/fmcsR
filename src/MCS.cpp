#include "../config.h"

#include "MCS.h"
#include "MCSMap.h"
#include "util.h"

#ifdef HAVE_LIBOPENBABEL

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace OpenBabel;

#endif

#include <ctime>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <R.h>

//#ifndef WINDOWS
//#include <signal.h>
//#endif

using namespace std;

namespace FMCS {
    
    bool timeoutStop = false;
    /*
    void process_alarm(int sig) {
        cout << "ATTENTION: time out... calculation may be incomplete.. terminate now.\n";
        cout << endl;
        timeoutStop = true;
    }*/

    MCS::MCS(const MCSCompound& compoundOne, const MCSCompound& compoundTwo,
             size_t userDefinedLowerBound, size_t substructureNumLimit,
             size_t atomMishmatchLower, size_t atomMismatchUpper,
             size_t bondMismatchLower, size_t bondMismatchUpper,
             MatchType matchType, RunningMode runningMode, int timeout)

        : compoundOne(compoundOne.size() > compoundTwo.size() ? compoundTwo : compoundOne),
          compoundTwo(compoundOne.size() > compoundTwo.size() ? compoundOne : compoundTwo),
    	  userDefinedLowerBound(userDefinedLowerBound), substructureNumLimit(substructureNumLimit),
    	  atomMismatchLowerBound(atomMishmatchLower), atomMismatchUpperBound(atomMismatchUpper),
    	  bondMismatchLowerBound(bondMismatchLower), bondMismatchUpperBound(bondMismatchUpper),
    	  matchType(matchType), runningMode(runningMode), _timeout(timeout),
    	  atomMismatchCurr(0), bondMismatchCurr(0), currSubstructureNum(0),
    	  timeUsed(0.0),startTime(0),  bestSize(0), identicalGraph(false), _isTimeout(false),
    	  haveBeenSwapped(compoundOne.size() > compoundTwo.size() ? true : false) {

    		  timeoutStop = false;
              ifstream ruleFile("rules");
              
              string line;
              stringstream ss;
              while (getline(ruleFile, line)) {
                  ss << line;
                  string atom1 = "", atom2 = "";
                  ss >> atom1 >> atom2;
                  if (atom1 != "" && atom2 != "") {
                      int atomType1 = MCSCompound::Atom::atomTypeMap[getUpper(atom1)];
                      int atomType2 = MCSCompound::Atom::atomTypeMap[getUpper(atom2)];
                      if (atomType1 != 0 && atomType2 != 0) {
                          rules[atomType1][atomType2] = true;
                      }
                      
                  }
              }
	}
    
    void MCS::calculate() {
        
        clearResult();
        /*
#ifndef WINDOWS
        signal(SIGALRM, process_alarm);
        alarm(_timeout);
#endif*/
        clock_t start = clock();
		  startTime=start;
#ifdef HAVE_LIBOPENBABEL
        if (compoundOne.getSmiString() == compoundTwo.getSmiString()) {
#else
        if (compoundOne.getSdfString() == compoundTwo.getSdfString()) {
#endif
            identicalGraph = true;
            if (runningMode == DETAIL) {
#ifdef HAVE_LIBOPENBABEL
                smiSet1.push_back(compoundOne.getSmiString());
                smiSet2.push_back(compoundTwo.getSmiString());
#endif
                sdfSet1.push_back(compoundOne.getSdfString());
                sdfSet2.push_back(compoundTwo.getSdfString());
                
                vector<size_t> tmpIdVector; 
                
               // originalIdArray1.push_back(idxOne);
                //originalIdArray2.push_back(idxTwo);
                
                const MCSCompound::Atom* atomsOne = compoundOne.getAtoms();
                size_t atomCountOne = compoundOne.getAtomCount();
                for (int i = 0; i < atomCountOne; ++i) {
                    tmpIdVector.push_back(atomsOne[i].originalId);
                }
                originalIdArray1.push_back(tmpIdVector);
                tmpIdVector.clear();
                
                const MCSCompound::Atom* atomsTwo = compoundTwo.getAtoms();
                size_t atomCountTwo = compoundTwo.getAtomCount();
                for (int i = 0; i < atomCountTwo; ++i) {
                    tmpIdVector.push_back(atomsTwo[i].originalId);
                }
                originalIdArray2.push_back(tmpIdVector);
            }
        } else {
            max();
        }
        
        clock_t end = clock();	
        
        timeUsed = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
        
        if(runningMode == DETAIL) {
            int resultCount = 0;
            for (std::list<MCSMap>::const_iterator iMap = bestList.begin(); iMap != bestList.end(); ++iMap) {
                ++resultCount;
                stringstream resultCountStringStream;
                resultCountStringStream << resultCount;
                string resultCountString = resultCountStringStream.str();
                stringstream resultStringStreamOne(
                		compoundOne.subgraph(iMap->getKeyList(), size(), string("fmcs_") + resultCountString));
                stringstream resultStringStreamTwo(
                		compoundTwo.subgraph(iMap->getValueList(), size(), string("fmcs_") + resultCountString));
                
                
                const size_t* idArrayOnePtr = iMap->getKeyList();
                const size_t* idArrayTwoPtr = iMap->getValueList();
                int length = size();
                
                vector<size_t> idxOne, idxTwo;
                for(int i = 0; i < length; ++i) {
                    idxOne.push_back(compoundOne.atoms[idArrayOnePtr[i]].originalId);
                    idxTwo.push_back(compoundTwo.atoms[idArrayTwoPtr[i]].originalId);
                }
                
                originalIdArray1.push_back(idxOne);
                originalIdArray2.push_back(idxTwo);

#ifdef HAVE_LIBOPENBABEL
                stringstream output1, output2;

                OBConversion conv1(&resultStringStreamOne, &output1);
                OBConversion conv2(&resultStringStreamTwo, &output2);
                
                if(conv1.SetInAndOutFormats("SDF","SMI"))
                {
                    OBMol mol;
                    conv1.Read(&mol);
                    conv1.Write(&mol);
                }
                
                if(conv2.SetInAndOutFormats("SDF","SMI"))
                {
                    OBMol mol;
                    conv2.Read(&mol);
                    conv2.Write(&mol);
                }
                smiSet1.push_back(output1.str());
                smiSet2.push_back(output2.str());
#endif
                sdfSet1.push_back(resultStringStreamOne.str());
                sdfSet2.push_back(resultStringStreamTwo.str());
            }
        }
    }
    
    void MCS::max() {
        MCSList<size_t> atomListOne = compoundOne.getAtomList();
        MCSList<size_t> atomListTwo = compoundTwo.getAtomList();
        grow(atomListOne, atomListTwo);
    }
    
    bool MCS::compatible(size_t atomOne, size_t atomTwo,
    		size_t& bondMisCount, bool& introducedNewComponent) const {
        
        MCSList<size_t> targetNeighborMapping;
        
        const MCSList<size_t>& atomOneNeighborList = compoundOne[atomOne];
        const size_t* atomOneNeighborsPtr = atomOneNeighborList.get();
        size_t atomOneNeighborSize = atomOneNeighborList.size();
        
        for (size_t i = 0; i < atomOneNeighborSize; ++i) {
            if(currentMapping.containsKey(atomOneNeighborsPtr[i])) {
                targetNeighborMapping.push_back(atomOneNeighborsPtr[i]);
            }
        }
        
        MCSList<size_t> currNeighborMapping;
        const MCSList<size_t>& atomTwoNeighborList = compoundTwo[atomTwo];
        const size_t* atomTwoNeighborsPtr = atomTwoNeighborList.get();
        size_t atomTwoNeighborSize = atomTwoNeighborList.size();
        for (size_t i = 0; i < atomTwoNeighborSize; ++i) {
            size_t k = currentMapping.getKey(atomTwoNeighborsPtr[i]);
            if (k != MCSMap::npos) {
                currNeighborMapping.push_back(k);
            }
        }
        if (!targetNeighborMapping.equals(currNeighborMapping)) {
            return false;
        }
        if (targetNeighborMapping.size() == 0) {
            introducedNewComponent = true;
        }
        size_t targetNeighborMappingSize = targetNeighborMapping.size();
        const size_t* targetNeighborMappingPtr = targetNeighborMapping.get();

        if (matchType == DEFAULT) {
            for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
                size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondOne = compoundOne.atoms[atomOne].getBond(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondTwo = compoundTwo.atoms[atomTwo].getBond(n);
                if (bondOne.bondType != bondTwo.bondType) {
                    ++bondMisCount;
                }
            }
        } else if (matchType == AROMATICITY_SENSETIVE){
            for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
                size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondOne = compoundOne.getAtom(atomOne).getBond(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondTwo = compoundTwo.getAtom(atomTwo).getBond(n);
                if (bondOne.isAromatic != bondTwo.isAromatic) {
                    ++bondMisCount;
                } else if (!bondOne.isAromatic && bondOne.bondType != bondTwo.bondType) {
                    ++bondMisCount;
                }
            }
        } else {
            for (size_t i = 0; i < targetNeighborMappingSize; ++i) {
                size_t n = currentMapping.getValue(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondOne = compoundOne.getAtom(atomOne).getBond(targetNeighborMappingPtr[i]);
                const MCSCompound::Bond bondTwo = compoundTwo.getAtom(atomTwo).getBond(n);
                if (bondOne.isInARing != bondTwo.isInARing) {
                    ++bondMisCount;
                } else {
                    if (bondOne.bondType != bondTwo.bondType) {
                        ++bondMisCount;
                    }
                }
            }
        }

        return true;
    }
    
    size_t MCS::top(MCSList<size_t>& atomList) {
        
        size_t bestCandidateAtom = atomList.front();
        size_t candidateAtom = static_cast<size_t>(-1);
        size_t i, bestIdx = 0, candidateIdx;
        size_t atomListSize = atomList.size();
        const size_t* atomPtr = atomList.get();
        for(i = 0; i < atomListSize; ++i) {
            if (compoundOne[atomPtr[i]].size() > compoundOne[bestCandidateAtom].size()) {
                bestCandidateAtom = atomPtr[i];
                bestIdx = i;
            }
            const MCSList<size_t>& neighborAtomList = compoundOne[atomPtr[i]];
            size_t neighborAtomListSize = neighborAtomList.size();
            const size_t* neighborAtomPtr = neighborAtomList.get();
            for (size_t j = 0; j < neighborAtomListSize; ++j) {
                if (currentMapping.containsKey(neighborAtomPtr[j])) {
                    if (candidateAtom == static_cast<size_t>(-1) || compoundOne[atomPtr[i]].size() > compoundOne[candidateAtom].size()) {
                        candidateAtom = atomPtr[i];
                        candidateIdx = i;
                        break;
                    }
                }
            }
        }
        
        if(candidateAtom == static_cast<size_t>(-1)) {
            atomList.eraseIdx(bestIdx);
            return bestCandidateAtom;
        }
		
        if (candidateAtom != static_cast<size_t>(-1)) {
            atomList.eraseIdx(candidateIdx);
        }
        return candidateAtom;
    }
    
    void MCS::boundary() {
        double diff = (double)(clock() - startTime) / CLOCKS_PER_SEC * 1000 ;
		  //printf("%f  at boundary. timeout: %d\n",diff,_timeout);
		  if(!timeoutStop && _timeout != 0 && diff >= _timeout){
			  warning("FMCS did not complete, timeout of %dms exceeded\n",_timeout);
			  timeoutStop = true;
		  }

        if (runningMode == FAST) {
            if (currentMapping.size() > bestSize) {
                
                if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
                    return;
                }
                bestSize = currentMapping.size();
            }
        } else {
            if (currentMapping.size() == size() ) {
                if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
                    return;
                }
                bestList.push_back(currentMapping);
            } else if (currentMapping.size() > size()) {
                
                if (atomMismatchCurr < atomMismatchLowerBound || bondMismatchCurr < bondMismatchLowerBound) {
                    return;
                }
                bestList.clear();
                bestList.push_back(currentMapping);
            }
        }
    }
    
    void MCS::grow(MCSList<size_t>& atomListOne, MCSList<size_t>& atomListTwo) {
#ifndef WINDOWS
        if (timeoutStop) {
            _isTimeout = true;
            return;
        }
#endif
        MCSList<size_t> atomListOneCopy = atomListOne;
        MCSList<size_t> atomListTwoCopy = atomListTwo;
        MCSList<size_t> atomListOneDegrees;
        MCSList<size_t> atomListTwoDegrees;
        
        size_t atomListOneSize = atomListOne.size();
        const size_t* atomListOnePtr = atomListOne.get();
        for (size_t i = 0; i < atomListOneSize; ++i) {
            if (!currentMapping.containsKey(atomListOnePtr[i])) {
                int degree = 0;
                const MCSList<size_t>& neighborAtomList = compoundOne.atoms[atomListOnePtr[i]].neighborAtoms;
                size_t neighborAtomListSize = neighborAtomList.size();
                const size_t* neighborAtomPtr = neighborAtomList.get();
                for (size_t j = 0; j < neighborAtomListSize; ++j) {
                    if (currentMapping.containsKey(neighborAtomPtr[j])) {
                        ++degree;
                    }
                }
                atomListOneDegrees.push_back(degree);
            }
        }
        
        size_t atomListTwoSize = atomListTwo.size();
        const size_t* atomListTwoPtr = atomListTwo.get();
        for (size_t i = 0; i < atomListTwoSize; ++i) {
            if (!currentMapping.containsValue(atomListTwoPtr[i])) {
                int degree = 0;
                const MCSList<size_t>& neighborAtomList = compoundTwo.atoms[atomListTwoPtr[i]].neighborAtoms;
                size_t neighborAtomListSize = neighborAtomList.size();
                const size_t* neighborAtomPtr = neighborAtomList.get();
                for (size_t j = 0; j < neighborAtomListSize; ++j) {
                    if (currentMapping.containsValue(neighborAtomPtr[j])) {
                        ++degree;
                    }
                }
                atomListTwoDegrees.push_back(degree);
            }
        }
        
        size_t currentBound = currentMapping.size();
        size_t atomListOneDegreesSize = atomListOneDegrees.size();
        const size_t* atomListOneDegreesPtr = atomListOneDegrees.get();
        for (size_t i = 0; i < atomListOneDegreesSize; ++i) {
            if (atomListTwoDegrees.contains(atomListOneDegreesPtr[i])) {
                ++currentBound;
                atomListTwoDegrees.erase(atomListOneDegreesPtr[i]);
            }
        }
        
        if (runningMode == FAST) {
            if(currentBound < userDefinedLowerBound || currentBound <= bestSize) {
                return;
            }
        } else {
            
            if(currentBound < userDefinedLowerBound || currentBound < size()) {
                return;
            }
        }
        
        
        while(true) {
            
            if (atomListOneCopy.empty() || atomListTwoCopy.empty()) { // atomListTwoCopy
                boundary();
                return;
            }
            
            size_t topCandidateAtom = top(atomListOneCopy);
            size_t atomListTwoSize = atomListTwoCopy.size(); // atomListTwoCopy
            const size_t* atomListTwoPtr = atomListTwoCopy.get(); // atomListTwoCopy
            for (size_t i = 0; i < atomListTwoSize; ++i) {
                
                bool atomMismatched = false;
                bool atomMismatchAllowed = true;
                int atom1 = compoundOne.getAtom(topCandidateAtom).atomType;
                int atom2 = compoundTwo.getAtom(atomListTwoPtr[i]).atomType;
                if (atom1 != atom2) {
                    
                    if (rules.count(atom1) > 0) {
                        if (!rules[atom1][atom2]) {
                            atomMismatchAllowed = false;
                        }
                    } else if (rules.count(atom2) > 0) {
                        if (!rules[atom2][atom1]) {
                            atomMismatchAllowed = false;
                        }
                    }
                    ++atomMismatchCurr;
                    atomMismatched = true;
                }
                
                if (!(atomMismatchCurr > atomMismatchUpperBound) && atomMismatchAllowed) {
                    
                    size_t bondMisCount = 0;
                    bool introducedNewComponent = false;
                    if (compatible(topCandidateAtom, atomListTwoPtr[i], bondMisCount, introducedNewComponent)) {
                        
                        if (!(bondMismatchCurr + bondMisCount > bondMismatchUpperBound)) {
                            
                            bondMismatchCurr = bondMismatchCurr + bondMisCount;
                            
                            if (introducedNewComponent) {
                                ++currSubstructureNum;
                            }
                            
                            if (!(currSubstructureNum > substructureNumLimit)) {
                                
                                currentMapping.push_back(topCandidateAtom, atomListTwoPtr[i]);
                                
                                atomListTwo.erase(atomListTwoPtr[i]);
                                
                                grow(atomListOneCopy, atomListTwo);
                                
                                atomListTwo.push_back(atomListTwoPtr[i]);
                                currentMapping.pop_back();
                            }
                            if (introducedNewComponent) {
                                --currSubstructureNum;
                            }
                            
                            bondMismatchCurr = bondMismatchCurr - bondMisCount;
                        }
                    }
                    
                }
                if (atomMismatched) {
                    --atomMismatchCurr;
                }
            }
        }
    }
    
    void MCS::clearResult() {
        
        bestSize = 0;
        
        bestList.clear();
        
        identicalGraph = false;
        currentMapping.clear();

#ifdef HAVE_LIBOPENBABEL
        smiSet1.clear();
        smiSet2.clear();
#endif
        sdfSet1.clear();
        sdfSet2.clear();
        
        timeoutStop = false;
        _isTimeout = false;
    }
    
}
