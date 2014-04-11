//#include "../config.h"

#include "MCSMap.h"
#include "MCSCompound.h"
#include "MCSRingDetector.h"
#include "util.h"

#ifdef HAVE_LIBOPENBABEL

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

using namespace OpenBabel;

#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <map>
#include <string>

using namespace std;

namespace FMCS {
    
    const char MCSCompound::Atom::elements[111][3] = { "Ru", "Re", "Rf", "Rg", "Ra", "Rb", "Rn", "Rh", "Be", "Ba", "Bh", "Bi", "Bk", "Br", "H", "P", "Os", "Ge", "Gd", "Ga", "Pr", "Pt", "Pu", "C", "Pb", "Pa", "Pd", "Cd", "Po", "Pm", "Hs", "Ho", "Hf", "Hg", "He", "Md", "Mg", "K", "Mn", "O", "Mt", "S", "W", "Zn", "Eu", "Zr", "Er", "Ni", "No", "Na", "Nb", "Nd", "Ne", "Np", "Fr", "Fe", "Fm", "B", "F", "Sr", "N", "Kr", "Si", "Sn", "Sm", "V", "Sc", "Sb", "Sg", "Se", "Co", "Cm", "Cl", "Ca", "Cf", "Ce", "Xe", "Tm", "Cs", "Cr", "Cu", "La", "Li", "Tl", "Lu", "Lr", "Th", "Ti", "Te", "Tb", "Tc", "Ta", "Yb", "Db", "Dy", "Ds", "At", "I", "U", "Y", "Ac", "Ag", "Ir", "Am", "Al", "As", "Ar", "Au", "Es", "In", "Mo" };
    
    map<string, int> MCSCompound::Atom::atomTypeMap;
    
    bool MCSCompound::Atom::atomTypeMapInitialized = atomTypeMapInit();
    bool MCSCompound::Atom::atomTypeMapInit() {
        for (int i = 0; i < 111; ++i) {
            stringstream symbolStringStream;
            //symbolStringStream.width(3);
            //symbolStringStream.setf(ios::left, ios::adjustfield);
            symbolStringStream << elements[i];
            string atomSymbol = getUpper(symbolStringStream.str());
            atomTypeMap[atomSymbol] = i+1;
        }
        return true;
    }
    
    
    
    MCSCompound::MCSCompound() 
        : bondCount(0), atomCount(0), atoms(NULL), bonds(NULL) {}
    
    MCSCompound::~MCSCompound() {
        if(atoms != NULL) {
            delete[] atoms;
            atoms = NULL;
        }
        if(bonds != NULL) {
            delete[] bonds;
            atoms = NULL;
        }
    }

#ifdef HAVE_LIBOPENBABEL

    MCSCompound::MCSCompound(const MCSCompound& other) 
        : bondCount(0), atomCount(0), atoms(NULL), bonds(NULL), SdfContentString(other.SdfContentString), SmiContentString(other.SmiContentString) {
            
        if (other.atoms != NULL) {
            atoms = new Atom[other.atomCount];
            memcpy(atoms, other.atoms, sizeof(Atom) * other.atomCount);
            atomCount = other.atomCount;
        }
        if (other.bonds != NULL) {
            bonds = new Bond[other.bondCount];
            memcpy(bonds, other.bonds, sizeof(Bond) * other.bondCount);
            bondCount = other.bondCount;
        }
    }

#else

    MCSCompound::MCSCompound(const MCSCompound& other)
            : bondCount(0), atomCount(0), atoms(NULL), bonds(NULL), SdfContentString(other.SdfContentString) {

            if (other.atoms != NULL) {
                atoms = new Atom[other.atomCount];
                memcpy(atoms, other.atoms, sizeof(Atom) * other.atomCount);
                atomCount = other.atomCount;
            }
            if (other.bonds != NULL) {
                bonds = new Bond[other.bondCount];
                memcpy(bonds, other.bonds, sizeof(Bond) * other.bondCount);
                bondCount = other.bondCount;
            }
        }

#endif
    
    const MCSCompound& MCSCompound::operator=(const MCSCompound& that) {
        if (this == &that) {
            return *this;
        }
        
        if(atoms != NULL) {
            delete[] atoms;
            atoms = NULL;
        }
        if(bonds != NULL) {
            delete[] bonds;
            bonds = NULL;
        }
        
        bondCount = 0;
        atomCount = 0;
        SdfContentString = that.SdfContentString;

#ifdef HAVE_LIBOPENBABEL
        SmiContentString = that.SmiContentString;
#endif
        
        if (that.atoms != NULL) {
            atoms = new Atom[that.atomCount];
            memcpy(atoms, that.atoms, sizeof(Atom) * that.atomCount);
            atomCount = that.atomCount;
        }
        if (that.bonds != NULL) {
            bonds = new Bond[that.bondCount];
            memcpy(bonds, that.bonds, sizeof(Bond) * that.bondCount);
            bondCount = that.bondCount;
        }
        
        return *this;
    }
#ifdef HAVE_LIBOPENBABEL
    void MCSCompound::read(const string& sdfString, ReadType type) {
        switch (type) {
            case SMI:
                parseSMI(sdfString);
                break;
            case SDF:
                parseSDF(sdfString);
                break;
        }
        
        for (int i = 0; i < bondCount; ++i) {
            atoms[bonds[i].firstAtom].neighborAtoms.push_back(bonds[i].secondAtom);
            atoms[bonds[i].firstAtom].neighborBonds.push_back(&bonds[i]);
            atoms[bonds[i].secondAtom].neighborAtoms.push_back(bonds[i].firstAtom);
            atoms[bonds[i].secondAtom].neighborBonds.push_back(&bonds[i]);
        }
        
    }
#else

    void MCSCompound::read(const std::string& sdfString) {
        parseSDF(sdfString);
        for (int i = 0; i < bondCount; ++i) {
            atoms[bonds[i].firstAtom].neighborAtoms.push_back(bonds[i].secondAtom);
            atoms[bonds[i].firstAtom].neighborBonds.push_back(&bonds[i]);
            atoms[bonds[i].secondAtom].neighborAtoms.push_back(bonds[i].firstAtom);
            atoms[bonds[i].secondAtom].neighborBonds.push_back(&bonds[i]);
        }
        
        MCSRingDetector ringDector(*this);
        ringDector.detect();

    }
    
#endif
    string MCSCompound::subgraph(const size_t* index, size_t indexLength, const string& newCompoundName) const {
        
        stringstream content(this->SdfContentString);
        string res;
        string line;
        vector<string> lines, subLines;
        
        while (getline(content, line)) {
            lines.push_back(line);
        }
        
        subLines.push_back(lines[0] + "_" + newCompoundName);
        subLines.push_back("FMCS substructure");
        
        subLines.push_back("Auto Generated from FMCS");
        subLines.push_back("counts line");
        
        stringstream atomNumss(lines[3].substr(0, 3));
        stringstream bondNumss(lines[3].substr(3, 3));
        int numAtoms, numbonds;
        atomNumss >> numAtoms;
        bondNumss >> numbonds;
        
        MCSMap atomMap;
        for (int i = 0; i < indexLength; ++i) {
            subLines.push_back(lines[4+index[i]]);
            atomMap.push_back(index[i], i);
        }
        int numResBonds = 0;
        for (int i = 0; i < numbonds; ++i) {
            string bondline = lines[4+numAtoms+i];
            stringstream beginAtomss(bondline.substr(0,3)), endAtomss(bondline.substr(3,3));
            int beginAtom, endAtom;
            beginAtomss >> beginAtom;
            endAtomss >> endAtom;
            if (atomMap.containsKey(beginAtom-1) && atomMap.containsKey(endAtom-1)) {
                size_t newBeginAtom = atomMap.getValue(beginAtom-1)+1;
                size_t newEndAtom = atomMap.getValue(endAtom-1)+1;
                stringstream resultStringStream;
                resultStringStream.width(3);
                resultStringStream << newBeginAtom;
                resultStringStream.width(3);
                resultStringStream << newEndAtom << bondline.substr(6);
                subLines.push_back(resultStringStream.str());
                ++numResBonds;
            }
        }
        
        stringstream resLine3ss;
        resLine3ss.width(3);
        resLine3ss << indexLength;
        resLine3ss.width(3);
        resLine3ss  << numResBonds << lines[3].substr(6);
        subLines[3] = resLine3ss.str();
        
        for (vector<string>::iterator i = subLines.begin(); i != subLines.end(); ++i) {
            res += *i;
            res += "\n";
        }
        res += "M END\n";
        res += "$$$$";
        return res;	
    }
    
    const MCSCompound::Bond& MCSCompound::Atom::getBond(int anotherAtom) const {
        return *neighborBonds[neighborAtoms.where(anotherAtom)];
    }
    
    size_t MCSCompound::getNeighborID(size_t e, size_t me) const {
        size_t other;
        if (bonds[e].firstAtom == me) {
            other = bonds[e].secondAtom;
        } else if (bonds[e].secondAtom == me) {
            other = bonds[e].firstAtom;
        } else {
            other = static_cast<size_t>(-1);
        }
        return other;
    }
    
    MCSList<size_t> MCSCompound::getAtomList() const {
        MCSList<size_t> l;
        for (size_t i = 0; i < atomCount; ++i) {
            l.push_back(i);
        }
        return l;
    }

#ifdef HAVE_LIBOPENBABEL
    void MCSCompound::parseSMI(const string& smiString) {

        stringstream smiss, outss;
        smiss << smiString;
        OBConversion conv(&smiss, &outss);
        if(conv.SetInAndOutFormats("SMI","SDF")) { 
            OBMol mol;
            conv.Read(&mol);
            conv.Write(&mol);
        }
        parseSDF(outss.str().c_str());
    }
#endif
    
#ifdef HAVE_LIBOPENBABEL
    void MCSCompound::parseSDF(const string& sdfString) {

        stringstream ss;
        stringstream ssSMI;
        stringstream ssSDF;
        
        ss << sdfString;
        ss >> compoundName;
        
        OBConversion conv(&ss, &ssSDF);
        
        if(conv.SetInAndOutFormats("SDF","SDF")) { 
            OBMol mol;
            if(conv.Read(&mol)) {
                mol.DeleteHydrogens();
                atoms = new Atom[mol.NumAtoms()];
                bonds = new Bond[mol.NumBonds()];
                int i = 0;
                FOR_ATOMS_OF_MOL(atom, mol) {
                    atoms[i] = Atom(atom->GetIdx()-1, atom->GetAtomicNum());
                    ++i;
                }
                i = 0;
                FOR_BONDS_OF_MOL(bond, mol) {
                    int bondType = 0;
                    if (bond->IsKSingle()) {
                        bondType = 1;
                    } else if (bond->IsKDouble()) {
                        bondType = 2;
                    } else if (bond->IsKTriple()) {
                        bondType = 3;
                    } else {
                        delete[] atoms;
                        atoms = NULL;
                        delete[] bonds;
                        atoms = NULL;
                        throw InvalidBondTypeException();
                    }
                    bonds[i] = Bond(bond->GetIdx(), bond->GetBeginAtomIdx()-1, bond->GetEndAtomIdx()-1, bondType, bond->IsAromatic(), bond->IsInRing());
                    ++i;
                }
                
                this->bondCount = mol.NumBonds();
                this->atomCount = mol.NumAtoms();
                
                conv.Write(&mol);
                SdfContentString = ssSDF.str();
            }
        }
        
        OBConversion conv2(&ssSDF, &ssSMI);
        if(conv2.SetInAndOutFormats("SDF","SMI")) { 
            OBMol mol;
            conv2.Read(&mol);
            conv2.Write(&mol);
            SmiContentString = ssSMI.str();
        }
    }
#else
    
    string MCSCompound::deleteHydrogens(const string& sdf, vector<size_t>& originalIds) {
    	stringstream originalStringStream;
    	originalStringStream << sdf;
    	string compoundNameLine;
    	string informationLine;
    	string commentLine;
    	getline(originalStringStream, compoundNameLine);
    	getline(originalStringStream, informationLine);
    	getline(originalStringStream, commentLine);

    	string oldCountsLine;
    	getline(originalStringStream, oldCountsLine);

    	string atomCountString = oldCountsLine.substr(0, 3);
    	string bondCountString = oldCountsLine.substr(3, 3);
    	int oldAtomCount = atoi(atomCountString.c_str());
    	int oldBondCount = atoi(bondCountString.c_str());

    	string newAtomBlock;
    	int *newAtomIndex = new int[oldAtomCount];
    	int newAtomCount = 0;
    	for(int i = 0; i < oldAtomCount; ++i) {
    		string atomBlockLine;
    		getline(originalStringStream, atomBlockLine);
            string atomSymbolRawString = atomBlockLine.substr(31, 3);
            stringstream rawStringStream(atomSymbolRawString);
            string atomSymbolString;
            rawStringStream >> atomSymbolString;
    		if (atomSymbolString!= "H") {
    			newAtomBlock += atomBlockLine;
    			newAtomBlock += "\n";
    			newAtomIndex[i] = newAtomCount+1;
                ++newAtomCount;
                originalIds.push_back(i+1);
    		} else {
    			newAtomIndex[i] = -1;
    		}
    	}
    	string newBondBlock;
    	int newBondCount = 0;
    	for (int i = 0; i < oldBondCount; ++i) {
    		string oldBondBlockLine;
    		getline(originalStringStream, oldBondBlockLine);
    		int oldFirstAtomIndex = atoi(oldBondBlockLine.substr(0, 3).c_str());
    		int oldSecondAtomIndex = atoi(oldBondBlockLine.substr(3, 3).c_str());
    		int newFirstAtomIndex = newAtomIndex[oldFirstAtomIndex-1];
    		int newSecondAtomIndex = newAtomIndex[oldSecondAtomIndex-1];
    		if (newFirstAtomIndex != -1 && newSecondAtomIndex != -1) {
    			stringstream newBondBlockLineStringStream;
    			newBondBlockLineStringStream.width(3);
    			newBondBlockLineStringStream << newFirstAtomIndex;
    			newBondBlockLineStringStream.width(3);
    			newBondBlockLineStringStream << newSecondAtomIndex;
    			newBondBlockLineStringStream << oldBondBlockLine.substr(6);
    			newBondBlock += newBondBlockLineStringStream.str();
    			newBondBlock += "\n";
    			++newBondCount;
    		}
    	}
        
    	stringstream newCountLineStringStream;
    	newCountLineStringStream.width(3);
    	newCountLineStringStream << newAtomCount;
    	newCountLineStringStream.width(3);
    	newCountLineStringStream << newBondCount;
    	newCountLineStringStream << oldCountsLine.substr(6);
    	string newCountLine = newCountLineStringStream.str();
    	delete newAtomIndex;
    	newAtomIndex = NULL;

    	string newSDFString = compoundNameLine + "\n"
    			+ informationLine + "\n"
    			+ commentLine + "\n"
    			+ newCountLine + "\n"
    			+ newAtomBlock
    			+ newBondBlock
    			+ "M END\n"
                + "$$$$";

    	return newSDFString;
    }

    void MCSCompound::parseSDF(const string& sdf) {
    	stringstream sdfStringStream;
        vector<size_t> originalIds;
    	SdfContentString = deleteHydrogens(sdf, originalIds);
    	sdfStringStream << SdfContentString;
    	string compoundNameLine;
    	string informationLine;
    	string commentLine;
    	getline(sdfStringStream, compoundNameLine);
    	getline(sdfStringStream, informationLine);
    	getline(sdfStringStream, commentLine);

    	string countsLine;
    	getline(sdfStringStream, countsLine);

    	string atomCountString = countsLine.substr(0, 3);
    	string bondCountString = countsLine.substr(3, 3);

    	atomCount = atoi(atomCountString.c_str());
    	bondCount = atoi(bondCountString.c_str());

    	atoms = new Atom[atomCount];
    	bonds = new Bond[bondCount];

    	for (size_t i = 0; i < atomCount; ++i) {
    		string atomBlockLine;
    		getline(sdfStringStream, atomBlockLine);
            
            string atomSymbolRawString = atomBlockLine.substr(31, 3);
            stringstream rawStringStream(atomSymbolRawString);
            string atomSymbol;
            rawStringStream >> atomSymbol;
            
    		atoms[i] = Atom(i, originalIds[i], MCSCompound::Atom::atomTypeMap[getUpper(atomSymbol)], atomSymbol);
    	}

    	for (size_t i = 0; i < bondCount; ++i) {
    		string bondBlockLine;
    		getline(sdfStringStream, bondBlockLine);
    		int firstAtom = -1, secondAtom = -1, bondType = -1;
    		firstAtom = atoi(bondBlockLine.substr(0, 3).c_str()) - 1;
    		secondAtom = atoi(bondBlockLine.substr(3, 3).c_str()) - 1;
    		bondType = atoi(bondBlockLine.substr(6, 3).c_str());
    		bonds[i] = Bond(i, firstAtom, secondAtom, bondType, false, false);
    	}
    }

#endif

    const MCSCompound::Bond* MCSCompound::getBond(size_t firstAtom, size_t secondAtom) const {
        
        for(int i = 0; i < bondCount; ++i) {
            if ((bonds[i].firstAtom == firstAtom && bonds[i].secondAtom == secondAtom)
                || (bonds[i].firstAtom == secondAtom && bonds[i].secondAtom == firstAtom)) {
                return bonds+i;
            }
        }
        
        return NULL;
    }
	
}
