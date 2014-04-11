
#ifndef _MCSCOMPOUND_H
#define _MCSCOMPOUND_H

#include "../config.h"

#include <iostream>
#include <string>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "MCSList.h"
#include "util.h"

namespace FMCS {
    
    class MCSCompound {
        friend class MCS;
        friend class MCSRingDetector;
        
        struct Bond;

        struct Atom {

            MCSList<size_t> neighborAtoms;
            MCSList<Bond*> neighborBonds;
#ifdef HAVE_LIBOPENBABEL
            unsigned int atomType;
#else
            int atomType;
            std::string atomSymbol; 
#endif

            size_t atomId;
            size_t originalId;

            Atom() : atomType(0), atomId(-1), originalId(-1) {}

#ifdef HAVE_LIBOPENBABEL
            Atom(size_t id, unsigned int type) 
                : atomId(id), atomType(type) {}
#else
            Atom(size_t id, size_t originalId, int atomType, std::string atomSymbol) 
                : atomId(id), originalId(originalId), atomType(atomType), atomSymbol(atomSymbol) {}
#endif

            size_t degree() const { return neighborAtoms.size(); }

            const Bond& getBond(int v2) const;
            static const char elements[111][3];
            static const char code[111];
            static std::map<std::string, int> atomTypeMap;
            static bool atomTypeMapInitialized;
            static bool atomTypeMapInit();

        };
        
        struct Bond {

            size_t bondId;
            
            size_t firstAtom;
            size_t secondAtom;
            size_t bondType;

            bool isAromatic;
            bool isInARing;

            Bond() : bondId(-1), firstAtom(-1), secondAtom(-1), bondType(0), isAromatic(false), isInARing(false) {}

            Bond(size_t bondId, size_t firstAtom, size_t secondAtom,
            		size_t bondType, bool isAromatic, bool isInARing)
                : bondId(bondId), firstAtom(firstAtom),
                  secondAtom(secondAtom), bondType(bondType),
                  isAromatic(isAromatic), isInARing(isInARing) {}

            bool isSingleBond() const { return bondType == 1; }
            bool isDoubleBond() const { return bondType == 2; }
            bool isTripleBond() const { return bondType == 3; }
            
            bool isAromaticBond() const { return isAromatic; }
            bool isRingBond() const { return isInARing; }
        };

        std::string SdfContentString; 

#ifdef HAVE_LIBOPENBABEL
        std::string SmiContentString;
#endif
        
        size_t bondCount;
        size_t atomCount;
        
        Atom *atoms;
        Bond *bonds;
#ifdef HAVE_LIBOPENBABEL
        void parseSMI(const std::string& sdfString);
#endif
        void parseSDF(const std::string& sdfString);
        
        std::string compoundName;
        
    public:
#ifdef HAVE_LIBOPENBABEL
        enum ReadType { SMI, SDF };
#endif
        MCSCompound();
        ~MCSCompound();
        
        MCSCompound(const MCSCompound&);
        const MCSCompound& operator=(const MCSCompound&);
#ifdef HAVE_LIBOPENBABEL
        void read(const std::string& sdfString, ReadType type);
#else
        void read(const std::string& sdfString);
#endif
        
        const std::string& getSdfString() const { return SdfContentString; }
#ifdef HAVE_LIBOPENBABEL
        const std::string& getSmiString() const { return SmiContentString; }
#endif
        
        std::string subgraph(const size_t* index, size_t indexLength, const std::string& newCompoundName) const;
#ifndef HAVE_LIBOPENBABEL
        std::string deleteHydrogens(const std::string& sdf, std::vector<size_t>& originalIds);
#endif
        const Bond* getBond(size_t firstAtom, size_t secondAtom) const;

        const Atom* getAtoms() const { return atoms; }
        const Bond* getBonds() const { return bonds; }
        
        void setRingBond(size_t bondPos) {
            bonds[bondPos].isInARing = true;
        }
        
        void setAromaticBond(size_t bondPos) {
            bonds[bondPos].isAromatic = true;
        }
        
        size_t getAtomCount() const { return atomCount; }
        size_t getBondCount() const { return bondCount; }
        
        inline size_t size() const {
            return atomCount;
        }
        
        std::string getCompoundName() const {
            return compoundName;
        }
        
        inline const Atom& getAtom(size_t atomPos) const {
        	return atoms[atomPos];
        }
        
        size_t getNeighborID(size_t bond, size_t atom) const;
        
        MCSList<size_t> getAtomList() const;
        
        inline const MCSList<size_t>& operator[](size_t atom) const {
            return atoms[atom].neighborAtoms;
        }
    };
}
#endif // _MCSCOMPOUND_H
