#ifndef _MCSMAP_H
#define _MCSMAP_H

#include "../config.h"

#include <iostream>

#include "MCSList.h"

namespace FMCS {
    
    class MCSMap {
    private:
        MCSList<size_t> keyList;
        MCSList<size_t> valueList;
        size_t length;
    public:
        static const size_t npos = static_cast<size_t>(-1);
        MCSMap();
        ~MCSMap();
        MCSMap(const MCSMap& m);
        const MCSMap& operator=(const MCSMap& m);
        
        bool containsKey(size_t key) const;
        bool containsValue(size_t value) const;
        void push_back(size_t, size_t);
        void pop_back();
        
        size_t size() const {
            return length;
        }
        void clear();
        
        size_t getKey(size_t value) const;
        size_t getValue(size_t key) const;
        
        const size_t* getKeyList() const;
        const size_t* getValueList() const;
    };
	
}

#endif // _MCSMAP_H
