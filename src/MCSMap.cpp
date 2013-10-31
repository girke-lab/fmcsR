#include "../config.h"

#include "MCSMap.h"

using namespace std;

namespace FMCS {
    
    MCSMap::MCSMap() : length(0) {}
    
    MCSMap::~MCSMap() {}
    
    MCSMap::MCSMap(const MCSMap& m) : keyList(m.keyList), valueList(m.valueList), length(m.length) {}
    
    const MCSMap& MCSMap::operator=(const MCSMap& m) { 
        
        if (this == &m) {
            return *this;
        }
        
        valueList = m.valueList;
        keyList = m.keyList;
        length = m.length;
        return *this;
    }
    
    bool MCSMap::containsKey(size_t key) const {
        return keyList.contains(key);
    }
    
    bool MCSMap::containsValue(size_t value) const {
        return valueList.contains(value);
    }
    
    void MCSMap::push_back(size_t a, size_t b) {
        keyList.push_back(a);
        valueList.push_back(b);
        length = keyList.size();
    }
    
    void MCSMap::pop_back() {
        keyList.pop_back();
        valueList.pop_back();
        length = keyList.size();
    }
    
    size_t MCSMap::getKey(size_t value) const {
        size_t n = valueList.where(value);
        if (n == MCSList<size_t>::npos) {
            return npos;
        }
        return keyList[n];
    }
    
    size_t MCSMap::getValue(size_t key) const {
        size_t n = keyList.where(key);
        if (n == MCSList<size_t>::npos) {
            return npos;
        }
        return valueList[n];
    }
    
    const size_t* MCSMap::getKeyList() const {
        return keyList.get();
    }
    
    const size_t* MCSMap::getValueList() const {
        return valueList.get();
    }
    
    void MCSMap::clear() {
        keyList.clear();
        valueList.clear();
        length = 0;
    }
    
}
