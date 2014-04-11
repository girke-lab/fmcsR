#ifndef _MCSLIST_H
#define _MCSLIST_H

#include "../config.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>

namespace FMCS {
    
    template<class T>
    class MCSList {
    private:
        T* data;
        size_t length;
        size_t maxLength;
        
        static const size_t initLength = 30;
        static const size_t growMultiplier = 5;
        static const size_t lengthLimit = 1000;
        
        void grow();
        
    public:
        
        MCSList();
        ~MCSList();
        MCSList(const MCSList<T>& l);
        const MCSList<T>& operator=(const MCSList<T>& l);
        
        void push_back(const T&);
        void pop_back();
        void clear();
        size_t size() const { return length; }
        
        const T* get() const { return data; };
        
        bool contains(const T&) const;
        
        bool empty() const;
        T front() const;
        T back() const;
        size_t erase(const T& d);
        void eraseIdx(size_t idx);
        
        size_t where(const T&) const;
        T& operator[](size_t) const;
        
        bool equals(const MCSList& l) const;
        
        static const size_t	npos = static_cast<size_t>(-1);
        
    };

    template<class T>
    MCSList<T>::MCSList() : data(NULL), length(0), maxLength(0) {}
    
    template<class T>
    MCSList<T>::~MCSList() {
        delete[] data;
        data = NULL;
    }
    
    
    template<class T>
    MCSList<T>::MCSList(const MCSList<T>& m) : data(NULL), length(0), maxLength(initLength) {
        if (m.data == NULL) {
            return;
        }
        clear();
        maxLength = m.maxLength;
        data = new T[m.maxLength];
        memcpy(data, m.data, sizeof(T) * m.length);
        length = m.length;
    }
    
    template<class T>
    const MCSList<T>& MCSList<T>::operator=(const MCSList<T>& m) {
        if (this == &m) {
            return *this;
        }
        
        clear();
        
        if (m.data == NULL) {
            return *this;
        }
        
        maxLength = m.maxLength;
        data = new T[m.maxLength];
        memcpy(data, m.data, sizeof(T) * m.length);
        length = m.length;	
        return *this;
    }
    
    template<class T>
    void MCSList<T>::grow() {
        if (maxLength == lengthLimit) {
            throw std::runtime_error(
                std::string("Length exceeds limit.."));
        }
        
        if(maxLength == 0) {
            maxLength = initLength;
        } else {
            maxLength = maxLength * growMultiplier;
        }
        if (maxLength > lengthLimit) {
            maxLength = lengthLimit;
        }
        
        T* n = new T[maxLength];
        memcpy(n, data, sizeof(T) * length);
        delete[] data;
        data = n;
    }
    
    template<class T>
    void MCSList<T>::push_back(const T& a) {
        if (length >= maxLength) {
            grow();
        }
        data[length] = a;
        ++length;
    }
    
    template<class T>
    void MCSList<T>::pop_back() {
        --length;
    }
    
    template<class T>
    void MCSList<T>::clear() {
        delete[] data;
        data = NULL;
        length = 0;
    }
    
    template<class T>
    bool MCSList<T>::contains(const T& d) const {
        for (int i = 0; i < length; ++i) {
            if (d == data[i]) {
                return true;
            }
        }
        return false;
    }
    
    template<class T>
    bool MCSList<T>::empty() const {
        return length == 0;
    }
    
    template<class T>
    T MCSList<T>::front() const {
        if (length == 0) {
            return T();
        }
        return data[0];
    }
    
    template<class T>
    T MCSList<T>::back() const {
        if (length == 0) {
            return T();
        }
        return data[length-1];
    }
    
    template<class T>
    size_t MCSList<T>::erase(const T& d) {
        for (size_t i = 0; i < length; ++i) {
            if (d == data[i]) {
                data[i] = data[length-1];
                --length;
                return i;
            }
        }
        return npos;
    }
    
    template<class T>
    void MCSList<T>::eraseIdx(size_t idx) {
        data[idx] = data[length-1];
        --length;
    }
    
    template<class T>
    size_t MCSList<T>::where(const T& t) const {
        for (size_t i = 0; i < length; ++i) {
            if (data[i] == t) {
                return i;
            }
        }
        return npos;
    }
    
    template<class T>
    T& MCSList<T>::operator[](size_t idx) const {
        return data[idx];
    }
    
    template<class T>
    bool MCSList<T>::equals(const MCSList<T>& l) const {
        if (length != l.length) {
            return false;
        }
        for (size_t i = 0; i < length; ++i) {
            if (!l.contains(data[i])) {
                return false;
            }
        }
        return true;
    }
    
	
}
#endif // _MCSLIST_H
