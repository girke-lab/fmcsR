#include "../config.h"

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

#include "MCS.h"

using namespace std;
using namespace FMCS;

extern "C" {

void fmcs_R_wrap(const char** structureStringOne, const char** structureStringTwo, 
    int *atomMismatchLowerBound, int *atomMismatchUpperBound, 
    int *bondMismatchLowerBound, int *bondMismatchUpperBound, 
    int *matchTypeInt, int* runningModeInt, 
    int *timeout, 
    const char** resultIdxOne, const char** resultIdxTwo,
    //const char** resultSdfOne, const char** resultSdfTwo,
    const char** sdfOneSize, const char** sdfTwoSize, const char** mcsSize) {

    if (*structureStringOne == NULL) {
		//cout << "input structure one cannot be NULL...\n";
		return;
	}

	if (*structureStringTwo == NULL) {
		//cout << "input structure two cannot be NULL...\n";
		return;
	}
	
	int substructureNumLimit = 1;
    int userDefinedLowerBound = 0;
	/*
    try {
    */
    	MCS::MatchType matchType;
        switch(*matchTypeInt) {
            case 0: matchType = MCS::DEFAULT; break;
            case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
            case 2: matchType = MCS::RING_SENSETIVE; break;
            default: 
                ;
        }


        MCS::RunningMode runningMode;
        switch(*runningModeInt) {
            case 0: runningMode = MCS::FAST; break;
            case 1: runningMode = MCS::DETAIL; break;
            default:
                ;
        }
        
        MCSCompound compoundOne, compoundTwo;
        
#ifdef HAVE_LIBOPENBABEL
        compoundOne.read(string(*structureStringOne), MCSCompound::SDF);
        compoundTwo.read(string(*structureStringTwo), MCSCompound::SDF);
#else
        compoundOne.read(string(*structureStringOne));
        compoundTwo.read(string(*structureStringTwo));
            
#endif           
        MCS mcs(compoundOne, compoundTwo,
            userDefinedLowerBound, substructureNumLimit,
            *atomMismatchLowerBound, *atomMismatchUpperBound,
            *bondMismatchLowerBound, *bondMismatchUpperBound,
            matchType, runningMode, *timeout);
       /*     
        cout << "Matching mode: ";
        if (matchType == MCS::DEFAULT) {
            cout << "Static\n";
        } else if (matchType == MCS::AROMATICITY_SENSETIVE) {
            cout << "Strict aromatic\n";
        } else if (matchType == MCS::RING_SENSETIVE) {
            cout << "Strict ring\n";
        } else {
            throw InvalidMatchTypeException();
        }


		cout << "Atom Mismatches Lower Bound: " << *atomMismatchLowerBound << "\n";
		cout << "Atom Mismatches Upper Bound: " << *atomMismatchUpperBound << "\n";
		cout << "Bond Mismatches Lower Bound: " << *bondMismatchLowerBound << "\n";
		cout << "Bond Mismatches Upper Bound: " << *bondMismatchUpperBound << "\n";
		
		cout << "\n";*/
/*	
		if (*timeout == 0) {
			cout << "Timeout: " << "turned off \n\n";
		} else {
			cout << "Timeout: " << *timeout << " seconds...\n\n";
        }
*/

/*
#ifdef HAVE_LIBOPENBABEL
        cout << "Srtructure 1: " << mcs.getCompoundOne().getSmiString();
        cout << "Size: " << mcs.getCompoundOne().size() << "\n\n";
        cout << "Srtructure 2: " << mcs.getCompoundTwo().getSmiString();
        cout << "Size: " << mcs.getCompoundTwo().size() << "\n\n";
            
#else
        cout << "Size of Structure 1: " << mcs.getCompoundOne().size() << "\n";
        cout << "Size of Structure 2: " << mcs.getCompoundTwo().size() << "\n";
        cout << endl;
#endif*/
            
        mcs.calculate();
/*
        cout << "FMCS size: " << mcs.size() << endl;
        cout << mcs.getTime() << " microseconds used..." << endl << endl;

        cout << "Tanimoto score: " 
            << (double)mcs.size() / (mcs.getCompoundOne().size() + mcs.getCompoundTwo().size() - mcs.size()) 
            << endl << endl;
        */
        static int cmpOneSize, cmpTwoSize, mSize;
    	cmpOneSize = mcs.getCompoundOne().size();
    	cmpTwoSize = mcs.getCompoundTwo().size();
    	mSize = mcs.size();
        
        if (runningMode == MCS::DETAIL) {
            
            //list<string> sdfList1 = mcs.getFirstSdfResultStringList();
            //list<string> sdfList2 = mcs.getSecondSdfResultStringList();
            
            list<vector<size_t> > index1 = mcs.getFirstOriginalIndice();
            list<vector<size_t> > index2 = mcs.getSecondOriginalIndice();
            
            //cout << index1.size() << " solution(s) found..." << endl;
            /*
            if (index1.size() > 0) {
                for (int i = 0; i < index1.begin()->size(); ++i) {
                    cout << (*index1.begin())[i] << ", ";
                }
                cout << endl;
            }*/
            
            //cout << mcs.getFirstSdfResultStringList().size() << " solution(s) found..." << endl;
            
            stringstream indexOneStringStream, indexTwoStringStream;
            for (list<vector<size_t> >::const_iterator i = index1.begin(); i != index1.end(); ++i) {
            	for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
            		indexOneStringStream << *j << " ";
            	}
            	indexOneStringStream << "\n";
            }
            
            for (list<vector<size_t> >::const_iterator i = index2.begin(); i != index2.end(); ++i) {
            	for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
            		indexTwoStringStream << *j << " ";
            	}
            	indexTwoStringStream << "\n";
            }
            
            static string indexOneString, indexTwoString;
            indexOneString = indexOneStringStream.str();
            indexTwoString = indexTwoStringStream.str();
            
            *resultIdxOne = indexOneString.c_str();
            *resultIdxTwo = indexTwoString.c_str();
            /*
            string sdfs1, sdfs2;
            for (list<string>::const_iterator i = sdfList1.begin(); i != sdfList1.end(); ++i) {
                    sdfs1 += *i;
                    sdfs1 += "\n";
            }
            
            for (list<string>::const_iterator i = sdfList2.begin(); i != sdfList2.end(); ++i) {
                    sdfs2 += *i;
                    sdfs2 += "\n";
            }

            *resultSdfOne = sdfs1.c_str();
            *resultSdfTwo = sdfs2.c_str();
            
            ofstream out("out.txt");
            out << *resultSdfOne << endl;*/
        }
        
        stringstream sizeStringStream;
        
        sizeStringStream << cmpOneSize;
        static string cmpOneSizeString;
        cmpOneSizeString= sizeStringStream.str();
        
        sizeStringStream.str("");
    	sizeStringStream << cmpTwoSize;
        static string cmpTwoSizeString;
        cmpTwoSizeString = sizeStringStream.str();
        
        sizeStringStream.str("");
        sizeStringStream << mSize;
        static string mSizeString;
        mSizeString = sizeStringStream.str();
        
        *sdfOneSize = cmpOneSizeString.c_str();
        *sdfTwoSize = cmpTwoSizeString.c_str();
        *mcsSize = mSizeString.c_str();
/*
	} catch (exception& e) {
		cerr << e.what() << endl;
	}*/
}
/*
void fmcsBatch_R_wrap(const char** query, const char** dataSet, const char** output, 
        int *al, int *au, int *bl, int *bu, int *matchMode, int* timeout) {
        
        const char* queryFile = *query;
        const char* searchSetFile = *dataSet;
        const char* outputFile = *output;
        int atomMismatchLowerBound = *al;
        int atomMismatchUpperBound = *au;
        int boundMismatchLowerBound = *bl;
        int boundMismatchUpperBound = *bu;
        int matchTypeInt = *matchMode;
        
        const int substructureNumLimit = 1;
        const int userDefinedLowerBound = 0;

        MCS::MatchType matchType;
        switch(matchTypeInt) {
            case 0: matchType = MCS::DEFAULT; break;
            case 1: matchType = MCS::AROMATICITY_SENSETIVE; break;
            case 2: matchType = MCS::RING_SENSETIVE; break;
            default: 
                ;
        }

        ifstream dbInstream(searchSetFile);
        ifstream queryInstream(queryFile);
        ofstream resOstream(outputFile);
                
        MCSCompound queryCompound;
        stringbuf sdfbuf;
                
        queryInstream.get(sdfbuf, '$');
        queryInstream.get(sdfbuf, '\n');
        while (sdfbuf.str().rfind("$$$$") ==  string::npos && queryInstream.good()) {
            queryInstream.get(sdfbuf, '$');
            queryInstream.get(sdfbuf, '\n');
        }
        queryCompound.read(sdfbuf.str());

                
        int cnt = 0;
        clock_t totalTimeUsed = 0;
        while (dbInstream.good()) {
            stringbuf sdfbuf;
            dbInstream.get(sdfbuf, '$');
            dbInstream.get(sdfbuf, '\n');
            dbInstream.ignore(1);
            while (sdfbuf.str().rfind("$$$$") ==  string::npos && dbInstream.good()) {
                dbInstream.get(sdfbuf, '$');
                dbInstream.get(sdfbuf, '\n');
                dbInstream.ignore(1);
            }
            if (sdfbuf.str().rfind("$$$$") !=  string::npos) {
                
                MCSCompound targetCompound;

                targetCompound.read(sdfbuf.str());

                MCS mcs(queryCompound, targetCompound,
                    userDefinedLowerBound, substructureNumLimit,
                    atomMismatchLowerBound, atomMismatchUpperBound,
                    boundMismatchLowerBound, boundMismatchUpperBound,
                    matchType, MCS::FAST, 0);

                clock_t start = clock();
                mcs.calculate();
                clock_t end = clock();
                
                totalTimeUsed += (end - start);
            
                resOstream << targetCompound.getCompoundName() << " " << targetCompound.size() << " " << queryCompound.size() << " " << mcs.size();
                if (mcs.isTimeout()) {
                    resOstream << " *";
                }
                resOstream << endl;
            }
        }
                
        double timeUsed = (double)totalTimeUsed / CLOCKS_PER_SEC;     
}*/

} // extern "C"
