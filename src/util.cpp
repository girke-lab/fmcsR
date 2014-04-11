
#include "util.h"

#include <cctype>
#include <string>

using namespace std;

string getUpper(const string& s) {
    string upper(s);
    for (int i = 0; i < upper.length(); ++i) {
        upper[i] = toupper(upper[i]);
    }
    return upper;
}

