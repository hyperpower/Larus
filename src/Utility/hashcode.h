#ifndef _HASHCODE_H
#define _HASHCODE_H

#include <string>

namespace Larus{

const int HASH_SEED = 5381;               // Starting point for first cycle
const int HASH_MULTIPLIER = 33;           // Multiplier for each cycle
const int HASH_MASK = unsigned(-1) >> 1;  // All 1 bits except the sign

/*
 * Function: hashCode
 * Usage: int hash = hashCode(key);
 * --------------------------------
 * Returns a hash code for the specified key, which is always a
 * nonnegative integer.  This function is overloaded to support
 * all of the primitive types and the C++ <code>string</code> type.
 */
inline int hashCode(bool key);
inline int hashCode(char key);
inline int hashCode(double key);
inline int hashCode(float key);
inline int hashCode(int key);
inline int hashCode(long key);
inline int hashCode(const char* str);
inline int hashCode(const std::string& str);
inline int hashCode(void* key);
/*
 * Implementation notes: hashCode
 * ------------------------------
 * This function takes a string key and uses it to derive a hash code,
 * which is a nonnegative integer related to the key by a deterministic
 * function that distributes keys well across the space of integers.
 * The general method is called linear congruence, which is also used
 * in random-number generators.  The specific algorithm used here is
 * called djb2 after the initials of its inventor, Daniel J. Bernstein,
 * Professor of Mathematics at the University of Illinois at Chicago.
 */
inline int hashCode(bool key) {
    return (int) key;
}

inline int hashCode(char key) {
    return key;
}

inline int hashCode(double key) {
    char* byte = (char*) &key;
    unsigned hash = HASH_SEED;
    for (int i = 0; i < (int) sizeof(double); i++) {
        hash = HASH_MULTIPLIER * hash + (int) *byte++;
    }
    return hash & HASH_MASK;
}

inline int hashCode(float key) {
    char* byte = (char*) &key;
    unsigned hash = HASH_SEED;
    for (int i = 0; i < (int) sizeof(float); i++) {
        hash = HASH_MULTIPLIER * hash + (int) *byte++;
    }
    return hash & HASH_MASK;
}

inline int hashCode(int key) {
    return key & HASH_MASK;
}

inline int hashCode(long key) {
    return int(key) & HASH_MASK;
}

inline int hashCode(const char* str) {
    unsigned hash = HASH_SEED;
    for (int i = 0; str && str[i] != 0; i++) {
        hash = HASH_MULTIPLIER * hash + str[i];
    }
    return int(hash & HASH_MASK);
}

inline int hashCode(const std::string& str) {
    unsigned hash = HASH_SEED;
    int n = str.length();
    for (int i = 0; i < n; i++) {
        hash = HASH_MULTIPLIER * hash + str[i];
    }
    return int(hash & HASH_MASK);
}

inline int hashCode(void* key) {
    return hashCode(reinterpret_cast<long>(key));
}

}


#endif // _hashcode_h
