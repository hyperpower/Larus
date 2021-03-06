/*
 * File: hashmap.h
 * ---------------
 * This file exports the <code>HashMap</code> class, which stores
 * a set of <i>key</i>-<i>value</i> pairs.
 * 
 * @version 2014/11/13
 * - added add() method as synonym for put()
 * - added template hashCode function
 * - moved hashCode functions to hashcode.h/cpp
 * @version 2014/10/29
 * - moved hashCode functions out to hashcode.h
 * @version 2014/10/10
 * - added comparison operators ==, !=
 */

#ifndef _hashmap_h
#define _hashmap_h

#include <cstdlib>
#include <map>
#include <string>
#include <sstream>

#include "hashcode.h"
#include "ArrayList.h"
//#include "vector.h"
namespace Larus {
/*
 * Class: HashMap<KeyType,ValueType>
 * ---------------------------------
 * This class implements an efficient association between
 * <b><i>keys</i></b> and <b><i>values</i></b>.  This class is
 * identical to the <a href="Map-class.html"><code>Map</code></a> class
 * except for the fact that it uses a hash table as its underlying
 * representation.  Although the <code>HashMap</code> class operates in
 * constant time, the iterator for <code>HashMap</code> returns the
 * values in a seemingly random order.
 */
template<typename KeyType, typename ValueType>
class HashMap {
public:
	/*
	 * Constructor: HashMap
	 * Usage: HashMap<KeyType,ValueType> map;
	 * --------------------------------------
	 * Initializes a new empty map that associates keys and values of
	 * the specified types.  The type used for the key must define
	 * the <code>==</code> operator, and there must be a free function
	 * with the following signature:
	 *
	 *<pre>
	 *    int hashCode(KeyType key);
	 *</pre>
	 *
	 * that returns a positive integer determined by the key.  This interface
	 * exports <code>hashCode</code> functions for <code>string</code> and
	 * the C++ primitive types.
	 */
	HashMap();

	/*
	 * Destructor: ~HashMap
	 * --------------------
	 * Frees any heap storage associated with this map.
	 */
	virtual ~HashMap();


	/*
	 * Method: clear
	 * Usage: map.clear();
	 * -------------------
	 * Removes all entries from this map.
	 */
	void Clear();

	/*
	 * Method: containsKey
	 * Usage: if (map.containsKey(key)) ...
	 * ------------------------------------
	 * Returns <code>true</code> if there is an entry for <code>key</code>
	 * in this map.
	 */
	bool HasKey(const KeyType& key) const;

	/*
	 * Method: equals
	 * Usage: if (map.equals(map2)) ...
	 * --------------------------------
	 * Returns <code>true</code> if the two maps contain exactly the same
	 * key/value pairs, and <code>false</code> otherwise.
	 */
	bool IsEqual(const HashMap& map2) const;

	/*
	 * Method: get
	 * Usage: ValueType value = map.get(key);
	 * --------------------------------------
	 * Returns the value associated with <code>key</code> in this map.
	 * If <code>key</code> is not found, <code>get</code> returns the
	 * default value for <code>ValueType</code>.
	 */
	ValueType Get(const KeyType& key) const;

	/*
	 * Method: isEmpty
	 * Usage: if (map.isEmpty()) ...
	 * -----------------------------
	 * Returns <code>true</code> if this map contains no entries.
	 */
	bool Empty() const;

	/*
	 * Method: keys
	 * Usage: arrayListT<KeyType> keys = map.keys();
	 * -----------------------------------------
	 * Returns a collection containing all keys in this map.
	 * Note that this implementation makes a deep copy of the keys,
	 * so it is inefficient to call on large maps.
	 */
	arrayListT<KeyType> keys() const;

	/*
	 * Method: mapAll
	 * Usage: map.mapAll(fn);
	 * ----------------------
	 * Iterates through the map entries and calls <code>fn(key, value)</code>
	 * for each one.  The keys are processed in an undetermined order.
	 */
	void mapAll(void (*fn)(KeyType, ValueType)) const;
	void mapAll(void (*fn)(const KeyType&, const ValueType&)) const;

	template<typename FunctorType>
	void mapAll(FunctorType fn) const;

	/*
	 * Method: put
	 * Usage: map.put(key, value);
	 * ---------------------------
	 * Associates <code>key</code> with <code>value</code> in this map.
	 * Any previous value associated with <code>key</code> is replaced
	 * by the new value.
	 */
	void Insert(const KeyType& key, const ValueType& value);

	/*
	 * Method: remove
	 * Usage: map.remove(key);
	 * -----------------------
	 * Removes any entry for <code>key</code> from this map.
	 * If the given key is not found, has no effect.
	 */
	void Remove(const KeyType& key);

	/*
	 * Method: putAll
	 * Usage: map.putAll(map2);
	 * ---------------------------
	 * Adds all key/value pairs from the given map to this map.
	 * If both maps contain a pair for the same key, the one from map2 will
	 * replace the one from this map.
	 * Returns a reference to this map.
	 */
	HashMap& InsertAll(const HashMap& map2);

	/*
	 * Method: removeAll
	 * Usage: map.removeAll(map2);
	 * ---------------------------
	 * Removes all key/value pairs from this map that are contained in the given map.
	 * If both maps contain the same key but it maps to different values, that
	 * mapping will not be removed.
	 * Returns a reference to this map.
	 */
	HashMap& RemoveAll(const HashMap& map2);

	/*
	 * Method: retainAll
	 * Usage: map.retainAll(map2);
	 * ---------------------------
	 * Removes all key/value pairs from this map that are not contained in the given map.
	 * If both maps contain the same key but it maps to different values, that
	 * mapping will be removed.
	 * Returns a reference to this map.
	 */
	HashMap& retainAll(const HashMap& map2);

	/*
	 * Method: size
	 * Usage: int nEntries = map.size();
	 * ---------------------------------
	 * Returns the number of entries in this map.
	 */
	int size() const;

	/*
	 * Method: toString
	 * Usage: string str = map.toString();
	 * -----------------------------------
	 * Converts the map to a printable string representation.
	 */
	std::string ToString() const;

	/*
	 * Method: values
	 * Usage: arrayListT<ValueType> values = map.values();
	 * -----------------------------------------------
	 * Returns a collection containing all values in this map.
	 * Note that this implementation makes a deep copy of the values,
	 * so it is inefficient to call on large maps.
	 */
	arrayListT<ValueType> values() const;

	/*
	 * Operator: []
	 * Usage: map[key]
	 * ---------------
	 * Selects the value associated with <code>key</code>.  This syntax
	 * makes it easy to think of a map as an "associative array"
	 * indexed by the key type.  If <code>key</code> is already present
	 * in the map, this function returns a reference to its associated
	 * value.  If key is not present in the map, a new entry is created
	 * whose value is set to the default for the value type.
	 */
	ValueType& operator [](const KeyType& key);
	ValueType operator [](const KeyType& key) const;

	/*
	 * Operator: ==
	 * Usage: if (map1 == map2) ...
	 * ----------------------------
	 * Compares two maps for equality.
	 */
	bool operator ==(const HashMap& map2) const;

	/*
	 * Operator: !=
	 * Usage: if (map1 != map2) ...
	 * ----------------------------
	 * Compares two maps for inequality.
	 */
	bool operator !=(const HashMap& map2) const;

	/*
	 * Operator: +
	 * Usage: map1 + map2
	 * ------------------
	 * Returns the union of the two maps, equivalent to a copy of the first map
	 * with addAll called on it passing the second map as a parameter.
	 * If the two maps both contain a mapping for the same key, the mapping
	 * from the second map is favored.
	 */
	HashMap operator +(const HashMap& map2) const;

	/*
	 * Operator: +=
	 * Usage: map1 += map2;
	 * --------------------
	 * Adds all key/value pairs from the given map to this map.
	 * Equivalent to calling addAll(map2).
	 */
	HashMap& operator +=(const HashMap& map2);

	/*
	 * Operator: -
	 * Usage: map1 - map2
	 * ------------------
	 * Returns the difference of the two maps, equivalent to a copy of the first map
	 * with removeAll called on it passing the second map as a parameter.
	 */
	HashMap operator -(const HashMap& map2) const;

	/*
	 * Operator: -=
	 * Usage: map1 -= map2;
	 * --------------------
	 * Removes all key/value pairs from the given map to this map.
	 * Equivalent to calling removeAll(map2).
	 */
	HashMap& operator -=(const HashMap& map2);

	/*
	 * Operator: *
	 * Usage: map1 * map2
	 * ------------------
	 * Returns the intersection of the two maps, equivalent to a copy of the first map
	 * with retainAll called on it passing the second map as a parameter.
	 */
	HashMap operator *(const HashMap& map2) const;

	/*
	 * Operator: *=
	 * Usage: map1 *= map2;
	 * ---------------------
	 * Removes all key/value pairs that are not found in the given map from this map.
	 * Equivalent to calling retainAll(map2).
	 */
	HashMap& operator *=(const HashMap& map2);

	/*
	 * Additional HashMap operations
	 * -----------------------------
	 * In addition to the methods listed in this interface, the HashMap
	 * class supports the following operations:
	 *
	 *   - Stream I/O using the << and >> operators
	 *   - Deep copying for the copy constructor and assignment operator
	 *   - Iteration using the range-based for statement and STL iterators
	 *
	 * The HashMap class makes no guarantees about the order of iteration.
	 */

	/* Private section */

	/**********************************************************************/
	/* Note: Everything below this point in the file is logically part    */
	/* of the implementation and should not be of interest to clients.    */
	/**********************************************************************/

	/*
	 * Implementation notes:
	 * ---------------------
	 * The HashMap class is represented using a hash table that uses
	 * bucket chaining to resolve collisions.
	 */
private:
	/* Constant definitions */
	static const int INITIAL_BUCKET_COUNT = 101;
	static const int MAX_LOAD_PERCENTAGE = 70;

	/* Type definition for cells in the bucket chain */
	struct Cell {
		KeyType key;
		ValueType value;
		Cell* next;
	};

	/* Instance variables */
	arrayListT<Cell*> buckets;
	int nBuckets;
	int numEntries;

	/* Private methods */

	/*
	 * Private method: createBuckets
	 * Usage: createBuckets(nBuckets);
	 * -------------------------------
	 * Sets up the vector of buckets to have nBuckets entries, each NULL.
	 * If asked to make empty vector, makes one bucket just to simplify
	 * handling elsewhere.
	 */
	void _CreateBuckets(int nBuckets) {
		if (nBuckets == 0) {
			nBuckets = 1;
		}
		buckets = arrayListT<Cell*>(nBuckets, NULL);
		this->nBuckets = nBuckets;
		numEntries = 0;
	}

	/*
	 * Private method: deleteBuckets
	 * Usage: deleteBuckets(buckets);
	 * ------------------------------
	 * Deletes all the cells in the linked lists contained in vector.
	 */
	void _DeleteBuckets(arrayListT<Cell*>& buckets) {
		for (int i = 0; i < buckets.size(); i++) {
			Cell* cp = buckets[i];
			while (cp != NULL) {
				Cell* np = cp->next;
				delete cp;
				cp = np;
			}
			buckets[i] = NULL;
		}
	}

	/*
	 * Private method: expandAndRehash
	 * Usage: expandAndRehash();
	 * -------------------------
	 * This method is used to increase the number of buckets in the map
	 * and then rehashes all existing entries and adds them into new buckets.
	 * This operation is used when the load factor (i.e. the number of cells
	 * per bucket) has increased enough to warrant this O(N) operation to
	 * enlarge and redistribute the entries.
	 */
	void _ExpandAndRehash() {
		arrayListT<Cell*> oldBuckets = buckets;
		_CreateBuckets(oldBuckets.size() * 2 + 1);
		for (int i = 0; i < oldBuckets.size(); i++) {
			for (Cell* cp = oldBuckets[i]; cp != NULL; cp = cp->next) {
				Insert(cp->key, cp->value);
			}
		}
		_DeleteBuckets(oldBuckets);
	}

	/*
	 * Private method: findCell
	 * Usage: Cell* cp = findCell(bucket, key);
	 *        Cell* cp = findCell(bucket, key, parent);
	 * ------------------------------------------------
	 * Finds a cell in the chain for the specified bucket that matches key.
	 * If a match is found, the return value is a pointer to the cell containing
	 * the matching key.  If no match is found, the function returns NULL.
	 * If the optional third argument is supplied, it is filled in with the
	 * cell preceding the matching cell to allow the client to splice out
	 * the target cell in the delete call.  If parent is NULL, it indicates
	 * that the cell is the first cell in the bucket chain.
	 */
	Cell* _FindCell(int bucket, const KeyType& key) const {
		Cell *dummy;
		return _FindCell(bucket, key, dummy);
	}

	Cell* _FindCell(int bucket, const KeyType& key, Cell*& parent) const {
		parent = NULL;
		Cell* cp = buckets.get(bucket);
		while (cp != NULL && key != cp->key) {
			parent = cp;
			cp = cp->next;
		}
		return cp;
	}

	void _DeepCopy(const HashMap& src) {
		_CreateBuckets(src.nBuckets);
		for (int i = 0; i < src.nBuckets; i++) {
			for (Cell* cp = src.buckets.get(i); cp != NULL; cp = cp->next) {
				Insert(cp->key, cp->value);
			}
		}
	}

public:
	/*
	 * Hidden features
	 * ---------------
	 * The remainder of this file consists of the code required to
	 * support deep copying and iteration.  Including these methods
	 * in the public interface would make that interface more
	 * difficult to understand for the average client.
	 */

	/*
	 * Deep copying support
	 * --------------------
	 * This copy constructor and operator= are defined to make a
	 * deep copy, making it possible to pass/return maps by value
	 * and assign from one map to another.
	 */
	HashMap& operator =(const HashMap& src) {
		if (this != &src) {
			Clear();
			_DeepCopy(src);
		}
		return *this;
	}

	HashMap(const HashMap& src) {
		_DeepCopy(src);
	}

	/*
	 * Iterator support
	 * ----------------
	 * The classes in the StanfordCPPLib collection implement input
	 * iterators so that they work symmetrically with respect to the
	 * corresponding STL classes.
	 */
	class iterator: public std::iterator<std::input_iterator_tag, KeyType> {
	private:
		const HashMap* mp; /* Pointer to the map           */
		int bucket; /* Index of current bucket      */
		Cell* cp; /* Current cell in bucket chain */

	public:
		iterator() :
				mp(NULL), bucket(0), cp(0) {
			/* Empty */
		}

		iterator(const HashMap* mp, bool end) {
			this->mp = mp;
			if (end) {
				bucket = mp->nBuckets;
				cp = NULL;
			} else {
				bucket = 0;
				cp = mp->buckets.get(bucket);
				while (cp == NULL && ++bucket < mp->nBuckets) {
					cp = mp->buckets.get(bucket);
				}
			}
		}

		iterator(const iterator& it) {
			mp = it.mp;
			bucket = it.bucket;
			cp = it.cp;
		}

		iterator& operator ++() {
			cp = cp->next;
			while (cp == NULL && ++bucket < mp->nBuckets) {
				cp = mp->buckets.get(bucket);
			}
			return *this;
		}

		iterator operator ++(int) {
			iterator copy(*this);
			operator++();
			return copy;
		}

		bool operator ==(const iterator& rhs) {
			return mp == rhs.mp && bucket == rhs.bucket && cp == rhs.cp;
		}

		bool operator !=(const iterator& rhs) {
			return !(*this == rhs);
		}

		KeyType& operator *() {
			return cp->key;
		}

		KeyType* operator ->() {
			return &cp->key;
		}

		friend class HashMap;
	};

	iterator begin() const {
		return iterator(this, false);
	}

	iterator end() const {
		return iterator(this, true);
	}
};

/*
 * Implementation notes: HashMap class
 * -----------------------------------
 * In this map implementation, the entries are stored in a hashtable.
 * The hashtable keeps a vector of "buckets", where each bucket is a
 * linked list of elements that share the same hash code (i.e. hash
 * collisions are resolved by chaining). The buckets are dynamically
 * allocated so that we can change the the number of buckets (rehash)
 * when the load factor becomes too high. The map should provide O(1)
 * performance on the put/remove/get operations.
 */
template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>::HashMap() {
	_CreateBuckets(INITIAL_BUCKET_COUNT);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>::~HashMap() {
	_DeleteBuckets(buckets);
}


template<typename KeyType, typename ValueType>
void HashMap<KeyType, ValueType>::Clear() {
	_DeleteBuckets(buckets);
	numEntries = 0;
}

template<typename KeyType, typename ValueType>
bool HashMap<KeyType, ValueType>::HasKey(const KeyType& key) const {
	return _FindCell(hashCode(key) % nBuckets, key) != NULL;
}

template<typename KeyType, typename ValueType>
bool HashMap<KeyType, ValueType>::IsEqual(
		const HashMap<KeyType, ValueType>& map2) const {
	// optimization: if literally same map, stop
	if (this == &map2) {
		return true;
	}

	if (size() != map2.size()) {
		return false;
	}

	// compare both ways; each must be subset of the other
	for (KeyType key : (*this)) {
		if (!map2.HasKey(key) || map2.Get(key) != Get(key)) {
			return false;
		}
	}
	for (KeyType key : map2) {
		if (!HasKey(key) || Get(key) != map2.Get(key)) {
			return false;
		}
	}
	return true;
}

template<typename KeyType, typename ValueType>
ValueType HashMap<KeyType, ValueType>::Get(const KeyType& key) const {
	Cell* cp = _FindCell(hashCode(key) % nBuckets, key);
	if (cp == NULL) {
		return ValueType();
	}
	return cp->value;
}

template<typename KeyType, typename ValueType>
bool HashMap<KeyType, ValueType>::Empty() const {
	return size() == 0;
}

template<typename KeyType, typename ValueType>
arrayListT<KeyType> HashMap<KeyType, ValueType>::keys() const {
	arrayListT<KeyType> keyset;
	for (KeyType key : *this) {
		keyset.push_back(key);
	}
	return keyset;
}

template<typename KeyType, typename ValueType>
void HashMap<KeyType, ValueType>::mapAll(void (*fn)(KeyType, ValueType)) const {
	for (int i = 0; i < buckets.size(); i++) {
		for (Cell* cp = buckets.get(i); cp != NULL; cp = cp->next) {
			fn(cp->key, cp->value);
		}
	}
}

template<typename KeyType, typename ValueType>
void HashMap<KeyType, ValueType>::mapAll(
		void (*fn)(const KeyType&, const ValueType&)) const {
	for (int i = 0; i < buckets.size(); i++) {
		for (Cell* cp = buckets.get(i); cp != NULL; cp = cp->next) {
			fn(cp->key, cp->value);
		}
	}
}

template<typename KeyType, typename ValueType>
template<typename FunctorType>
void HashMap<KeyType, ValueType>::mapAll(FunctorType fn) const {
	for (int i = 0; i < buckets.size(); i++) {
		for (Cell* cp = buckets.get(i); cp != NULL; cp = cp->next) {
			fn(cp->key, cp->value);
		}
	}
}

template<typename KeyType, typename ValueType>
void HashMap<KeyType, ValueType>::Insert(const KeyType& key,
		const ValueType& value) {
	(*this)[key] = value;
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::InsertAll(
		const HashMap& map2) {
	for (KeyType key : map2) {
		Insert(key, map2.Get(key));
	}
	return *this;
}

template<typename KeyType, typename ValueType>
void HashMap<KeyType, ValueType>::Remove(const KeyType& key) {
	int bucket = hashCode(key) % nBuckets;
	Cell *parent;
	Cell* cp = _FindCell(bucket, key, parent);
	if (cp != NULL) {
		if (parent == NULL) {
			buckets[bucket] = cp->next;
		} else {
			parent->next = cp->next;
		}
		delete cp;
		numEntries--;
	}
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::RemoveAll(
		const HashMap& map2) {
	for (KeyType key : map2) {
		if (HasKey(key) && Get(key) == map2.Get(key)) {
			Remove(key);
		}
	}
	return *this;
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::retainAll(
		const HashMap& map2) {
	arrayListT<KeyType> toRemove;
	for (typename arrayListT<Cell*>::const_iterator iter = buckets.begin();
			iter != buckets.end(); ++iter) {
		KeyType key = (*iter)->key;
		if (!map2.HasKey(key) || Get(key) != map2.Get(key)) {
			toRemove.add(key);
		}
	}
	for (KeyType key : toRemove) {
		Remove(key);
	}
	return *this;
}

template<typename KeyType, typename ValueType>
int HashMap<KeyType, ValueType>::size() const {
	return numEntries;
}

template<typename KeyType, typename ValueType>
std::string HashMap<KeyType, ValueType>::ToString() const {
	std::ostringstream os;
	os << *this;
	return os.str();
}

template<typename KeyType, typename ValueType>
arrayListT<ValueType> HashMap<KeyType, ValueType>::values() const {
	arrayListT<ValueType> values;
	for (typename arrayListT<Cell*>::const_iterator iter = buckets.begin();
			iter != buckets.end(); ++iter) {
		KeyType key = (*iter)->key;
		//for (KeyType key : *this) {
		values.add(this->Get(key));
	}
	return values;
}

template<typename KeyType, typename ValueType>
ValueType& HashMap<KeyType, ValueType>::operator [](const KeyType& key) {
	int bucket = hashCode(key) % nBuckets;
	Cell* cp = _FindCell(bucket, key);
	if (cp == NULL) {
		if (numEntries > MAX_LOAD_PERCENTAGE * nBuckets / 100.0) {
			_ExpandAndRehash();
			bucket = hashCode(key) % nBuckets;
		}
		cp = new Cell;
		cp->key = key;
		cp->value = ValueType();
		cp->next = buckets[bucket];
		buckets[bucket] = cp;
		numEntries++;
	}
	return cp->value;
}

template<typename KeyType, typename ValueType>
ValueType HashMap<KeyType, ValueType>::operator [](const KeyType& key) const {
	return Get(key);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType> HashMap<KeyType, ValueType>::operator +(
		const HashMap& map2) const {
	HashMap<KeyType, ValueType> result = *this;
	return result.InsertAll(map2);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::operator +=(
		const HashMap& map2) {
	return InsertAll(map2);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType> HashMap<KeyType, ValueType>::operator -(
		const HashMap& map2) const {
	HashMap<KeyType, ValueType> result = *this;
	return result.RemoveAll(map2);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::operator -=(
		const HashMap& map2) {
	return RemoveAll(map2);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType> HashMap<KeyType, ValueType>::operator *(
		const HashMap& map2) const {
	HashMap<KeyType, ValueType> result = *this;
	return result.retainAll(map2);
}

template<typename KeyType, typename ValueType>
HashMap<KeyType, ValueType>& HashMap<KeyType, ValueType>::operator *=(
		const HashMap& map2) {
	return retainAll(map2);
}

template<typename KeyType, typename ValueType>
bool HashMap<KeyType, ValueType>::operator ==(const HashMap& map2) const {
	return IsEqual(map2);
}

template<typename KeyType, typename ValueType>
bool HashMap<KeyType, ValueType>::operator !=(const HashMap& map2) const {
	return IsEqual(map2);
}

/*
 * Template hash function for hash maps.
 * Requires the key and value types in the HashMap to have a hashCode function.
 */
template<typename K, typename V>
int hashCode(const HashMap<K, V>& map) {
	int code = HASH_SEED;
	for (K k : map) {
		code = HASH_MULTIPLIER * code + hashCode(k);
		V v = map[k];
		code = HASH_MULTIPLIER * code + hashCode(v);
	}
	return int(code & HASH_MASK);
}

/*
 * Implementation notes: << and >>
 * -------------------------------
 * The insertion and extraction operators use the template facilities in
 * strlib.h to read and write generic values in a way that treats strings
 * specially.
 */
template<typename KeyType, typename ValueType>
std::ostream& operator <<(std::ostream& os,
		const HashMap<KeyType, ValueType>& map) {
	os << "{";
	typename HashMap<KeyType, ValueType>::iterator begin = map.begin();
	typename HashMap<KeyType, ValueType>::iterator end = map.end();
	typename HashMap<KeyType, ValueType>::iterator it = begin;
	while (it != end) {
		if (it != begin) {
			os << ", ";
		}
		os << (*it);
		os << ":";
		os << map[*it];
		++it;
	}
	return os << "}";
}

template<typename KeyType, typename ValueType>
std::istream& operator >>(std::istream& is, HashMap<KeyType, ValueType>& map) {
	char ch = '\0';
	is >> ch;
	if (ch != '{') {
		std::cerr << "HashMap::operator >>: Missing {\n";
	}
	map.clear();
	is >> ch;
	if (ch != '}') {
		is.unget();
		while (true) {
			KeyType key;
			is >> key;
			is >> ch;
			if (ch != ':') {
				std::cerr
						<< ("HashMap::operator >>: Missing colon after key\n");
			}
			ValueType value;
			is >> value;
			map[key] = value;
			is >> ch;
			if (ch == '}') {
				break;
			}
			if (ch != ',') {
				std::cerr
						<< (std::string(
								"HashMap::operator >>: Unexpected character ")
								+ ch + std::string("\n"));
			}
		}
	}
	return is;
}
}

#endif
