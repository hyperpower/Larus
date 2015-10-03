/************************
 //  \file   ListT.h
 //  \brief
 // 
 //  \author czhou
 //  \date   21 juil. 2014 
 ***********************/
#ifndef _LIST_H_
#define _LIST_H_

#include <stddef.h>
#include <assert.h>
#include <iostream>
#include "Iterator.h"

namespace Larus {

template<typename TYPE>
class ListNode {
public:
	ListNode* prev;
	ListNode* next;
	TYPE node;

	ListNode() :
			node() {
		prev = NULL;
		next = NULL;
	}
	ListNode(const TYPE& _x) :
			node(_x) {
		prev = NULL;
		next = NULL;
	}
	ListNode(const ListNode<TYPE>& _x) {
		prev = _x.prev;
		next = _x.next;
		node = _x.node;
	}
	ListNode<TYPE>& operator=(const ListNode<TYPE>& _x) {
		if (this == &_x) {
			return *this;
		} else {
			prev = _x.prev;
			next = _x.next;
			node = _x.node;
			return *this;
		}
	}
};


template<class _Tp, class _Ref, class _Ptr>
class _List_iterator {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef bidirectional_iterator_tag iterator_category;

	typedef _List_iterator<_Tp, _Tp&, _Tp*> iterator;
	typedef _List_iterator<_Tp, const _Tp&, const _Tp*> const_iterator;
	typedef _List_iterator<_Tp, _Ref, _Ptr> _Self;

	typedef _Tp value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;
	typedef ListNode<_Tp> _Node;

	ListNode<_Tp>* _ptr;

	_List_iterator() {
		_ptr = NULL_PTR;
	}
	_List_iterator(ListNode<_Tp>* _x) {
		this->_ptr = _x;
	}
	_List_iterator(const iterator& _x) {
		this->_ptr = _x._ptr;
	}

	void _incr() {
		_ptr = _ptr->next;
	}
	void _decr() {
		_ptr = _ptr->prev;
	}

	bool operator==(const _Self& _x) const {
		return _ptr == _x._ptr;
	}
	bool operator!=(const _Self& _x) const {
		return _ptr != _x._ptr;
	}

	reference operator*() const {
		return ((_Node*) _ptr)->node;
	}

	pointer operator->() const {
		return &(operator*());
	}

	_Self& operator++() {
		this->_incr();
		return *this;
	}
	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}
	_Self& operator--() {
		this->_decr();
		return *this;
	}
	_Self operator--(int) {
		_Self __tmp = *this;
		this->_decr();
		return __tmp;
	}
};

template<typename TYPE>
class ListT {
protected:
	typedef ListNode<TYPE>* pNode;
	typedef ListNode<TYPE> Node;
	//typedef TYPE& Ref_ND;
	ListNode<TYPE>* _mNode;        //!< The _mNode item of the list.
	//ListNode<TYPE>* first;       //!< The first item of the list.
	LarusDef::size_type count;   //!< Zero-based count of items within the list.
public:
	// type definitions
	typedef TYPE value_type;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;

	typedef _List_iterator<TYPE, TYPE&, TYPE*> iterator;
	typedef _List_iterator<TYPE, const TYPE&, const TYPE*> const_iterator;

	ListT();
	~ListT();
	ListT(const ListT<TYPE>&);
	ListT<TYPE>& operator=(const ListT<TYPE>&);

	//Iterators============================
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	//Capacity=============================
	size_type size() const;
	bool empty() const;

	//Element access=======================
	reference front();
	const_reference front() const;
	reference back();
	const_reference back() const;
	reference at(size_type);
	const_reference at(size_type) const;
	value_type get(size_type index) const;
	//Modifiers============================
	void push_back(const TYPE& listItem);
	void push_front(const TYPE& listItem);
	void pop_back();
	void pop_front();

	int findIdx(TYPE value);
	pNode findPnt(TYPE value);

	iterator erase(iterator iter);
	iterator erase_p(iterator iter);
	iterator insert(iterator pln, const TYPE& value);

	void clear();
};

template<typename TYPE>
ListT<TYPE>::ListT() {
	count = 0;
	_mNode = new Node();
	_mNode->next = _mNode;
	_mNode->prev = _mNode;
}

template<typename TYPE>
ListT<TYPE>::ListT(const ListT<TYPE>& a) {
	count = 0;
	_mNode = new Node();
	_mNode->next = _mNode;
	_mNode->prev = _mNode;
	if (a.size() == 0) {
		return;
	} else {
		for (const_iterator iter = a.begin(); iter != a.end(); iter++) {
			this->push_back(iter._ptr->node);
		}
	}
}

template<typename TYPE>
ListT<TYPE>& ListT<TYPE>::operator=(const ListT<TYPE>& a) {
	if (this == &a) {
		return *this;
	} else {
		clear();
		if (a.count == 0) {
			return *this;
		} else {
			for (const_iterator iter = a.begin(); iter != a.end(); iter++) {
				this->push_back(iter._ptr->node);
			}
		}
		return *this;
	}
}

template<typename TYPE>
ListT<TYPE>::~ListT() {
	if (count > 0) {
		clear();    //Clear the List in order to free memory
	}
	delete _mNode;
}

template<typename TYPE>
_List_iterator<TYPE, TYPE&, TYPE*> ListT<TYPE>::insert(
		_List_iterator<TYPE, TYPE&, TYPE*> pln, const TYPE& value) {
	pNode newItem = new Node(value);
	newItem->next = pln._ptr;
	newItem->prev = pln._ptr->prev;
	pln._ptr->prev->next = newItem;
	pln._ptr->prev = newItem;
	count++;
	return newItem;
}

template<typename TYPE>
void ListT<TYPE>::push_back(const TYPE& li) {
	insert(end(), li);
}

template<typename TYPE>
void ListT<TYPE>::push_front(const TYPE& li) {
	insert(begin(), li);
}

template<typename TYPE>
bool ListT<TYPE>::empty() const {
	if (count <= 0) {
		return true;
	} else {
		return false;
	}
}

template<typename TYPE>
_List_iterator<TYPE, TYPE&, TYPE*> ListT<TYPE>::erase(
		_List_iterator<TYPE, TYPE&, TYPE*> iter) {
	pNode prev_node = iter._ptr->prev;
	pNode next_node = iter._ptr->next;
	delete iter._ptr;
	prev_node->next = next_node;
	next_node->prev = prev_node;
	count--;
	return next_node;
}

template<typename TYPE>
_List_iterator<TYPE, TYPE&, TYPE*> ListT<TYPE>::erase_p(
		_List_iterator<TYPE, TYPE&, TYPE*> iter) {
	pNode prev_node = iter._ptr->prev;
	pNode next_node = iter._ptr->next;
	delete iter._ptr;
	prev_node->next = next_node;
	next_node->prev = prev_node;
	count--;
	return prev_node;
}

template<typename TYPE>
void ListT<TYPE>::pop_back() {
	if (empty()) {
		return;
	} else {
		iterator _tmp = end();
		erase(--_tmp);
	}
}

template<typename TYPE>
void ListT<TYPE>::pop_front() {
	if (empty()) {
		return;
	} else {
		erase(begin());
	}
}

template<typename TYPE>
TYPE ListT<TYPE>::get(LarusDef::size_type index) const {
	ListNode<TYPE>* pLItem = _mNode->next;
	if (index < count && index >= 0) {
		int i = 0;
		while (i < index) {
			pLItem = pLItem->next;
			i++;
		}
	} else {
		std::cerr << "!> index=" << index << " out of range, max index= "
				<< count - 1 << std::endl;
		assert(false);
	}
	return pLItem->node;
}

template<typename TYPE>
TYPE& ListT<TYPE>::at(LarusDef::size_type index) {
	ListNode<TYPE>* pLItem = _mNode->next;
	if (index < count && index >= 0) {
		int i = 0;
		while (i < index) {
			pLItem = pLItem->next;
			i++;
		}
	} else {
		std::cerr << "!> index=" << index << " out of range, max index= "
				<< count - 1 << std::endl;
		assert(false);
	}
	return pLItem->node;
}

template<typename TYPE>
const TYPE& ListT<TYPE>::at(LarusDef::size_type index) const {
	ListNode<TYPE>* pLItem = _mNode->next;
	if (index < count && index >= 0) {
		int i = 0;
		while (i < index) {
			pLItem = pLItem->next;
			i++;
		}
	} else {
		std::cerr << "!> index=" << index << " out of range, max index= "
				<< count - 1 << std::endl;
		assert(false);
	}
	return pLItem->node;
}

template<typename TYPE>
TYPE& ListT<TYPE>::front() {
	return _mNode->next->node;
}

template<typename TYPE>
const TYPE& ListT<TYPE>::front() const {
	return _mNode->next->node;
}

template<typename TYPE>
TYPE& ListT<TYPE>::back() {
	return _mNode->node;
}

template<typename TYPE>
const TYPE& ListT<TYPE>::back() const {
	return _mNode->node;
}

template<typename TYPE>
typename ListT<TYPE>::iterator ListT<TYPE>::begin() {
	return _mNode->next;
}

template<typename TYPE>
typename ListT<TYPE>::iterator ListT<TYPE>::end() {
	return _mNode;
}

template<typename TYPE>
typename ListT<TYPE>::const_iterator ListT<TYPE>::begin() const {
	return _mNode->next;
}

template<typename TYPE>
typename ListT<TYPE>::const_iterator ListT<TYPE>::end() const {
	return _mNode;
}

template<typename TYPE>
int ListT<TYPE>::size() const {
	return count;
}

template<typename TYPE>
void ListT<TYPE>::clear() {
	pNode _tmpPtr = _mNode->next;
	while (count > 0) {
		ListNode<TYPE>* _tmpItem = _tmpPtr;
		_tmpPtr = _tmpPtr->next;
		delete _tmpItem;
		count--;
	}
	_mNode->next = _mNode;
	_mNode->prev = _mNode;
}

template<typename TYPE>
int ListT<TYPE>::findIdx(TYPE value) {
	if (count <= 0) {
		return -1;
	} else {
		ListNode<TYPE>* pLItem = _mNode->next;
		int i = 0;
		while (i < count) {
			if (pLItem->node == value) {
				return i;
			} else {
				pLItem = pLItem->next;
				i++;
			}
		}
		return -1;
	}
}

template<typename TYPE>
ListNode<TYPE>* ListT<TYPE>::findPnt(TYPE value) {
	if (count <= 0) {
		return -1;
	} else {
		ListNode<TYPE>* pLItem = _mNode->next;
		int i = 0;
		while (i < count) {
			if (pLItem->node == value) {
				return pLItem;
			} else {
				pLItem = pLItem->next;
				i++;
			}
		}
		return -1;
	}
}

}
#endif /* LISTT_H_ */
