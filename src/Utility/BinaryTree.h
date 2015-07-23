/*
 * BinaryTree.h
 *
 *  Created on: Dec 22, 2014
 *      Author: zhou
 */

#ifndef _BINARYTREE_H_
#define _BINARYTREE_H_

#include "../TypeDef.h"
#include "Iterator.h"

namespace Larus {

template<typename TYPE>
class BinaryTreeNode: public ObjectBase {
protected:
	typedef BinaryTreeNode<TYPE>* pTreeNode;
public:
	TYPE m_value;
	pTreeNode lchild;
	pTreeNode rchild;
	pTreeNode father;
	BinaryTreeNode();
	BinaryTreeNode(const TYPE &data);
	BinaryTreeNode(const TYPE &data, pTreeNode lc, pTreeNode rc, pTreeNode fc);
	pTreeNode leftmost();
	pTreeNode rightmost();
};

template<typename TYPE>
BinaryTreeNode<TYPE>::BinaryTreeNode() {
	m_value = TYPE();
	lchild = NULL_PTR;
	rchild = NULL_PTR;
	father = NULL_PTR;
}

template<typename TYPE>
BinaryTreeNode<TYPE>::BinaryTreeNode(const TYPE &data) {
	m_value = data;
	lchild = NULL_PTR;
	rchild = NULL_PTR;
	father = NULL_PTR;
}

template<typename TYPE>
BinaryTreeNode<TYPE>::BinaryTreeNode(const TYPE &data, pTreeNode lc,
		pTreeNode rc, pTreeNode fc) {
	m_value = data;
	lchild = lc;
	rchild = rc;
	father = fc;
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BinaryTreeNode<TYPE>::leftmost() {
	pTreeNode ptn = this;
	while (ptn->lchild != NULL_PTR)
		ptn = ptn->lchild;
	return ptn;
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BinaryTreeNode<TYPE>::rightmost() {
	pTreeNode ptn = this;
	while (ptn->rchild != NULL_PTR)
		ptn = ptn->rchild;
	return ptn;
}

//This is the end of Class BinaryTreeNode
//-----------------------------------------------

template<class _Tp, class _Ref, class _Ptr>
class _BinaryTree_iterator {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef bidirectional_iterator_tag iterator_category;

	typedef _BinaryTree_iterator<_Tp, _Tp&, _Tp*> iterator;
	typedef _BinaryTree_iterator<_Tp, const _Tp&, const _Tp*> const_iterator;
	typedef _BinaryTree_iterator<_Tp, _Ref, _Ptr> _Self;

	typedef _Tp  value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;
	typedef BinaryTreeNode<_Tp> _Node;

	BinaryTreeNode<_Tp>* _ptr;

	_BinaryTree_iterator() {
		_ptr = NULL_PTR;
	}
	_BinaryTree_iterator(BinaryTreeNode<_Tp>* _x) {
		this->_ptr = _x;
	}
	_BinaryTree_iterator(const iterator& _x) {
		this->_ptr = _x._ptr;
	}

	void _incr() {
		BinaryTreeNode<_Tp>* fn;
		if (_ptr != NULL_PTR) {
			if (_ptr->rchild != NULL_PTR) {
				_ptr = _ptr->rchild->leftmost();
				return;
			}
			fn = _ptr->father;
			while (fn && _ptr == fn->rchild) {
				_ptr = fn;
				fn = fn->father;
			}
			_ptr = fn;
		}
	}

	void _decr() {
		BinaryTreeNode<_Tp>* fn;
		if (_ptr != NULL) {
			if (_ptr->lchild) {
				_ptr = _ptr->lchild->rightmost();
				return;
			}
			fn = _ptr->father;
			while (fn && _ptr == fn->lchild) {
				_ptr = fn;
				fn = _ptr->father;
			}
		}
	}

	bool operator==(const _BinaryTree_iterator& _x) const {
		return _ptr == _x._ptr;
	}
	bool operator!=(const _BinaryTree_iterator& _x) const {
		return _ptr != _x._ptr;
	}

	reference operator*() const {
		return ((_Node*) _ptr)->m_value;
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

	bool isExist(){
		return _ptr!=NULL_PTR;
	}
};

//===============================================
template<typename TYPE>
class BinaryTree: public ObjectBase {
public:
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef _BinaryTree_iterator<TYPE, TYPE&, TYPE*> iterator;
	typedef _BinaryTree_iterator<TYPE, const TYPE&, const TYPE*> const_iterator;

	typedef LarusDef::size_type difference_type;
	typedef LarusDef::size_type size_type;
protected:
	typedef BinaryTreeNode<TYPE>* pTreeNode;
	typedef void (*pFun_BinaryTree)(pTreeNode, utPointer);

	pTreeNode _root;

	void _preorder(pTreeNode, pFun_BinaryTree, utPointer);
	void _inorder(pTreeNode, pFun_BinaryTree, utPointer);
	void _postorder(pTreeNode, pFun_BinaryTree, utPointer);

	void _destory(pTreeNode&);
	void _copy(pTreeNode&, const pTreeNode&);
	void _reverse(pTreeNode);
	size_type _height(pTreeNode) const;
	size_type _size(const pTreeNode&) const;

	iterator root();
	const_iterator root() const;
public:
	//constructor================================
	BinaryTree();
	BinaryTree(const BinaryTree<TYPE>&);
	//destructor ================================
	~BinaryTree();
	//operator= =================================
	BinaryTree<TYPE>& operator=(const BinaryTree<TYPE>&);
	//Traversal =================================
	void PreOrder(pFun_BinaryTree, utPointer);
	void InOrder(pFun_BinaryTree, utPointer);
	void PostOrder(pFun_BinaryTree, utPointer);
	//iterator===================================
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	//===========================================

	bool empty() const;
	size_type size() const;
	void reverse();
	void clear();
	size_type height() const;
};

template<typename TYPE>
BinaryTree<TYPE>::BinaryTree() {
	_root = new BinaryTreeNode<TYPE>();
}

template<typename TYPE>
BinaryTree<TYPE>::BinaryTree(const BinaryTree<TYPE>& a) {
	this->_root = NULL_PTR;
	_copy(this->_root, a._root);
}

template<typename TYPE>
BinaryTree<TYPE>::~BinaryTree() {
	_destory(_root);
}

template<typename TYPE>
BinaryTree<TYPE>& BinaryTree<TYPE>::operator=(
		const BinaryTree<TYPE>& original) {
	_destory(this->_root);
	this->_root=NULL_PTR;
	_copy(this->_root, original._root);
	return *this;
}

template<typename TYPE>
void BinaryTree<TYPE>::_destory(pTreeNode& Current) {
	if (Current != NULL_PTR) {
		_destory(Current->lchild);
		_destory(Current->rchild);
		delete Current;
		Current = NULL_PTR;
	}
}

template<class TYPE>
void BinaryTree<TYPE>::_copy(pTreeNode& Current, const pTreeNode& original) {
	if (Current == NULL_PTR) {
		Current = new BinaryTreeNode<TYPE>(original->m_value);
	}
	if (original->lchild != NULL_PTR) {
		Current->lchild = new BinaryTreeNode<TYPE>(original->lchild->m_value);
		Current->lchild->father = Current;
		_copy(Current->lchild, original->lchild);
	}
	if (original->rchild != NULL_PTR) {
		Current->rchild = new BinaryTreeNode<TYPE>(original->rchild->m_value);
		Current->rchild->father = Current;
		_copy(Current->rchild, original->rchild);
	}
}

template<typename TYPE>
void BinaryTree<TYPE>::_preorder(pTreeNode Current, pFun_BinaryTree visit,
		utPointer utp) {
	if (Current != NULL_PTR) {
		(*visit)(Current->m_value, utp);
		_preorder(Current->lchild, visit, utp);
		_preorder(Current->rchild, visit, utp);
	}
}

template<class TYPE>
void BinaryTree<TYPE>::PreOrder(pFun_BinaryTree visit, utPointer utp) {
	_preorder(_root->lchild, visit, utp);
}

template<typename TYPE>
void BinaryTree<TYPE>::_postorder(pTreeNode Current, pFun_BinaryTree visit,
		utPointer utp) {
	if (Current != NULL_PTR) {
		_postorder(Current->lchild, visit, utp);
		_postorder(Current->rchild, visit, utp);
		(*visit)(Current->m_value, utp);
	}
}

template<class TYPE>
void BinaryTree<TYPE>::PostOrder(pFun_BinaryTree visit, utPointer utp) {
	_postorder(_root->lchild, visit, utp);
}

template<typename TYPE>
void BinaryTree<TYPE>::_inorder(pTreeNode Current, pFun_BinaryTree visit,
		utPointer utp) {
	if (Current != NULL_PTR) {
		_inorder(Current->lchild, visit, utp);
		(*visit)(Current->m_value, utp);
		_inorder(Current->rchild, visit, utp);
	}
}

template<class TYPE>
void BinaryTree<TYPE>::InOrder(pFun_BinaryTree visit, utPointer utp) {
	_inorder(_root->lchild, visit, utp);
}

template<class TYPE>
_BinaryTree_iterator<TYPE, const TYPE&, const TYPE*> BinaryTree<TYPE>::begin() const {
	BinaryTreeNode<TYPE> *pnode = _root->lchild;
	if (pnode == NULL_PTR) {
		return _root;
	}
	if (pnode->lchild) {
		while (pnode->lchild)
			pnode = pnode->lchild;
	}
	return pnode;
}

template<class TYPE>
_BinaryTree_iterator<TYPE, TYPE&, TYPE*> BinaryTree<TYPE>::begin() {
	BinaryTreeNode<TYPE> *pnode = _root->lchild;
	if (pnode == NULL_PTR) {
		return _root;
	}
	if (pnode->lchild) {
		while (pnode->lchild)
			pnode = pnode->lchild;
	}
	return pnode;
}
template<class TYPE>
_BinaryTree_iterator<TYPE, TYPE&, TYPE*> BinaryTree<TYPE>::root() {
	return _root->lchild;
}
template<class TYPE>
_BinaryTree_iterator<TYPE, const TYPE&, const TYPE*> BinaryTree<TYPE>::root() const {
	return _root->lchild;
}

template<class TYPE>
_BinaryTree_iterator<TYPE, const TYPE&, const TYPE*> BinaryTree<TYPE>::end() const {
	return _root;
}

template<class TYPE>
_BinaryTree_iterator<TYPE, TYPE&, TYPE*> BinaryTree<TYPE>::end() {
	return _root;
}



template<typename TYPE>
bool BinaryTree<TYPE>::empty() const {
	return _root->lchild == NULL_PTR;
}

template<typename TYPE>
LarusDef::size_type BinaryTree<TYPE>::_size(const pTreeNode& Current) const {
	if (Current == NULL_PTR) {
		return 0;
	} else {
		return _size(Current->lchild) + _size(Current->rchild) + 1;
	}
}

template<typename TYPE>
LarusDef::size_type BinaryTree<TYPE>::size() const {
	return _size(_root->lchild);
}

template<class TYPE>
void BinaryTree<TYPE>::reverse() {
	_reverse(_root->lchild);
}

template<class TYPE>
void BinaryTree<TYPE>::_reverse(pTreeNode Current) {
	if (Current != NULL_PTR) {
		pTreeNode temp = Current->lchild;
		Current->lchild = Current->rchild;
		Current->rchild = temp;
		_reverse(Current->lchild);
		_reverse(Current->rchild);
	}
}

template<class TYPE>
void BinaryTree<TYPE>::clear() {
	_destory(_root->lchild);
	//_root = NULL_PTR;
}

template<class TYPE>
int BinaryTree<TYPE>::_height(pTreeNode Current) const {
	if (Current == NULL_PTR)
		return 0;
	else
		return 1 + max(_height(Current->lchild), _height(Current->rchild));
}

template<class TYPE>
LarusDef::size_type BinaryTree<TYPE>::height() const {
	return _height(_root->lchild);
}

}

#endif /* UTILITY_BINARYTREE_H_ */
