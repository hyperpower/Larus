/************************
 //  \file   BSTree.h
 //  \brief
 // 
 //  \author czhou
 //  \date   19 févr. 2014 
 ***********************/
#ifndef BTREE_H_
#define BTREE_H_

#include "../TypeDef.h"
#include "BinaryTree.h"
#include "ArrayList.h"
#include <queue>
#include <iostream>

namespace Larus {

template<typename TYPE>
class BSTree: public BinaryTree<TYPE> {
public:
	typedef typename BinaryTree<TYPE>::reference reference;
	typedef typename BinaryTree<TYPE>::const_reference const_reference;
	typedef typename BinaryTree<TYPE>::pointer pointer;
	typedef typename BinaryTree<TYPE>::const_pointer const_pointer;
	typedef typename BinaryTree<TYPE>::iterator iterator;
	typedef typename BinaryTree<TYPE>::const_iterator const_iterator;

	typedef typename BinaryTree<TYPE>::difference_type difference_type;
	typedef typename BinaryTree<TYPE>::size_type size_type;
protected:
	typedef BinaryTreeNode<TYPE>* pTreeNode;

	void _insert(pTreeNode &Current, const TYPE &data);
	void _erase(pTreeNode &Current, const TYPE &data);
	void _delete(pTreeNode &Current);

	void _printtree(pTreeNode Current, int layer);
	pTreeNode _find(pTreeNode& Current, const TYPE& data);
public:

	BSTree() :
			BinaryTree<TYPE>() {
	}

	BSTree(TYPE* data, int n);
	BSTree(const arrayListT<TYPE>&);

	bool hasValue(const TYPE &data);
	void Insert(const TYPE &x);
	void Insert(const arrayListT<TYPE>&);
	void erase(const TYPE &x);
	void erase(iterator &x);

	void PrintTree();

	TYPE MinVal(pTreeNode Current) const;
	TYPE MaxVal(pTreeNode Current) const;
	pTreeNode MinValNode(pTreeNode Current) const;
	pTreeNode MaxValNode(pTreeNode Current) const;

	TYPE MinVal() const;
	TYPE MaxVal() const;
	pTreeNode MinValNode() const;
	pTreeNode MaxValNode() const;

	pTreeNode Find_p(const TYPE &data);
	iterator  Find(const TYPE &data);

	TYPE Next(const TYPE &data) const;
	TYPE Next(pTreeNode &Current) const;
	pTreeNode NextNode(BinaryTreeNode<TYPE> *&Current) const;

	TYPE Previous(const TYPE& data) const;
	TYPE Previous(pTreeNode &Current) const;
	pTreeNode PreviousNode(BinaryTreeNode<TYPE> *&Current) const;

	const BinaryTreeNode<TYPE>* getRoot() const;
};
template<typename TYPE>
BSTree<TYPE>::BSTree(TYPE* data, int n) :
		BinaryTree<TYPE>() {
	for (int i = 0; i < n; i++)
		Insert(data[i]);
}
template<typename TYPE>
BSTree<TYPE>::BSTree(const arrayListT<TYPE>& arr) :
		BinaryTree<TYPE>() {
	Insert(arr);
}

template<typename TYPE>
TYPE BSTree<TYPE>::MinVal(BinaryTreeNode<TYPE>* Current) const {
	return MinValNode(Current)->m_value;
}

template<typename TYPE>
TYPE BSTree<TYPE>::MaxVal(BinaryTreeNode<TYPE>* Current) const {
	return MaxValNode(Current)->m_value;
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::MinValNode(
		BinaryTreeNode<TYPE>* Current) const {
	BinaryTreeNode<TYPE> *pnode = Current;
	if (pnode->lchild) {
		while (pnode->lchild)
			pnode = pnode->lchild;
	}
	return pnode;
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::MaxValNode(
		BinaryTreeNode<TYPE>* Current) const {
	BinaryTreeNode<TYPE> *pnode = Current;
	if (pnode->rchild) {
		while (pnode->rchild)
			pnode = pnode->rchild;
	}
	return pnode;
}

template<typename TYPE>
TYPE BSTree<TYPE>::MinVal() const {
	if (this->_root->lchild == NULL_PTR) {
		return TYPE();
	}
	return MinVal(this->_root);
}

template<typename TYPE>
TYPE BSTree<TYPE>::MaxVal() const {
	if (this->_root->lchild == NULL_PTR) {
		return TYPE();
	}
	return MaxVal(this->_root);
}
template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::MinValNode() const {
	if (this->_root->lchild == NULL_PTR) {
		return NULL_PTR;
	}
	return MinValNode(this->_root);
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::MaxValNode() const {
	if (this->_root->lchild == NULL_PTR) {
		return NULL_PTR;
	}
	return MaxValNode(this->_root);
}

template<typename TYPE>
TYPE BSTree<TYPE>::Next(const TYPE &data) const {
	BinaryTreeNode<TYPE>* pnode = Find(data);
	return Next(pnode);
}

template<typename TYPE>
TYPE BSTree<TYPE>::Next(BinaryTreeNode<TYPE> *&Current) const {
	BinaryTreeNode<TYPE>* fn;
	if (Current != NULL) {
		if (Current->rchild)
			return MinVal(Current->rchild);
		fn = Current->father;
		while (fn && Current == fn->rchild) {
			Current = fn;
			fn = fn->father;
		}
		if (fn)
			return fn->m_value;
		else
			return TYPE();
	}
	return TYPE();
}

template<typename TYPE>
BinaryTreeNode<TYPE> * BSTree<TYPE>::NextNode(
		BinaryTreeNode<TYPE> *&Current) const {
	if (Current == this->_root) {
		return this->_root;
	}
	BinaryTreeNode<TYPE>* fn;
	if (Current != NULL) {
		if (Current->rchild)
			return MinValNode(Current->rchild);
		fn = Current->father;
		while (fn && Current == fn->rchild) {
			Current = fn;
			fn = fn->father;
		}
		if (fn)
			return fn;
		else
			return NULL;
	}
	return NULL;
}

template<class TYPE>
TYPE BSTree<TYPE>::Previous(const TYPE& data) const {
	BinaryTreeNode<TYPE>* Current = Find(data);
	return Previous(Current);
}

template<class TYPE>
TYPE BSTree<TYPE>::Previous(BinaryTreeNode<TYPE> *&Current) const {
	BinaryTreeNode<TYPE>* fn;
	if (Current != NULL_PTR) {
		if (Current->lchild)
			return MaxVal(Current->lchild);
		fn = Current->father;
		while (fn && Current == fn->lchild) {
			Current = fn;
			fn = Current->father;
		}
		if (fn)
			return fn->m_value;
		else
			return TYPE();
	}
	return TYPE();
}

template<class TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::PreviousNode(
		BinaryTreeNode<TYPE> *&Current) const {
	BinaryTreeNode<TYPE>* fn;
	if (Current != NULL) {
		if (Current->lchild)
			return MaxVal(Current->lchild);
		fn = Current->father;
		while (fn && Current == fn->lchild) {
			Current = fn;
			fn = Current->father;
		}
		if (fn)
			return fn;
		else
			return NULL;
	}
	return NULL;
}

template<typename TYPE>
void BSTree<TYPE>::_printtree(BinaryTreeNode<TYPE> *Current, int layer) {
	int i;
	if (Current == NULL) {
		return;
	}
	_printtree(Current->rchild, layer + 1);
	for (i = 0; i < layer; i++) {
		if (layer == 1) {
			std::cout << "  R>";
		} else {
			std::cout << "   ";
		}
	}
	std::cout << "O" << '\n';
	_printtree(Current->lchild, layer + 1);
}

template<typename TYPE>
void BSTree<TYPE>::PrintTree() {
	_printtree(this->_root->lchild, 1);
}

template<typename TYPE>
void BSTree<TYPE>::_insert(BinaryTreeNode<TYPE> *&Current, const TYPE &data) {
	if (Current == NULL_PTR) {
		Current = new BinaryTreeNode<TYPE>(data);
		assert(Current != NULL_PTR);
	} else if (data < Current->m_value) {
		_insert(Current->lchild, data);
		Current->lchild->father = Current;
	} else {   //data >= Current->m_value
		_insert(Current->rchild, data);
		Current->rchild->father = Current;
	}
}

template<class TYPE>
void BSTree<TYPE>::Insert(const TYPE &x) {
	if (this->_root->lchild == NULL_PTR) {
		_insert(this->_root->lchild, x);
		this->_root->lchild->father = this->_root;
	} else {
		_insert(this->_root->lchild, x);
	}
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::_find(BinaryTreeNode<TYPE>* &Current,
		const TYPE &data) {
	if (Current != NULL) {
		BinaryTreeNode<TYPE> *temp = Current;
		while (temp != NULL) {
			if (temp->m_value == data)
				return temp;
			if (temp->m_value < data)
				temp = temp->rchild;
			else
				temp = temp->lchild;
		}
	}
	return NULL;
}
template<typename TYPE>
void BSTree<TYPE>::_delete(pTreeNode &Current) {
	if (Current->lchild && Current->rchild) {    //the node has two children
		BinaryTreeNode<TYPE>* temp = Current->rchild;
		while (temp->lchild != NULL) {
			temp = temp->lchild;
		}
		//把右子树中最小节点的值赋值给本节点
		Current->m_value = temp->m_value;
		_delete(temp); //删除右子树中最小值的节点
		Current = NextNode(Current);
		//Current->rchild = temp;
	} else { //this node has 1 or 0 children
		BinaryTreeNode<TYPE>* temp = Current;
		if (Current->lchild == NULL_PTR) {
			if (Current->rchild != NULL_PTR) {  //has left son
				Current->rchild->father = Current->father;
				if (Current == Current->father->rchild) {
					Current->father->rchild = Current->rchild;
				} else {
					Current->father->lchild = Current->rchild;
				}
			} else {  //no child
				if (Current == Current->father->rchild) {
					Current->father->rchild = NULL_PTR;
				} else {
					Current->father->lchild = NULL_PTR;
				}
				Current->rchild = NULL_PTR;
			}
		} else if (Current->rchild == NULL_PTR) { //has left son
			if (Current->lchild != NULL_PTR) {
				Current->lchild->father = Current->father;
				if (Current == Current->father->rchild) {
					Current->father->rchild = Current->lchild;
				} else {
					Current->father->lchild = Current->lchild;
				}
			}
		}
		Current = NextNode(Current);
		delete (temp);
	}
}

template<typename TYPE>
void BSTree<TYPE>::_erase(BinaryTreeNode<TYPE>* &Current, const TYPE &data) {
	if (Current == NULL_PTR)
		return; //not find node data equal to data;
	if (data < Current->m_value) {
		_erase(Current->lchild, data); //data less than Current, delete in left tree
	} else if (data > Current->m_value) {
		_erase(Current->rchild, data); //data less than Current, delete in right tree
	} else {                        //if data equal, this is the node for delete
		_delete(Current);
	}
	return;
}
//删除接口
template<typename TYPE>
void BSTree<TYPE>::erase(const TYPE &x) {
	_erase(this->_root->lchild, x);
}
template<typename TYPE>
void BSTree<TYPE>::erase(BSTree<TYPE>::iterator &iter) {
	if (iter.isExist()) {
		_delete(iter._ptr);
	}
}
template<typename TYPE>
void BSTree<TYPE>::Insert(const arrayListT<TYPE>& data) {
	for (typename arrayListT<TYPE>::size_type i = 0; i < data.size(); i++)
		Insert(data[i]);
}

template<typename TYPE>
bool BSTree<TYPE>::hasValue(const TYPE &data) {
	BinaryTreeNode<TYPE>* temp = _find(data, this->_root->lchild);
	if (temp == NULL)
		return false;
	else
		return true;
}

template<typename TYPE>
BinaryTreeNode<TYPE>* BSTree<TYPE>::Find_p(const TYPE &data) {
	return _find(this->_root->lchild, data);
}
template<typename TYPE>
typename BSTree<TYPE>::iterator BSTree<TYPE>::Find(const TYPE &data) {
	return BSTree<TYPE>::iterator(_find(this->_root->lchild, data));
}

template<typename TYPE>
const BinaryTreeNode<TYPE>* BSTree<TYPE>::getRoot() const {
	return this->_root->lchild;
}

} //namespace Larus
#endif /* BTREE_H_ */
