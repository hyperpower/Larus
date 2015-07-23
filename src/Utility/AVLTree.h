/************************
 //  \file   AVLTree.h
 //  \brief
 // 
 //  \author zhou
 //  \date   23 févr. 2014 
 ***********************/
#ifndef AVLTREE_H_
#define AVLTREE_H_

#include "../TypeDef.h"
#include "BSTree.h"
#include "BinaryTree.h"

namespace Larus {

template<typename ND>
class AVLTree: public BSTree<ND> {
protected:
    typedef BinaryTreeNode<ND>* pTreeNode;
	void _insert(pTreeNode &Current, const ND &data);
	void _delete(pTreeNode &Current, const ND &data);
	void _singrotateL(pTreeNode  &k2); //左左情况下的旋转
	void _singrotateR(pTreeNode  &k2); //右右情况下的旋转
	void _doubrotateLR(pTreeNode &k3); //左右情况下的旋转
	void _doubrotateRL(pTreeNode &k3); //右左情况下的旋转
public:
	AVLTree() :
			BSTree<ND>() {
	}
	void Insert(const ND &x); //插入接口
	void Delete(const ND &x); //删除接口
};

template<typename ND>
void AVLTree<ND>::_insert(pTreeNode &Current, const ND &data) {
	if (Current == NULL) { //如果节点为空,就在此节点处加入x信息
		Current = new BinaryTreeNode<ND>(data);
	} else if (data < Current->m_value) { //如果x小于节点的值,就继续在节点的左子树中插入x
		_insert(Current->lchild, data);
		Current->lchild->father = Current;
		if (2 == (_height(Current->lchild) - _height(Current->rchild))) {
			if (data < Current->lchild->m_value) {
				_singrotateL(Current);
			} else {
				_doubrotateLR(Current);
			}
		}

	} else { //如果x大于节点的值,就继续在节点的右子树中插入x
		_insert(Current->rchild, data);
		Current->rchild->father = Current;
		if (2 == (_height(Current->rchild) - _height(Current->lchild))) { //如果高度之差为2的话就失去了平衡,需要旋转
			if (data > Current->rchild->m_value) {
				_singrotateR(Current);
			} else {
				_doubrotateRL(Current);
			}
		}

	}
}

template<typename ND>
void AVLTree<ND>::Insert(const ND &x) {
	_insert(BSTree<ND>::root->lchild, x);
}

template<typename ND>
void AVLTree<ND>::_delete(pTreeNode &Current, const ND &data) {
	if (Current == NULL) {
		return;	    //没有找到值是x的节点
	}
	if (data < Current->m_value) {
		_delete(Current->lchild, data);	    //如果x小于节点的值,就继续在节点的左子树中删除x
		if (2 == (_height(Current->rchild) - _height(Current->lchild))) {
			if (Current->rchild->lchild != NULL
					&& (_height(Current->rchild->lchild)
							> _height(Current->rchild->rchild))) {
				_doubrotateRL(Current);
			} else {
				_singrotateR(Current);
			}
		}
	}

	else if (data > Current->m_value) {
		_delete(Current->rchild, data);	    //如果x大于节点的值,就继续在节点的右子树中删除x
		if (2 == (_height(Current->lchild) - _height(Current->rchild))) {
			if (Current->lchild->rchild != NULL
					&& (_height(Current->lchild->rchild)
							> _height(Current->lchild->lchild))) {
				_doubrotateLR(Current);
			} else {
				_singrotateL(Current);
			}
		}
	}

	else	    //如果相等,此节点就是要删除的节点
	{
		if (Current->lchild && Current->rchild) {    //此节点有两个儿子
			BinaryTreeNode<ND>* temp = Current->rchild;	    //temp指向节点的右儿子
			while (temp->lchild != NULL) {
				temp = temp->lchild;  //找到右子树中值最小的节点
			}
			//把右子树中最小节点的值赋值给本节点
			Current->m_value = temp->m_value;
			//Current->freq = temp->freq;
			_delete(Current->rchild, temp->m_value);	          //删除右子树中最小值的节点
			if (2 == _height(Current->lchild) - _height(Current->rchild)) {
				if (Current->lchild->rchild != NULL
						&& (_height(Current->lchild->rchild)
								> _height(Current->lchild->lchild))) {
					_doubrotateLR(Current);
				} else {
					_singrotateL(Current);
				}
			}
		} else	            //此节点有1个或0个儿子
		{
			BinaryTreeNode<ND>* temp = Current;
			if (Current->lchild == NULL) {	            //有右儿子或者没有儿子
				if (Current->rchild != NULL)
					Current->rchild->father = Current->father;
				Current = Current->rchild;
			} else if (Current->rchild == NULL) {	            //有左儿子
				if (Current->lchild != NULL)
					Current->lchild->father = Current->father;
				Current = Current->lchild;
			}
			delete (temp);
			temp = NULL;
		}
	}
	if (Current == NULL)
		return;
	//Current->hgt = Max(height(Current->lchild), height(Current->rchild)) + 1;
	//return;
}

template<typename ND>
void AVLTree<ND>::Delete(const ND &x) {
	_delete(BSTree<ND>::root->lchild, x);
}

template<typename ND>
void AVLTree<ND>::_singrotateL(pTreeNode &k2) {
	BinaryTreeNode<ND>* k1;
	k1 = k2->lchild;
	k2->lchild = k1->rchild;
	k1->rchild = k2;
	k1->father = k2->father;
	k2->father = k1;
	if (k1->rchild)
		k1->rchild->father = k2;
	k2 = k1;
	//---

	//k2->hgt=Max(height(k2->lson),height(k2->rson))+1;
	//k1->hgt=Max(height(k1->lson),k2->hgt)+1;
}
template<typename ND>
void AVLTree<ND>::_singrotateR(pTreeNode &k2) {
	BinaryTreeNode<ND>* k1;
	k1 = k2->rchild;
	k2->rchild = k1->lchild;
	k1->lchild = k2;
	k1->father = k2->father;
	k2->father = k1;
	if (k1->lchild)
		k1->lchild->father = k2;
	k2 = k1;
}
template<typename ND>
void AVLTree<ND>::_doubrotateLR(pTreeNode &k3) {
	_singrotateR(k3->lchild);
	_singrotateL(k3);
}
template<typename ND>
void AVLTree<ND>::_doubrotateRL(pTreeNode &k3) {
	_singrotateL(k3->rchild);
	_singrotateR(k3);
}

}

#endif /* AVLTREE_H_ */
