// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/mytree/tree.hpp"
// Created: "Qui, 28 Jul 2011 20:31:25 -0300 (kassick)"
// Updated: "Qui, 28 Jul 2011 23:53:24 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  tree.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  28-07-2011 20:31:25 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#ifndef __TREE_HPP_H__
#define __TREE_HPP_H__

#include <iostream>
#include <list>
#include <string>

using namespace std;

template<typename T>
class TreeNode {
  public:
    typedef T value_type;

    typedef TreeNode<T> self_type;
  private:
     
    typedef list< TreeNode<T> * > children_t;
    
    T val_;
    children_t children;

  protected:
     TreeNode<T> * parent_;

  public:
    typedef typename children_t::iterator iterator;
    typedef typename children_t::const_iterator const_iterator;

    bool operator==( const self_type& that ) const {
      return this->getVal() == that.getVal();
    }


    iterator begin()
    {
      return iterator( children.begin() );
    }
   
    iterator end()
    {
      return iterator (children.end() );
    }

    const_iterator begin() const
    {
      return const_iterator( children.begin() );
    }

    const_iterator end() const
    {
      return const_iterator( children.end() );
    }
   

     TreeNode(const T& val)
     {
       val_ = val;
       parent_ = NULL;
     }

    ~TreeNode()
    {
      if (parent_)
        parent_->removeChild(this);

      //cout <<"Free!\n" <<flush;
    }

     const T& getVal() const {
       return(val_);
     }

     void setVal(const T& val) {
       val_ = val;
     }

     void removeChild(TreeNode<T>* p ) {
      removeChild(p->getVal());
       /*
       iterator it = find(children.begin(),children.end(),p);
       if (it != children.end()) {
         children.remove(it);
         p->parent = NULL;
       }*/

     }

     void removeChild(const T& val) {
       typename children_t::iterator it;

       for (it = children.begin(); it != children.end(); ++it)
       {
         if ((*it)->getVal() == val)
           break;
       }

       if (it != children.end())
       {
         children.erase(it);
         (*it)->parent_ = NULL;
       }
     }


     void addChild(TreeNode<T>* p) {
       p->parent_ = this;
       children.push_back(p);
     }

};


template <typename T>
void print_tree_( T val, typename TreeNode<T>::iterator it1, typename TreeNode<T>::iterator it2, int level) {
  int i = level;

  while (i)
  {
    cout << " ";
    --i;
  }
  cout << *val <<endl;
  while (it1 != it2) {
    print_tree_((*it1)->getVal(), (*it1)->begin(), (*it1)->end(),level+1);
    ++it1;
  }


}



template <typename T>
void print_tree(TreeNode<T> * t)
{
  print_tree_<T>(t->getVal(), t->begin(),t->end(),0);
}

template <typename T>
void free_subtree_(TreeNode<T> * t,bool free_vals)
{
  typename TreeNode<T>::const_iterator begin, end;

  //cout << "called for " << *(t->getVal()) <<endl;
  for (begin = t->begin(); begin != t->end(); )
  {
    TreeNode<T> *tmp = *begin;
    ++begin;
    free_subtree_(tmp,free_vals);
    delete tmp;
  }

  if (free_vals)
  {
    delete t->getVal();
    t->setVal(NULL);
  }

}

template <typename T>
void free_subtree(TreeNode<T> * t,bool free_vals)
{
  free_subtree_(t,free_vals);

  delete t;
}

#endif

