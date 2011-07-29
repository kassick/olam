// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/mytree/tree.cpp"
// Created: "Qui, 28 Jul 2011 20:31:25 -0300 (kassick)"
// Updated: "Qui, 28 Jul 2011 23:52:23 -0300 (kassick)"
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

#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include "tree.hpp"

using namespace std;


class thing_t {

  public:
    int _id;
    thing_t(int id) {
      _id = id;
    }
    
    ~thing_t() {
      cout << "Freeing " <<this->toString() << endl <<flush;
    }


    const string toString() const {
      stringstream t;
      t << "id is " << _id;
      return t.str();
    }
};

    template <typename CharT, typename Traits>
    basic_ostream<CharT, Traits>& operator<<(
    basic_ostream<CharT, Traits>& out, const thing_t& r)
    {
    return out<< r.toString();
    }


int main( ) {
  thing_t *a, *b, *c, *d, *e,*f,*g,*h;

  a = new thing_t(1);
  b = new thing_t(2);
  c = new thing_t(3);
  d = new thing_t(4);
  e = new thing_t(5);
  f = new thing_t(6);
  g = new thing_t(7);
  h = new thing_t(8);

   TreeNode<thing_t*> * node1 = new TreeNode<thing_t*>(a);
   TreeNode<thing_t*> * node2 = new TreeNode<thing_t*>(b);
   TreeNode<thing_t*> * node3 = new TreeNode<thing_t*>(c);
   TreeNode<thing_t*> * node4 = new TreeNode<thing_t*>(d);
   TreeNode<thing_t*> * node5 = new TreeNode<thing_t*>(e);
   TreeNode<thing_t*> * node6 = new TreeNode<thing_t*>(f);
   TreeNode<thing_t*> * node7 = new TreeNode<thing_t*>(g);
   TreeNode<thing_t*> * node8 = new TreeNode<thing_t*>(h);

   node1->addChild(node2);
   node1->addChild(node3);
   node2->addChild(node4);
   node2->addChild(node5);
   node4->addChild(node6);
   node4->addChild(node7);
   node4->addChild(node8);

   print_tree<thing_t*>(node1);

   //delete node2;
   //node1->removeChild(node2);
   free_subtree<thing_t*>(node2,true);
   cout <<"removed" <<endl <<flush;
   cout << "now printing from 4" <<endl;
   //print_tree<thing_t*>(node4);

   //cout << *d <<endl;
     cout <<"out" <<endl;


   print_tree<thing_t*>(node1);

   cout << "ok" <<endl <<flush;

   free_subtree(node1,true);

  cout << *a <<endl <<flush;
  cout << *e <<endl;


}
