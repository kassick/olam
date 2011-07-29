// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.cc"
// Created: "Qua, 27 Jul 2011 11:07:19 -0300 (kassick)"
// Updated: "Qui, 28 Jul 2011 18:36:19 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  container.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27-07-2011 11:07:19 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "tree.hh"
#include <iostream>
#include <map>
#include "container.hh"

namespace Paje {

  string Container::toPaje() {
    return "PAjeCreateContainer " + this->typeName;

  } Container::Container() {
    triggerParent = destroyChildren = false;
  }


  Container::Container(attribs_t::iterator it) {
    while (attributes->is_valid(it)) {
      cout << "attrib has id " << (*it)->id << endl;
      it++;
    }


    // here is a case treating every parameter...

  }

}


//Globals
attribs_t *attributes;
hierarchy_t *toplevel_hierarchy;




//Utility functions

struct semantic_attribute *new_semantic_attribute()
{
  return new semantic_attribute;
}


void reparent_containers(Paje::Container * c, attribs_t::iterator it)
{
}

void init_desc_parser()
{
  toplevel_hierarchy = new tree < Paje::Container * >();
  attributes = new tree < struct semantic_attribute *>();
}

void print_tree(const attribs_t * tr, attribs_t::iterator it,
                attribs_t::iterator end)
{
  if (!tr->is_valid(it)) {
    cout << "Iterator is not valid!" << endl;
    return;
  }
  
  cout << "Begin has depth " << tr->depth(it) <<endl;
  cout << "end has depth " << tr->depth(end) <<endl;

  int rootdepth = tr->depth(it);
  std::cout << "-----" << std::endl;
  while (it != end) {
    cout << tr->depth(it) <<endl;
    for (int i = 0; i < tr->depth(it) - rootdepth; ++i)
      std::cout << "  ";
    std::cout << (*it) << std::endl << std::flush;
    ++it;
  }
  std::cout << "-----" << std::endl;
}
