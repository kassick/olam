// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/baseeventtype.cc"
// Created: "Ter, 04 Out 2011 12:03:25 -0300 (kassick)"
// Updated: "Qui, 06 Out 2011 15:19:33 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  baseeventtype.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 12:03:25 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "baseeventtype.hh"
#include "container.hh"
#include "paje_functions.hh"

#include <string>
#include <iostream>


void Paje::BaseEventType::do_header(ostream &out)
{
  //pajeDefineStateType(typeName, container->typeName, typeName, out);
  // Here it should call a pajeDefineEventType that is currently not
  // implemented
}


Paje::BaseEventType::BaseEventType(const string & typeName) {
  this->typeName = typeName;
  this->container = NULL;
}

Paje::BaseEventType::BaseEventType(const string & typeName,Paje::Container * c) {
  this->typeName = typeName;
  this->container = c;

}

#if 0
// this is unused for now
Paje::BaseEventType::BaseEventType(string & typeName, attribs_t * attribs)
{
  this->typeName = typeName;

  // walks the parent node (a Container, we hope...)
  walk_tree_up_to_root(attribs,[&](attrib_t * n, int level) {
      SemanticAttribute * attr = n->getVal();
      if (attr->id == ID_CONTAINER) {
        if (! container_type_names->count(attr->vals.name))
        {
          cerr << "Error: BaseEventType " << typeName << "can not find it's parent container " << attr->vals.name << endl;
          exit(1);
        }

        this->container = (*container_type_names)[attr->vals.name];
        return true; // stop now, otherwise it will get to the first containers!
      }
    });
}
#endif
