// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.hh"
// Created: "Qua, 27 Jul 2011 11:08:49 -0300 (kassick)"
// Updated: "Qui, 28 Jul 2011 18:11:24 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  container.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  27-07-2011 11:08:49 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */



#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <iostream>
#include <map>
#include <string>
#include "tree.hh"

using namespace std;


namespace Paje {
  class Container;
}


typedef tree<Paje::Container*> hierarchy_t;


enum _attrib_ids {
  ID_NAME,
  ID_FORMAT_NAME,
  ID_CREATE_EVENT,
  ID_CREATE_PARENT,
  ID_DESTROY_EVENT,
  ID_DESTROY_CHILDREN,
  ID_CONTAINER,

  ID_NOP
};

struct semantic_attribute {
  int id;
  union _atrib_vals {
    char *name;
    char *format_name;
    char *create_event;
    char *destroy_event;
    bool create_parent;
    bool destroy_children;
    Paje::Container *container;
  } vals;
};


typedef tree<struct semantic_attribute*> attribs_t;


namespace Paje {

  using namespace std;

  class Container {

    public:

      Container();
      Container(attribs_t::iterator);
      string typeName;
      string formatName;
      bool triggerParent,destroyChildren;
      string triggerEvent;
      void * eventTypes;

      string toPaje();

  }   ;



}



extern attribs_t * attributes;
extern hierarchy_t * toplevel_hierarchy;


struct semantic_attribute * new_semantic_attribute();
void reparent_containers(Paje::Container * c, attribs_t::iterator it);
void init_desc_parser();
void print_tree(const attribs_t * tr, attribs_t::iterator it,
                attribs_t::iterator end);


#endif
