// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.hh"
// Created: "Qua, 27 Jul 2011 11:08:49 -0300 (kassick)"
// Updated: "Sex, 29 Jul 2011 18:52:10 -0300 (kassick)"
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
#include "tree.hpp"

using namespace std;


namespace Paje {
  class Container;
}


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

class SemanticAttribute {
  public:
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


    const string toString() const;
};

//typedef SemanticAttribute semantic_attribute;



typedef TreeNode<Paje::Container*> hierarchy_t;
typedef TreeNode<SemanticAttribute*> attribs_t;



namespace Paje {

  using namespace std;

  class Container {
    private:
      void fill_fields_from_cb(SemanticAttribute * attr);

    public:

      Container();
      Container(attribs_t * attr_head);
      string typeName;
      string formatName;
      string createEvent;
      string destroyEvent;
      bool triggerParent,destroyChildren;
      string triggerEvent;
      void * eventTypes;

      string toPaje();

      const string toString()const;

  }   ;



}



extern attribs_t * attributes;
extern hierarchy_t * toplevel_hierarchy;


SemanticAttribute * new_semantic_attribute();

//void reparent_containers(Paje::Container * c, attribs_t::iterator it);

void init_desc_parser();
hierarchy_t * attr_to_container_hierarchy(attribs_t * attr, hierarchy_t *top);

template <typename CharT, typename Traits>
basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& out, const Paje::Container& r)
{
 return out<< r.toString();
}

template <typename CharT, typename Traits>
basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& out, const SemanticAttribute& r)
{
 return out<< r.toString();
}


#endif
