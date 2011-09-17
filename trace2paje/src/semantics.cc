// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.cc"
// Created: "Seg, 01 Ago 2011 15:34:08 -0300 (kassick)"
// Updated: "Sex, 16 Set 2011 16:24:59 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  semantics.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01-08-2011 15:34:08 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "paje.hh"
#include "tree.hpp"
#include <iostream>
#include <map>
#include <sstream>
#include "attributes.hh"
#include "container.hh"
#include "semantics.hh"


//Globals
attribs_t *attributes;
hierarchy_t *toplevel_hierarchy;
container_type_names_t * container_type_names;

event_type_name_map_t * eventtype_names;
event_name_map_t      * event_names;
event_id_map_t        * event_ids;
queue<string> files_to_parse;


const string SemanticAttribute::toString() const {
  stringstream s;
  switch (id) {
    case ID_NAME:
      s << "FormatName: " <<vals.format_name;
      break;
    case ID_CREATE_EVENT:
      s << "Create on Event "<<vals.create_event;
      break;
    case ID_CREATE_PARENT:
      if (vals.create_parent)
        s << "Create on Parent";
      else
        s << "Do not create on parent -- huh!?";
      break;
    case ID_DESTROY_CHILDREN:
      if (vals.destroy_children)
        s << "Destroy on Children";
      else
        s << "Do not destroy on children -- huh!?";
      break;
    case ID_DESTROY_EVENT:
      s << "Destroy on Event " << vals.destroy_event;
      break;
    case ID_CONTAINER:
      s << "Container Type " << vals.container_type;
      break;

    case ID_ACCEPT_LIST:
      s << "Accept event " << vals.identifier_name;
      break;

    case ID_EVENT_TYPE:
      s << "Event type " << vals.name;

  }

  return s.str();
}

//Utility functions

// Create a new semantic attribute
SemanticAttribute *new_semantic_attribute()
{
  return new SemanticAttribute();
}


// Create a tree of Paje::Container* from a tree of SemanticAttribute
hierarchy_t * attr_to_container_hierarchy(attribs_t * attr, hierarchy_t *top)
{
  hierarchy_t * top_, *h;
  attribs_t::iterator it;

  h = NULL;

  if (!attr) return NULL;

  if (attr->getVal()->id == ID_CONTAINER) {
    Paje::Container * c;
    char * container_type = attr->getVal()->vals.container_type;


    c = new Paje::Container(container_type, attr);

    top_ = new hierarchy_t(c);
    top->addChild(top_);

    h = top_;
  } else {
    top_ = top;
  }

  for(it = attr->begin(); it != attr->end(); ++it)
  {
    hierarchy_t * child = attr_to_container_hierarchy(*it,top_);
    /*if (child) {
      top_->addChild(child);
    }*/
  }

  return h;

}





// Fill in container ids and check for unique types
void check_unique_container_types()
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level) {
        Paje::Container * c = h->getVal();
        if (container_type_names->find(c->typeName) != container_type_names->end())
        {
          cerr << "Error: A container with type " << c->typeName << " has already been defined" <<endl;
          return true;
        }
        (*container_type_names)[c->typeName] = h;
        return false;
      }) )
  {
    exit(1);
  }
}

// Fill in container ids and check for unique types
void check_unique_event_types()
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level)
      {
        Paje::Container * c = h->getVal();


        list<string>::iterator it;
        for (it = c->event_types.begin(); it != c->event_types.end(); ++it)
        {
          if (eventtype_names->find(*it) != eventtype_names->end())
          {
            cerr << "Error: An event type named " << *it << " has already been defined" <<endl;
            return true;
          }

          EventType &t = (*eventtype_names)[*it];
          t.typeName = *it;
          t.container = c;
        }
        return false;
      }) )
  {
    exit(1);
  }
}


void check_unique_types() 
{
  check_unique_container_types();
  check_unique_event_types();
}



// walk the hierarchy dump paje formated output
void hierarchy_to_paje(ostream &out)
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level) {
        Paje::Container * c = h->getVal();
        Paje::Container * parent;

        if (c->typeName == "0")
          return false; // skip container 0

        parent = h->getParent()->getVal();

        pajeDefineContainerType(c->typeName, parent->typeName, c->typeName,out);

        return false;
      })) {
    cerr << "Error while dumping paje hierarchy. What on earth!? " <<endl;
  }
}

// Walk the hierarchy and dump event types
void event_types_to_paje(ostream &out)
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level) {
        Paje::Container * c = h->getVal();

        list<string>::iterator it;
        for (it = c->event_types.begin(); it != c->event_types.end(); ++it)
        {
          pajeDefineStateType(*it, c->typeName, *it, out);
        }

        return false;
      })) {
    cerr << "Error while dumping paje event types. What on earth!? " << endl;
  }
}











// Initialization Function -- Parser wide

void init_desc_parser()
{
  using namespace Paje;
  Paje::Container * zero;
  

  zero = new Container("0");

  toplevel_hierarchy = new TreeNode < Paje::Container * >(zero);
  attributes = new TreeNode < SemanticAttribute *>(NULL);
  container_type_names = new container_type_names_t();// must happen before creating a container

  eventtype_names = new event_type_name_map_t();
  event_names     = new event_name_map_t();
  event_ids       = new event_id_map_t();
}


