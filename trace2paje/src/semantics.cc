// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.cc"
// Created: "Seg, 01 Ago 2011 15:34:08 -0300 (kassick)"
// Updated: "Seg, 19 Set 2011 19:31:07 -0300 (kassick)"
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

attribs_t * late_parse_tree; // treee to hold ids, values and types statements that must be parsed AFTER all the rest of the parsing


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
      break;
          
    case ID_EVENT_START:
      s << "Event ``" << vals.name << "'' with id ";
      break;
    case ID_STATE_START:
      s << "State ``" << vals.name << "'' starts on ";
      break;
    case ID_STATE_END:
      s << "State ``" << vals.name << "'' ends on   ";
      break;
    case ID_EVENT_ID:
      s << "Eventy Id: " << vals.event_id;
      break;
    case ID_NOP:
      s << "NOP";
      break;

    default:
      s << "Unknown evt " << id;
      break;

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


Paje::Event * get_event_or_warn(char * evt_name) {
  if (! event_names->count(evt_name) ) {
    cerr << "Warning: Event " << evt_name << "has id but no definition, ignoring" << endl;
    return NULL;
  }

  Paje::Event * evt = (*event_names)[evt_name];
  return evt;
}



void parse_late_tree()
{
  Paje::event_id_t evt_id = 0;

  //ids are on the leaves; their parents are always the event/state name
  //do a deep search, whenever you get an event_start, state_start or
  //state_end then there is already an evt_id set up :)
  walk_tree_depth_first(late_parse_tree, [&](attribs_t *n, int level)
      {
        SemanticAttribute * attr = n->getVal();
        Paje::Event *evt;

        switch (attr->id) {
          case ID_EVENT_START:
            if (!(evt = get_event_or_warn(attr->vals.name))) {
              return false; // ignore this one
            }
            evt->set_trigger_id(EVENT_TRIGGER, evt_id);
            break;

          case ID_STATE_START:
            if (!(evt = get_event_or_warn(attr->vals.name))) {
              return false; // ignore this one
            }
            evt->set_trigger_id(EVENT_START, evt_id);
            break;

          case ID_STATE_END:
            if (!(evt = get_event_or_warn(attr->vals.name))) {
              return false; // ignore this one
            }
            evt->set_trigger_id(EVENT_END, evt_id);
            break;

          case ID_EVENT_ID:
            evt_id = attr->vals.event_id;
            break;
          case ID_NOP:
            break;
          default:
            // ignore
            cerr << "ignored evt:" << attr->vals.event_id << endl;
            break;
         }

        return false; // not found -- aka go all over the tree


      } );
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

  SemanticAttribute * attr = new SemanticAttribute();
  attr->id = ID_NOP;
  late_parse_tree = new attribs_t(attr);
}


