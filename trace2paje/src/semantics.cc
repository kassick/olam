// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.cc"
// Created: "Seg, 01 Ago 2011 15:34:08 -0300 (kassick)"
// Updated: "Qua, 21 Set 2011 21:53:32 -0300 (kassick)"
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
#include <algorithm>
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

attribs_t * late_parse_tree, // tree to hold ids, values and types statements that must be parsed AFTER all the rest of the parsing
          * early_parse_tree; // three that holds the hierarchy



///************************************************************************************
// Class functions
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



//******************************************************************************
//Utility functions

// ****************************************
// Create a new semantic attribute
SemanticAttribute *new_semantic_attribute()
{
  return new SemanticAttribute();
}





// **********************************************************************
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
    if (container_type_names->count(container_type))
    {
      cerr << "Error: A container with type " << container_type << " has already been defined" <<endl;
      exit(1);
    }



    c = new Paje::Container(container_type, attr);

    c->parent = top->getVal(); // the tree takes care of this, but do_header needs to know who is the parent... matter of isolation

    top_ = new hierarchy_t(c);
    (*container_type_names)[container_type] = top_;
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


void attr_to_event_types(attribs_t * attribs)
{
  string container_name;
  walk_tree_head_first(attribs,[&](attribs_t * n, int level) {
    SemanticAttribute * attr = n->getVal();
    if (attr->id == ID_CONTAINER)
    {
      container_name = attr->vals.name;
    } else if (attr->id == ID_EVENT_TYPE) 
    {
      if (eventtype_names->count(attr->vals.name))
      {
        cerr << "Error: EventType " << attr->vals.name << "has already been defined" << endl;
        exit(1);
      }

      if (container_type_names->count(container_name) == 0) {
        cerr << "Error: Event Type " << attr->vals.name << "can not find it's container type " << container_name << endl;
        exit(1);
      }

      Paje::Container * c = (*container_type_names)[container_name]->getVal();
      string tn = attr->vals.name;

      Paje::EventType * evttype = new Paje::EventType(tn,c);
      (*eventtype_names)[attr->vals.name] = evttype;

      return true; // can stop this subtree
    }
    return false;
  });
}

void attr_to_link_types(attribs_t * attribs)
{
  string container_name;
  walk_tree_head_first(attribs,[&](attribs_t * n, int level) {
    SemanticAttribute * attr = n->getVal();
    if (attr->id == ID_CONTAINER)
    {
      container_name = attr->vals.name;
    } else if (attr->id == ID_LINK_TYPE) 
    {
      if (eventtype_names->count(attr->vals.name))
      {
        cerr << "Error: LinkType " << attr->vals.name << "has already been defined" << endl;
        exit(1);
      }


      if (container_type_names->count(container_name) == 0) {
        cerr << "Error: Link Type " << attr->vals.name << "can not find it's container type " << container_name << endl;
        exit(1);
      }

      Paje::Container * c = (*container_type_names)[container_name]->getVal();
      string tn = attr->vals.name;

      Paje::LinkType * evttype = new Paje::LinkType(tn,c,n);
      (*eventtype_names)[attr->vals.name] = evttype;

      return true; // can stop this subtree
    }
    return false;
  });
}

void map_accept_attrs(attribs_t * attribs)
{
  string container_name;
  string event_type;
  walk_tree_head_first(attribs,[&](attribs_t * n, int level) {
    SemanticAttribute * attr = n->getVal();
    
    if (attr->id == ID_EVENT_TYPE)
    {
      event_type = attr->vals.name;
      //cerr << "event type " << event_type << endl;
      attribs_t::iterator it;
      for(it = n->begin(); it != n->end(); ++it)
      {
        SemanticAttribute *list_attr = (*it)->getVal();
        if (list_attr->id == ID_ACCEPT_LIST) {
          walk_tree_head_first(*it,[&](attribs_t * n1, int level1) {
              SemanticAttribute * attr1 = n1->getVal();
              if (attr1->id == ID_ACCEPT_LIST) {
                //cerr << "  accepts " << attr1->vals.name << endl;
                if (!event_names->count(attr1->vals.name)) {
                  cerr << "Warning: Event type " << event_type << " accepts undefined event " << attr1->vals.name << endl;
                } else {
                  Paje::Event * evt = (*event_names)[attr1->vals.name];
                  evt->eventType = (*eventtype_names)[event_type];
                }
              }
              return false; // visit all the children of event_type
            });
        }
      }
    }
    return false;
  });
}

void check_events_have_type()
{
  for_each(event_names->begin(), event_names->end(), [&](pair<string, Paje::Event *> p) {
      if (p.second->eventType == NULL)
      {
        cerr << "Error: Event " << p.first << " has no defined type" << endl;
        exit(1);
      }
    });
}



#if 0
// ************************************************************
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

#endif





#if 0
// ****************************************
// Fill in container ids and check for unique types
void check_unique_event_types()
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level)
      {
        Paje::Container * c = h->getVal();


        list<string>::iterator it;
        for (it = c->event_types.begin(); it != c->event_types.end(); ++it)
        {
          if (eventtype_names->count(*it))
          {
            cerr << "Error: An event type named " << *it << " has already been defined" <<endl;
            return true;
          }

          EventType *t = new EventType();
          t->typeName = *it;
          t->container = c;
          (*eventtype_names)[*it] = t;
        }
        return false;
      }) )
  {
    exit(1);
  }
}
#endif

#if 0
void check_unique_types() 
{
  check_unique_container_types();
  check_unique_event_types();
}
#endif



//************************************************************
// walk the hierarchy dump paje formated output
void hierarchy_to_paje(ostream &out)
{
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level) {
        Paje::Container * c = h->getVal();
        Paje::Container * parent;

        if (c->typeName == "0")
          return false; // skip container 0

        parent = h->getParent()->getVal();

        //pajeDefineContainerType(c->typeName, parent->typeName, c->typeName,out);
        c->do_header(out);

        return false;
      })) {
    cerr << "Error while dumping paje hierarchy. What on earth!? " <<endl;
  }
}

// Walk the hierarchy and dump event types
void event_types_to_paje(ostream &out)
{
#if 0
  if (walk_tree_head_first(toplevel_hierarchy, [&](hierarchy_t * h, int level) {
        Paje::Container * c = h->getVal();

        list<string>::iterator it;
        for (it = c->event_types.begin(); it != c->event_types.end(); ++it)
        {
          //pajeDefineStateType(*it, c->typeName, *it, out);
          if (! eventtype_names->count(*it) ) {
            cerr << "Internal Error: no event type with name " << *it << endl;
            exit(1);
          }
          Paje::EventType * t = (*eventtype_names)[*it];
          t->do_header(out);
        }

        return false;
      })) {
    cerr << "Error while dumping paje event types. What on earth!? " << endl;
  }
#endif
  for_each(eventtype_names->begin(), eventtype_names->end(), 
      [&](pair<string,Paje::EventType *> p) {
        p.second->do_header(out);
      });
}


//****************************
//Get an event name or wanr that it does not exist
Paje::Event * get_event_or_warn(char * evt_name) {
  if (! event_names->count(evt_name) ) {
    cerr << "Warning: Event " << evt_name << "has id but no definition, ignoring" << endl;
    return NULL;
  }

  Paje::Event * evt = (*event_names)[evt_name];
  return evt;
}


//*****************************
//Parse the ID, VALUE and TYPE tokens after all event types have been
//already defined

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





//****************************************
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

  SemanticAttribute * attr1 = new SemanticAttribute();
  SemanticAttribute * attr2 = new SemanticAttribute();
  attr1->id = ID_NOP;
  attr2->id = ID_NOP;
  late_parse_tree = new attribs_t(attr1);
  early_parse_tree = new attribs_t(attr2);
}


