// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 13:57:10 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  event.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02-09-2011 15:23:14 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */


#include "attributes.hh"
#include "paje.hh"
#include "event.hh"
#include "symbols.hh"
#include <assert.h>
#include <sstream>
#include <string>
#include <map>
#include <stack>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "semantic_types.hh"

extern "C" {
#include <rastro.h>
}

using namespace std;
  
set<pair<string,string>> Paje::container_unique_names;








/*******************************************************************************
 * Paje State functions
 ******************************************************************************/

Paje::State::State(string& name, attribs_t * attribs) {
  this->name = name;
  this->type_identifier = "State";

  this->fill_from_attr(attribs);
}

bool Paje::State::do_start(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName, eventValue;
  
  containerName = format_values(eventType->container->formatName, symbols);

  eventValue = format_values(this->formatValue, symbols);

  pajePushState(timestamp,
                containerName,
                eventType->typeName,
                eventValue, out);

  return true;
}

bool Paje::State::do_end(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName;

  containerName = format_values(eventType->container->formatName, symbols);

    pajePopState(timestamp,
        containerName,
        eventType->typeName,
        out);

  return true;
}


void Paje::State::fill_from_attr(attribs_t * attrs)
{
  Paje::BaseEvent::fill_from_attr(attrs);
  return;
}


bool Paje::State::has_ids() const {
  return ((this->start_id != 0) && (this->end_id != 0));
}

void Paje::State::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_START,*base_id);
  (*base_id)++;

  this->set_trigger_id(EVENT_END, *base_id);
  (*base_id)++;
}

/*******************************************************************************
 * Paje::Link class
 ******************************************************************************/

Paje::Link::Link(string &name, attribs_t * attribs) {
  this->name = name;

  this->fill_from_attr(attribs);
}


string Paje::Link::toString() {
  stringstream out;

  out << Paje::BaseEvent::toString();
  out << "   " << "Key format: " << format_key << endl;

  return out.str();
}




/////
// Trigger functions
bool Paje::Link::do_start(double timestamp,
          symbols_table_t * symbols, ostream &out) {

  LinkType * lt = (LinkType * ) this->eventType;

  string thisContainer, sourceContainer, thisValue, key;

  thisContainer = format_values(lt->container->formatName, symbols);
  sourceContainer = format_values(lt->source->formatName, symbols);
  key  = format_values(this->format_key, symbols);

  pajeStartLink(timestamp,
                   thisContainer,
                   lt->typeName,
                   sourceContainer,
                   thisValue,
                   key,
                   out);

  return false;
}

bool Paje::Link::do_end(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  LinkType * lt = (LinkType * ) this->eventType;

  string thisContainer, destContainer, thisValue, key;

  thisContainer = format_values(lt->container->formatName, symbols);
  destContainer = format_values(lt->source->formatName, symbols);
  key  = format_values(this->format_key, symbols);
  pajeEndLink(timestamp,
                 thisContainer,
                 lt->typeName,
                 destContainer,
                 thisValue,
                 key,
                 out);
  return false;
}

bool Paje::Link::do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  // Actually, Link should call PajeLink....
  cerr << "Error: Class Link has no trigger action" << endl;
  return false;
}

void Paje::Link::fill_from_attr(attribs_t * attrs) {
  Paje::BaseEvent::fill_from_attr(attrs);
  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
        SemanticAttribute * attr = n->getVal();
        switch (attr->id) {
          case ID_LINK:
            break; // name has already been set
          case ID_KEY_FORMAT:
            this->format_key = attr->vals.name;
            break;
          default:
            break;
        }


        return false;
      });
}


bool Paje::Link::has_ids() const {
  return ((this->start_id != 0) && (this->end_id != 0));
}

void Paje::Link::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_START,*base_id);
  (*base_id)++;

  this->set_trigger_id(EVENT_END, *base_id);
  (*base_id)++;
}





//******************************************************************************
// Class ContainerTrigger
// called to execute PajeCreateContainer/PajeDestroyContainer when the
// correct id is seen in the rastro


Paje::ContainerCreateTrigger::ContainerCreateTrigger(Paje::Container * c,hierarchy_t * n)
{
  this->container = c;
  this->hierarchy = n;
  push_timestamp(CONTAINER_CREATE_PRIO); // can leave the EVENT_DEFAULT_PRIO on the stack, won't hurt
}


bool Paje::ContainerCreateTrigger::do_start(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName, parentName;
  
  if (this->container->parent)
    containerName =  format_values(container->parent->formatName, symbols);
  else
    containerName = PAJE_ROOT_CONTAINER;


  walk_tree_head_first(hierarchy,[&](hierarchy_t * n, int level)
      {
        Paje::Container * c = n->getVal();
        parentName = containerName;

        //cerr << "Creating " << c->typeName << endl;

        if (c->triggerParent || (c == this->container) ) {
          // do this for the current container and all it's create on parent children
          containerName = format_values(c->formatName, symbols);

          // do we already have a container with this name, of this type?
          pair<string,string> p(containerName,c->typeName);
          if (Paje::container_unique_names.count(p))
          {
              return true; // ignore
              // can happen, if a container is triggered by event1, and
              // event1 happens several times
          }

          Paje::container_unique_names.insert(p);

          pajeCreateContainer(timestamp,
                                 containerName,
                                 c->typeName,
                                 parentName,
                                 containerName,
                                 out);
          return false;
        }
        return true; // does not create on parent, this one and it's children shall not be created here, there will be an event for it

      });

  return true;
}


bool Paje::ContainerCreateTrigger::do_end(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName;
  
  // needs to destroy all children
  walk_tree_depth_first(this->hierarchy,[&](hierarchy_t * n, int level) 
      {
        Paje::Container * c = n->getVal();
        containerName = format_values(c->formatName, symbols);
        pajeDestroyContainer(timestamp,
                              c->typeName,
                              containerName,
                              out);
        return false;
      }
    );

  return true;
}


void Paje::ContainerCreateTrigger::push_timestamp(double timestamp)
{
  timestamp_stack.push(CONTAINER_DESTROY_PRIO);
}


