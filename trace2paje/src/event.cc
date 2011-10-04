// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:06:01 -0300 (kassick)"
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


