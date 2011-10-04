// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/containertrigger.cc"
// Created: "Ter, 04 Out 2011 14:07:13 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:25:48 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  containertrigger.cc
 *
 *    Description:  Class ContainerTrigger
 *                  called to execute PajeCreateContainer/PajeDestroyContainer when the
 *                  correct id is seen in the rastro
 *
 *        Version:  1.0
 *        Created:  04-10-2011 14:07:13 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "containertrigger.hh"
#include "paje_functions.hh"
#include <set>

set<pair<string,string>> Paje::container_unique_names;


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


