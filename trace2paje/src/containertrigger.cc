// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/containertrigger.cc"
// Created: "Ter, 04 Out 2011 14:07:13 -0300 (kassick)"
// Updated: "Seg, 10 Out 2011 17:20:31 -0300 (kassick)"
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
#include "semantic_types.hh"
#include "state.hh"
#include <map>

map<string,Paje::unique_container_name_t> Paje::container_unique_names;


Paje::ContainerCreateTrigger::ContainerCreateTrigger(Paje::Container * c,hierarchy_t * n)
{
  this->container = c;
  this->hierarchy = n;
}

bool Paje::ContainerCreateTrigger::load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols)
{
  if ( (id == this->end_id) && (event_names->count(this->container->destroyEvent)) )
  {
    // When destroying an event, load the symbols from it's associated
    // event
    Paje::BaseEvent * destroy_evt = (*event_names)[this->container->destroyEvent];
    destroy_evt->load_symbols(destroy_evt->end_id, event, symbols);

    //Paje::BaseEvent * create_evt  = (*event_names)[this->container->createEvent];
  }
  
  if ( (id == this->start_id) && (event_names->count(this->container->createEvent)) )
  {
    Paje::BaseEvent * create_evt = (*event_names)[this->container->createEvent];
    create_evt->load_symbols(create_evt->start_id, event, symbols);
  }

}






bool Paje::ContainerCreateTrigger::do_start(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  string containerName, parentName;

  // forcefully load symbols from it's create event
  

  if (this->container->parent)
    containerName =  format_values(container->parent->formatName, symbols);
  else
  { // Should not happen, 0 is not associated with any event and every container needs to be child of 0 (or someone under it)
    containerName = PAJE_ROOT_CONTAINER;
  }

  walk_tree_head_first(hierarchy,[&](hierarchy_t * n, int level)
      {
        Paje::Container * c = n->getVal();
        parentName = containerName;

        //cerr << "Creating " << c->typeName << endl;

        if (c->triggerParent || (c == this->container) ) {
          // do this for the current container and all it's create on parent children
          containerName = format_values(c->formatName, symbols);

          // do we already have a container with this name, of this type?

          if (Paje::container_unique_names.count(containerName))
          {
              /*cerr << "already contains "
                    << p.containerName << ","
                    << p.parentName << ","
                    << p.typeName << ")"  << endl; */
              return false; // ignore
              // can happen, if a container is triggered by event1, and
              // event1 happens several times
          }
          
          //Add the name to the map
          unique_container_name_t *p = &(Paje::container_unique_names[containerName]);

          p->containerName = containerName;
          p->typeName      = c->typeName;
          p->parentName    = parentName;
          p->nchild = 0;
            
          /*
          cerr << "Inserted? " << containerName << "("
                    << p->containerName << ","
                    << p->parentName << ","
                    << p->typeName << ","
                    << p->nchild << ")"  
                    << endl;
                    */


          pajeCreateContainer(timestamp,
                                 containerName,
                                 c->typeName,
                                 parentName,
                                 containerName,
                                 out);

          // increase the number of children of the parent
          p = &(Paje::container_unique_names[parentName]);
          p->nchild++;
          /*
          cerr << "incremented? " << parentName << "("
                    << p->containerName << ","
                    << p->parentName << ","
                    << p->typeName << ","
                    << p->nchild << ")"  
                    << endl;
                    */


          return false;
        }
        return true; // does not create on parent, this one and it's children shall not be created here, there will be an event for it

      });

  *priority = CONTAINER_CREATE_PRIO;

  return true;
}


bool Paje::ContainerCreateTrigger::do_end(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  string containerName;
  
  // needs to destroy all children
  walk_tree_depth_first(this->hierarchy,[&](hierarchy_t * n, int level)
      {
        Paje::Container * thatcontainer = n->getVal();
        string parentName = format_values(thatcontainer->parent->formatName,symbols);

       
        /*
        cerr << "Container type " << thatcontainer->typeName << "child of " 
              << thatcontainer->parent->typeName << "(" << parentName << ")"
              << endl;
              */

        map<string,Paje::unique_container_name_t>::iterator it;

        it = container_unique_names.begin();
        while (it != container_unique_names.end())
        {

        /*
            cerr << "Kill? ("
                    << it->second.containerName << ","
                    << it->second.parentName << ","
                    << it->second.typeName << ","
                    << it->second.nchild << ")"  
                    << endl;
                    */

            if ( (thatcontainer->typeName == it->second.typeName) &&
                 (parentName == it->second.parentName) &&
                 (it->second.nchild == 0 ) )
            {
              close_pending_states(timestamp,it->second.containerName, thatcontainer, out);
              pajeDestroyContainer(timestamp,
                                    it->second.typeName,
                                    it->second.containerName,
                                    out);

              // tell parent it has less one child
              Paje::unique_container_name_t &p = container_unique_names[parentName];
              if (!p.nchild)
                cerr << "[Warning] Internal Error: Container " 
                     << parentName 
                     << " already has 0 children; how come are we destroying it yet?" 
                     << endl;
              else
                p.nchild--;
              
              it = container_unique_names.erase(it);
            } else 
              ++it;
        }
        return false;
      }
    );





  *priority = CONTAINER_DESTROY_PRIO;
  return true;
}




void destroy_missing_containers(double timestamp, ostream & out)
{
  walk_tree_depth_first(toplevel_hierarchy,
      [&timestamp,&container_unique_names,&out](hierarchy_t * n, int level) 
      {
        Paje::Container * thisContainer = n->getVal();

        if (thisContainer->parent == NULL)
          return true;  // got to the parent

        /*
        cerr << "Container type " << thatcontainer->typeName << "child of " 
              << thatcontainer->parent->typeName << "(" << parentName << ")"
              << endl;*/

        map<string,Paje::unique_container_name_t>::iterator it;

        it = container_unique_names.begin();
        while (it != container_unique_names.end())
        {
            // pair of containername, typename

            if ( (thisContainer->typeName == it->second.typeName))
            {
              cerr << "[Warning] Destroying missing container " << it->second.containerName << endl;
              pajeDestroyContainer(timestamp,
                                    it->second.typeName,
                                    it->second.containerName,
                                    out);
              
              
              Paje::unique_container_name_t &p = container_unique_names[it->second.parentName];
              if (!p.nchild)
                cerr << "[Warning] Internal Error: Container " 
                     << it->second.parentName 
                     << " already has 0 children; how come are we destroying it yet?" 
                     << endl;
              else
                p.nchild--;
              
              it = container_unique_names.erase(it);
            } else 
              ++it;
        }
        return false;
      }
    );
}
