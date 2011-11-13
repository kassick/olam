// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/containertrigger.hh"
// Created: "Ter, 04 Out 2011 14:07:16 -0300 (kassick)"
// Updated: "Dom, 13 Nov 2011 02:35:19 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  containertrigger.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 14:07:16 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __CONTAINERTRIGGER_H__
#define __CONTAINERTRIGGER_H__


#define CONTAINER_CREATE_PRIO (-10)
#define CONTAINER_DESTROY_PRIO (-1)

#include "baseevent.hh"
#include "container.hh"
#include <iostream>
#include <string>
#include <map>



namespace Paje {

  typedef struct _unique_container_name_t {
    string containerName;
    string typeName;
    string parentName;
    Paje::Container * container;
    unsigned long nchild;

    _unique_container_name_t() {
      nchild = 0;
    }

    /*
    bool operator<(const struct _unique_container_name_t &other) const
    {
      return ( (this->containerName < other.containerName) );
    } */

    /*
    bool operator==(const struct _unique_container_name_t &other) const
    {
      return (( this->containerName == other.containerName) );
    }*/

    /*

    struct _unique_container_name_t & operator=(const struct _unique_container_name_t &other)
    {
      this->containerName = other.containerName;
      this->typeName = other.typeName;
      this->parentName = other.parentName;
      this->nchild = other.nchild;
      return *this;

    }*/

  } unique_container_name_t;

  extern map<string,unique_container_name_t> container_unique_names;



  class ContainerCreateTrigger: public BaseEvent {
    public:
      ContainerCreateTrigger(Paje::Container * c, hierarchy_t * n);

      virtual bool do_start(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out);
     
      virtual double get_priority(event_id_t evt_id, double timestamp, symbols_table_t ** symbols);

      // need to overload these so that we load from the right event -- the
      // create and destroy events
      virtual bool load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols);
      virtual void push_symbols(event_id_t id,
                                   symbols_table_t * from,
                                   symbols_table_t * to);
      
      virtual void map_symbols(Paje::event_id_t evt_id, symbols_table_t ** symbols, user_defined_maps_t & usermaps);
      virtual void get_symbols_from_map (Paje::event_id_t evt_id, symbols_table_t * dest_symbols, symbols_table_t ** symbols, user_defined_maps_t & usermaps);

    protected:
      Container * container;
      hierarchy_t * hierarchy;

  };



} // Paje namespace

void destroy_missing_containers(double timestamp,ostream &out);

#endif
