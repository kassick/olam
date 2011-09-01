// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Sex, 02 Set 2011 18:05:52 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  event.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  03-08-2011 16:14:50 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __EVENT_H__
#define __EVENT_H__

#include "paje.hh"
#include "event.hh"
#include "symbols.hh"
#include "container.hh"
#include <sstream>
#include <string>
#include <map>
#include <list>
#include <iostream>
#include <rastro.h>

namespace Paje {

  typedef unsigned int event_id_t;

  typedef struct _identifier_entry_t {
    string field_name;
    Paje::basic_types_t type;
  } identifier_entry_t;

  typedef list<identifier_entry_t> identifier_list_t;

  class EventType {
    public:
      string typeName;
      Container * container;
  };

  class Event {
    public:
      EventType *eventType;

      string formatValue;
      string name;

      identifier_list_t identifier_names;

      bool trigger(event_id_t evt_id, double timestamp,
          symbols_table_t * symbols, ostream &out);

  
     bool load_symbols(rst_event_t *event, symbols_table_t * symbols);

  };


  class StateEvent : Event {
    public:
      event_id_t start_id, end_id;

      StateEvent(event_id_t _start_id, event_id_t _end_id): start_id{_start_id}, end_id{_end_id} {};

      bool trigger(event_id_t evt_id, double timestamp,
          symbols_table_t * symbols, ostream &out);
  } ;
}

typedef std::map<std::string      , Paje::Event *>      event_name_map_t;
typedef std::map<Paje::event_id_t , Paje::Event *>      event_id_map_t;
typedef std::map<std::string      , Paje::EventType *>  event_type_name_map_t;



#define CASE_TYPE(SHORT_NAME , PAJE_NAME) \
        case SHORT_NAME:                                                             \
          if (count_ ## SHORT_NAME < event->ct.n_ ## PAJE_NAME) {              \
            symbol = new Paje::Symbol();                                         \
            symbol->set_value( event->v_ ## PAJE_NAME [ count_ ## SHORT_NAME ] );\
            count_ ## SHORT_NAME ++  ;                                      \
          } else {                                                          \
            cerr << "Not enought values of type " << #SHORT_NAME  <<endl;   \
            return false;                                                   \
          }                                                                 \
          break; // trick?




#endif
