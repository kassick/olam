// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Qua, 21 Set 2011 21:51:53 -0300 (kassick)"
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

#include "attributes.hh"
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
#include "rastro_helper.hh"

namespace Paje {

  typedef enum _trigger_id_t {
    EVENT_TRIGGER,   // one time event
    EVENT_START,     // do a start action
    EVENT_END,       // do an end action
    EVENT_NOP,       // warn of bad implementation
  } trigger_id_t;


  typedef unsigned int event_id_t;

  typedef struct _identifier_entry_t {
    string field_name;
    rastro_basic_types_t type;
  } identifier_entry_t;

  typedef list<identifier_entry_t> identifier_list_t;

  //************************************************
  //Class: Paje::EventType
  class EventType: public PajeElement {
    public:
      string typeName;
      Container * container;

      EventType(string &typeName);
      EventType(string &typeName,Paje::Container * c);
      //EventType(string &typeName, attribs_t * attribs);

      virtual void do_header(ostream &out);

  };
  
  //************************************************
  //Class: Paje::LinkType
  class LinkType: public EventType {
    public:
      Container *source, *dest;

      LinkType(string& typeName, Paje::Container * c, attribs_t * t);

      virtual void do_header(ostream &out);
  };



  //************************************************
  //Class: Paje::Event
  class Event: public PajeElement {
    public:
      EventType *eventType;

      string formatValue;
      string name;

      event_id_t  start_id, end_id, trigger_id;

      identifier_list_t identifier_names;

      Event();

      virtual bool trigger(event_id_t evt_id, double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual void fill_from_attr(attribs_t * attrs);

      virtual void set_trigger_id(trigger_id_t trigger_id, event_id_t id);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out);

  
     bool load_symbols(rst_event_t *event, symbols_table_t * symbols);

  };


  class State: public Event {
    public:
      event_id_t start_id, end_id;

      State(string &_name, event_id_t _start_id, event_id_t _end_id):
        start_id{_start_id}, end_id{_end_id} 
      {
        this->name = name;
      };

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);
      
      void fill_from_attr(attribs_t * attrs);
  } ;
}




#define CASE_TYPE(SHORT_NAME , PAJE_NAME) \
        case SHORT_NAME:                                                             \
          if (count_ ## SHORT_NAME < event->ct.n_ ## PAJE_NAME) {              \
            symbol.set_value( event->v_ ## PAJE_NAME [ count_ ## SHORT_NAME ] );\
            count_ ## SHORT_NAME ++  ;                                      \
          } else {                                                          \
            cerr << "Not enought values of type " << #SHORT_NAME  <<endl;   \
            return false;                                                   \
          }                                                                 \
          break; // trick?




#endif
