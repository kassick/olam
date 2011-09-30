// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Sex, 30 Set 2011 18:36:22 -0300 (kassick)"
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
#include <set>
#include <stack>
#include <list>
#include <iostream>
#include <limits>
#include "rastro_helper.hh"
extern "C" {
#include <rastro.h>
}



#define DEFAULT_EVENT_PRIO (-(numeric_limits<double>::max()) )
#define CONTAINER_CREATE_PRIO (-10)
#define CONTAINER_DESTROY_PRIO (-1)



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


  extern set<pair<string,string>> container_unique_names;

  //************************************************
  //Class: Paje::EventType
  class EventType: public PajeElement {
    public:
      string typeName;
      Container * container;

      EventType(const string &typeName);
      EventType(const string &typeName,Paje::Container * c);
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

      string type_identifier;

      event_id_t  start_id, end_id, trigger_id;

      identifier_list_t identifier_names;
      stack<double> timestamp_stack;

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

  
      bool load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols);
      void add_symbol_from_tree(attribs_t * attrs);

      virtual string toString();

      bool operator<(const Event * e) const;

      virtual void push_timestamp(const double timestamp);
      virtual double pop_timestamp(const double timestamp);
      virtual double get_priority() const;

  };


  class DummyEvent: public Event {
      public: 
        DummyEvent(EventType * evt_type);
        virtual bool do_start(double timestamp,
            symbols_table_t * symbols, ostream &out);

        virtual bool do_end(double timestamp,
            symbols_table_t * symbols, ostream &out);

        virtual bool do_trigger(double timestamp,
            symbols_table_t * symbols, ostream &out);
  };



  class State: public Event {
    public:

      State(string &_name, attribs_t * attrs) ;

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);
      
      virtual void fill_from_attr(attribs_t * attrs);
  } ;



  class Link: public Event {
    public:
      string format_key;
      Link(string &name, attribs_t * attribs);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);
      
      virtual bool do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out);
      
      virtual void fill_from_attr(attribs_t * attrs);
      virtual string toString() ;
  };



  class ContainerCreateTrigger: public Event {
    public:
      ContainerCreateTrigger(Paje::Container * c, hierarchy_t * n);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);
      

      virtual void push_timestamp(double timestamp);

    protected:
      Container * container;
      hierarchy_t * hierarchy;

  };



} // Paje namespace


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




// sadly, #SHORT_NAME[0] does not work as the compiler does not recognize
// str_const[int_const] as a char in compilation time
#define CASE_TYPE1(CHAR_NAME,SHORT_NAME,VALUE) \
  case CHAR_NAME : \
    VALUE = SHORT_NAME;   \
    break;                \


#endif
