// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/baseevent.hh"
// Created: "Ter, 04 Out 2011 11:50:46 -0300 (kassick)"
// Updated: "Qui, 06 Out 2011 15:49:53 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  baseevent.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 11:50:46 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __BASEEVENT_H__
#define __BASEEVENT_H__

extern "C" {
#include <rastro.h>
}
#include "rastro_helper.hh"

#include "symbols.hh"
#include "attributes.hh"
#include "container.hh"

#include <sstream>
#include <string>
#include <map>
#include <set>
#include <stack>
#include <list>
#include <iostream>
#include <limits>


#include "baseeventtype.hh"

#define DEFAULT_EVENT_PRIO (-(numeric_limits<double>::max()) )
#define DEFAULT_EVENT_PRIO_TRIGG (0)

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
  //Class: Paje::Event
  class BaseEvent: public PajeElement {
    public:
      BaseEventType *eventType;

      string formatValue;
      string name;

      string type_identifier;

      event_id_t  start_id, end_id, trigger_id;

      identifier_list_t identifier_names;
      map<string, stack<double>> timestamp_map;

      BaseEvent();

      virtual bool trigger(event_id_t evt_id, 
          double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);

      virtual void fill_from_attr(attribs_t * attrs);

      virtual void set_trigger_id(trigger_id_t trigger_id, event_id_t id);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, 
          double * priotity,
          ostream &out);


      virtual bool do_end(double timestamp,
          symbols_table_t * symbols,
          double * priotity,
          ostream &out);

      virtual bool do_trigger(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);

      virtual void set_event_type(BaseEventType * evt_type);

  
      virtual bool load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols);
      void add_symbol_from_tree(attribs_t * attrs);

      virtual string toString();

      //bool operator<(const BaseEvent * e) const;

      virtual void push_timestamp(const string & containerName, const double timestamp);
      virtual double pop_timestamp(const string & containerName, const double timestamp);

      //virtual double get_priority(const string & containerName) const;

      virtual bool has_ids() const;

      virtual void gen_auto_ids(long int * base_id);

      virtual void gen_ids_description(ostream & out);

      void get_rst_function_signature(set<string> & signatures);
      void gen_fort_header(ostream & out,
          const string & start_suffix,
          const string & end_suffix,
          const string & evt_suffix);
      void gen_c_header(ostream & out, 
          const string & start_suffix,
          const string & end_suffix,
          const string & evt_suffix);

      virtual bool fits_in_event_type(const BaseEventType* evt_type) const;

    private:
      stringstream rst_function_signature_buf;

  };

} // Paje namespace





#endif
