// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Ter, 06 Set 2011 15:42:28 -0300 (kassick)"
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


#include "paje.hh"
#include "event.hh"
#include "symbols.hh"
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <rastro.h>

using namespace std;

bool Paje::Event::trigger(event_id_t evt_id, double timestamp,
          symbols_table_t * symbols, ostream &out)
{
  out << "Abstract class, nothing to see here" << endl;
  return false;
}


bool Paje::StateEvent::trigger(event_id_t evt_id, double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName, eventValue;
  
  containerName = format_values(eventType->container->formatName, symbols);


  if (evt_id == start_id) {
    eventValue = format_values(this->formatValue, symbols);

    pajePushState(timestamp,
                  containerName,
                  eventType->typeName,
                  eventValue, out);

  } else {
    pajePopState(timestamp,
        containerName,
        eventType->typeName,
        out);
  }

  return true;
}

bool Paje::Event::load_symbols(rst_event_t *event, symbols_table_t * symbols)
{
  int count_c, count_w, count_i, count_l, count_f, count_d, count_s;
  Paje::Symbol *symbol;


  Paje::identifier_list_t::iterator it;
  for (it = identifier_names.begin(); it != identifier_names.end(); ++it)
  {
    string &field_name = (*it).field_name;
    // automatically creates an entry and allocate symbol, I don't need to
    // worry yay! \o/
    Paje::Symbol &symbol = (*symbols)[field_name];

    switch ((*it).type) {
      CASE_TYPE( c , uint8);
      CASE_TYPE( w , uint16);
      CASE_TYPE( i , uint32);
      CASE_TYPE( l , uint64);
      CASE_TYPE( f , float);
      CASE_TYPE( d , double);
      CASE_TYPE( s , string);

      default:
        cerr << "Unknown type !?!?!? " <<endl;
        return false;
        break;
    }

  }

  (*symbols)["EVT_NAME"].set_value(this->name.c_str());

}
