// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Seg, 19 Set 2011 19:29:57 -0300 (kassick)"
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
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <rastro.h>

using namespace std;

/*******************************************************************************
 Paje::Event functions

*******************************************************************************/

/////
// Trigger functions
bool Paje::Event::do_start(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  cerr << "Error: Class Event has no start action" << endl;
  return false;
}

bool Paje::Event::do_end(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  cerr << "Error: Class Event has no end  action" << endl;
  return false;
}

bool Paje::Event::do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  // Actually, Event should call PajeEvent....
  cerr << "Error: Class Event has no trigger action" << endl;
  return false;
}




/////
//Sets the trigger field
void Paje::Event::set_trigger_id(trigger_id_t trigger_field, event_id_t id)
{
  switch (trigger_field) {
    case EVENT_TRIGGER:
      this->trigger_id = id;
      break;
    case EVENT_START:
      this->start_id = id;
      break;
    case EVENT_END:
      this->end_id = id;
      break;
    default:
      cerr << "Invalid value at set_trigger_id: " << trigger_id << endl;
      break;
  }
}

////
//Calls the correct function based on the id

bool Paje::Event::trigger(event_id_t evt_id, double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  if (evt_id == trigger_id)
    return do_trigger(timestamp,symbols,out);

  if (evt_id == start_id)
    return do_start(timestamp,symbols,out);

  if (evt_id == end_id)
    return do_end(timestamp,symbols,out);
}

///
void Paje::Event::fill_from_attr(attribs_t * attrs)
{
  cerr << "Abstract class, dumbass!" << endl;
}


////
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



/*******************************************************************************
 * Paje State functions
 ******************************************************************************/

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
  return;
}
