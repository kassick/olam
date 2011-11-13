// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Dom, 13 Nov 2011 02:38:16 -0200 (kassick)"
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



#include "event.hh"
#include "paje_functions.hh"


using namespace std;




Paje::Event::Event(const string &name, attribs_t * attribs) {
  this->name = name;

  this->fill_from_attr(attribs);
}


bool Paje::Event::do_start(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) 
{
  return true;
}



bool Paje::Event::do_end(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) 
{
  return true;
}

bool Paje::Event::do_trigger(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) {

  // here should there be something important....
  // like creating an event!!! (I'll just ignore it for now)
  *priority = 0;
  return true;
}

double Paje::Event::get_priority(event_id_t evt_id, double timestamp, symbols_table_t ** symbols)
{
  return 0;
}



bool Paje::Event::has_ids() const {
  return ((this->trigger_id != 0));
}

void Paje::Event::gen_auto_ids(long int * base_id, set<event_id_t> & unique_ids)
{
  Paje::BaseEvent::gen_auto_ids(base_id,unique_ids);
  unique_ids.insert(*base_id);
  this->set_trigger_id(EVENT_TRIGGER,*base_id);
}

