// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/state.cc"
// Created: "Ter, 04 Out 2011 13:59:53 -0300 (kassick)"
// Updated: "Qua, 05 Out 2011 15:12:44 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  state.cc
 *
 *    Description:  Paje State functions
 *
 *        Version:  1.0
 *        Created:  04-10-2011 13:59:53 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */


#include "state.hh"

Paje::State::State(string& name, attribs_t * attribs) {
  this->name = name;
  this->type_identifier = "State";

  this->fill_from_attr(attribs);
}

bool Paje::State::do_start(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  string containerName, eventValue;

  Paje::BaseEvent::do_start(timestamp, symbols, priority, out);
  
  containerName = format_values(eventType->container->formatName, symbols);

  eventValue = format_values(this->formatValue, symbols);

  pajePushState(timestamp,
                containerName,
                eventType->typeName,
                eventValue, out);

  return true;
}

bool Paje::State::do_end(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  string containerName;

  Paje::BaseEvent::do_end(timestamp, symbols, priority, out);
  containerName = format_values(eventType->container->formatName, symbols);

  if (*priority   > 0) // happened somewhere in time
  {

    pajePopState(timestamp,
                  containerName,
                  eventType->typeName,
                  out);
  } else {
    cerr << "[Warning] Can not pop state " << this->name 
         << " in container " << containerName
         << ": no previous state pushed that we can pop"
         << endl;
  }

  return true;
}


void Paje::State::fill_from_attr(attribs_t * attrs)
{
  Paje::BaseEvent::fill_from_attr(attrs);
  return;
}


bool Paje::State::has_ids() const {
  return ((this->start_id != 0) && (this->end_id != 0));
}

void Paje::State::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_START,*base_id);
  (*base_id)++;

  this->set_trigger_id(EVENT_END, *base_id);
  (*base_id)++;
}

