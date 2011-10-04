// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/dummyevent.cc"
// Created: "Ter, 04 Out 2011 13:48:25 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 16:13:32 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  dummyevent.cc
 *
 *    Description:  Paje Dummy Event
 *                  This one is created for referenced events with id but that have no EVENT
 *                  or STATE defined. They act as nicknames for some IDs to help sometimes
 *                  -- like having a INIT and a FINALIZE id that does nothing but to mark
 *                   where the rastro begins or ends for some given container
 *
 *        Version:  1.0
 *        Created:  04-10-2011 13:48:25 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "dummyevent.hh"



Paje::DummyEvent::DummyEvent(const string & name, Paje::BaseEventType * evt_type)
{
  BaseEvent();
  this->name = name;
  this->eventType = evt_type;
}

bool Paje::DummyEvent::do_start(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  *priority = 0;
  return true;
}


bool Paje::DummyEvent::do_end(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  *priority = 0;
  return true;
}


bool Paje::DummyEvent::do_trigger(double timestamp,
    symbols_table_t * symbols,
    double * priority,
    ostream &out)
{
  *priority = 0;
  return true;
}

bool Paje::DummyEvent::has_ids() const {
  return (( this->trigger_id | this->start_id | this->end_id) != 0);
}

void Paje::DummyEvent::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_TRIGGER,*base_id);
  (*base_id)++;
  
  this->set_trigger_id(EVENT_START,*base_id);
  (*base_id)++;

  this->set_trigger_id(EVENT_END, *base_id);
  (*base_id)++;
}

