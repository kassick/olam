// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/state.cc"
// Created: "Ter, 04 Out 2011 13:59:53 -0300 (kassick)"
// Updated: "Sex, 07 Out 2011 15:21:58 -0300 (kassick)"
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
#include "statetype.hh"
#include "semantic_types.hh"
#include <algorithm>
#include <typeinfo>

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



bool Paje::State::fits_in_event_type(
    const Paje::BaseEventType * evt_type) const
{
  const Paje::StateType * s = dynamic_cast<const Paje::StateType *>(evt_type);
  //cerr << "Check if is state type" << endl;
  return (s != NULL);
}







//-----------------------------------------------

struct pending_state {
  string container;
  Paje::BaseEvent* evt;
  double priority;
};

void check_opened_states(ostream & out)
{
  vector<struct pending_state> pending;
  for_each(event_names->begin(), event_names->end(), 
      [&](pair<string,Paje::BaseEvent*> p)
      {
        Paje::BaseEvent * evt = p.second;
        if (!evt->eventType)
          return;

        for (auto it = evt->timestamp_map.begin(); it != evt->timestamp_map.end(); ++it)
        {
          const string &container = it->first;
          stack<double> &ts_stack = it->second;

            while (ts_stack.size())
            {
              struct pending_state s;
              s.container = container;
              s.evt = evt;
              s.priority = ts_stack.top();
              pending.push_back(s);
              ts_stack.pop();
            }
          }
      });

  sort(pending.begin(),pending.end(),
      [](const struct pending_state e1, const struct pending_state e2)
      {
        return (e1.priority > e2.priority);
      }
    );

  for_each(pending.begin(),pending.end(),
      [&out](struct pending_state s) {
        cerr << "[Warning] Popping missed state " << s.evt->name << " on " << s.container << endl;
        pajePopState(s.priority,
                  s.container,
                  s.evt->eventType->typeName,
                  out);
      }
    );
}




int close_pending_states(double timestamp,string containerName, Paje::Container * c, ostream & out)
{
  

  vector<struct pending_state> pending;
  for_each(event_names->begin(), event_names->end(), 
      [&](pair<string,Paje::BaseEvent*> p)
      {
        Paje::BaseEvent * evt = p.second;

        // filter out events that are not of the type we are using
        if (( !evt->eventType) ||
            ( evt->eventType->container != c) )
          return;

        // no event pushed in this container
        if (!(evt->timestamp_map.count(containerName)))
          return;

        stack<double> &ts_stack = evt->timestamp_map[containerName];

        while (ts_stack.size())
        {
          struct pending_state s;
          s.container = containerName;
          s.evt = evt;
          s.priority = ts_stack.top();
          pending.push_back(s);
          ts_stack.pop();
        }
      });

  sort(pending.begin(),pending.end(),
      [](const struct pending_state e1, const struct pending_state e2)
      {
        return (e1.priority > e2.priority);
      }
    );

  for_each(pending.begin(),pending.end(),
      [&out,&timestamp](struct pending_state s) {
        cerr << "[Warning] Popping missed state " << s.evt->name << " on " << s.container << endl;
        pajePopState(timestamp,
                  s.container,
                  s.evt->eventType->typeName,
                  out);
      }
    );

  return pending.size();
}


