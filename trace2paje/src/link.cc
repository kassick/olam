// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/link.cc"
// Created: "Ter, 04 Out 2011 14:03:18 -0300 (kassick)"
// Updated: "Qua, 09 Nov 2011 18:41:55 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  link.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 14:03:18 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */


#include "link.hh"
#include "paje_functions.hh"


#include <map>
#include <list>
#include <algorithm>
#include <typeinfo>

typedef pair <string,double> ts_list_entry_t;
typedef list<ts_list_entry_t> timestamp_list_t;
typedef struct _link_key_t {
  string key,
         typeName,
         containerName;

  bool operator<(const struct _link_key_t & other) const
  {
    return (this->key < other.key);
  }

  bool operator==(const struct _link_key_t & other) const
  {
    return((this->key           == other.key) &&
           (this->typeName      == other.typeName) &&
           (this->containerName == other.containerName) );
  }
} link_key_t;

typedef map<link_key_t, timestamp_list_t  > key_to_timestamp_map_t;

key_to_timestamp_map_t pending_links, started_links;


void check_missing_links() {

  for_each(pending_links.begin(), pending_links.end(),
      [](pair<link_key_t, timestamp_list_t> p)
      {
        if (p.second.size() > 0) {
          cerr << "[Warning] Link with type " << p.first.typeName
               << ", key " << p.first.key 
               << " on " << p.first.containerName
               << " has " << p.second.size() << " lost arrows (pending)" << endl;
        }
      } );
  
  for_each(started_links.begin(), started_links.end(),
      [](pair<link_key_t, timestamp_list_t> p)
      {
        if (p.second.size() > 0) {
          cerr << "[Warning] Link with type " << p.first.typeName
               << ", key " << p.first.key 
               << " on " << p.first.containerName
               << " has " << p.second.size() << " lost arrows (started)" << endl;
        }
      } );
}



Paje::Link::Link(string &name, attribs_t * attribs) {
  this->name = name;

  this->fill_from_attr(attribs);
}


string Paje::Link::toString() {
  stringstream out;

  out << Paje::BaseEvent::toString();
  out << "   " << "Key format: " << format_key_start << "|" << format_key_end << endl;

  return out.str();
}




/////
// Trigger functions
bool Paje::Link::do_start(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out)
{
  link_key_t link_key;
  LinkType * lt = (LinkType * ) this->eventType;

  //Paje::BaseEvent::do_start(timestamp, symbols, priority, out);
  *priority = 0;

  string sourceContainer, thisValue;

  link_key.containerName = format_values(lt->container->formatName, symbols);
  link_key.key           = format_values(this->format_key_start, symbols);
  link_key.typeName      = lt->typeName;
  sourceContainer = format_values(lt->source->formatName, symbols);
  thisValue       = format_values(this->formatValue_start,symbols);

  pajeStartLink(timestamp,
                 link_key.containerName,
                 lt->typeName,
                 sourceContainer,
                 thisValue,
                 link_key.key,
                 out);

  // Special case: if start == end, we must call the end here and exit in
  // do_end
  if (this->start_id == this->end_id) {
    string destContainer   = format_values(lt->dest->formatName, symbols);
    string otherValue = format_values(this->formatValue_end,symbols);
    pajeEndLink(timestamp,
                   link_key.containerName,
                   lt->typeName,
                   destContainer,
                   otherValue,
                   link_key.key,
                   out);
    return true;
  }
  
 
  bool add_to_started = true;
  if (pending_links.count(link_key)) {
    timestamp_list_t &ts_list = pending_links[link_key];

    // check if we have arrows waiting in the queue to be finalized
    // normal case: there is no conflicting key on the trace
    if (ts_list.size() > 0)
    {
      // Some arrow ended before it began
      ts_list_entry_t entry = ts_list.front();
      ts_list.pop_front();
      add_to_started = false;

      if (entry.second < timestamp) {
        cerr << "[Warning] Arrow of type " << lt->typeName << " key " << link_key.key << " ends before it begins" << endl;
        // this should always happen in this block, the if is just so we make
        // it clear

      }

      string destContainer = format_values(lt->dest->formatName,symbols);

      pajeEndLink(timestamp,
                     link_key.containerName,
                     lt->typeName,
                     entry.first,
                     thisValue,
                     link_key.key,
                     out);
    }
  }

  if (add_to_started)
  {
    ts_list_entry_t entry;
    entry.first = sourceContainer;
    entry.second = timestamp;
    started_links[link_key].push_back(entry);
  }


  return true;
}

bool Paje::Link::do_end(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) 
{
  link_key_t link_key;
  LinkType * lt = (LinkType * ) this->eventType;
  
  //Paje::BaseEvent::do_end(timestamp, symbols, priority, out);
  *priority = 0;

  string destContainer, thisValue;
  
  // do not process the end event if ids are equal
  if (this->start_id == this->end_id) {
    return true;
  }

  link_key.containerName = format_values(lt->container->formatName, symbols);
  link_key.key           = format_values(this->format_key_end, symbols);
  link_key.typeName      = lt->typeName;
  destContainer   = format_values(lt->dest->formatName, symbols);
  thisValue       = format_values(this->formatValue_end,symbols);


  bool add_to_pending = true;
  if (started_links.count(link_key)) {
    timestamp_list_t &ts_list = started_links[link_key];
    
    // Normal case: we have one arrow with that key
    if (ts_list.size() > 0)
    {
      // if there are several arrows with the same key, take the oldest one
      ts_list_entry_t entry = ts_list.front();
      ts_list.pop_front(); // discard oldest reference
      add_to_pending = false;


      if (entry.second > timestamp)
      {
        cerr << "[Warning] Arrow of type " << lt->typeName << " key " << link_key.key << " ends before it begins" << endl;
        timestamp = entry.second; // flaten arrow
      }

      pajeEndLink(timestamp,
                     link_key.containerName,
                     lt->typeName,
                     destContainer,
                     thisValue,
                     link_key.key,
                     out);
    }
  }

  if (add_to_pending)
  {
    // we are finishing an arrow but there is none started that we can finish -- add it to
    // the list
    ts_list_entry_t entry;
    entry.first = destContainer;
    entry.second = timestamp;
    cerr << "adding to pending " << link_key.key << " " << link_key.typeName << " " << link_key.containerName << endl;
    started_links[link_key].push_back(entry);
  }
  
  return true;
}

bool Paje::Link::do_trigger(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) {
  // Actually, Link should call PajeLink....
  cerr << "Error: Class Link has no trigger action" << endl;
  *priority = 0;
  return false;
}

void Paje::Link::fill_from_attr(attribs_t * attrs) {

  Paje::BaseEvent::fill_from_attr(attrs);
  this->formatValue_start = this->formatValue;
  this->formatValue_end   = this->formatValue;
  
  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
        SemanticAttribute * attr = n->getVal();


        switch (attr->id) {
          case ID_LINK:
            break; // name has already been set
          case ID_KEY_FORMAT:
            this->format_key_start = attr->vals.name;
            this->format_key_end   = attr->vals.name;
            break;
          case ID_KEY_FORMAT_START:
            this->format_key_start = attr->vals.name;
            break;
          case ID_KEY_FORMAT_END:
            this->format_key_end = attr->vals.name;
            break;
          case ID_VALUE_FORMAT_START:
            this->formatValue_start = attr->vals.name;
            //cerr << "format val start is " << this->formatValue_start << endl;
            break;
          case ID_VALUE_FORMAT_END:
            this->formatValue_end   = attr->vals.name;
            break;
          default:
            break;
        }


        return false;
      });
}


bool Paje::Link::has_ids() const {
  return ((this->start_id != 0) && (this->end_id != 0));
}

void Paje::Link::gen_auto_ids(long int * base_id, set<event_id_t> & unique_ids)
{
  if (!this->start_id)
  {
    while (unique_ids.count(*base_id))
    {
      (*base_id)++;
    }
    unique_ids.insert(*base_id);
    this->set_trigger_id(EVENT_START,*base_id);
  }

  
  
  if (!this->end_id)
  {
    while (unique_ids.count(*base_id))
    {
      (*base_id)++;
    }
    unique_ids.insert(*base_id);
    this->set_trigger_id(EVENT_END, *base_id);
  }
}



bool Paje::Link::fits_in_event_type(
    const Paje::BaseEventType * evt_type) const
{
  const Paje::LinkType * s = dynamic_cast<const Paje::LinkType *>(evt_type);
  return (s != NULL);
}
