// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Ter, 27 Set 2011 17:09:45 -0300 (kassick)"
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
#include <assert.h>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include "semantics.hh"

extern "C" {
#include <rastro.h>
}

using namespace std;


/*******************************************************************************
 Paje::EventType functions
*******************************************************************************/

void Paje::EventType::do_header(ostream &out)
{
  pajeDefineStateType(typeName, container->typeName, typeName, out);
}


Paje::EventType::EventType(string & typeName) {
  this->typeName = typeName;
}

Paje::EventType::EventType(string & typeName,Paje::Container * c) {
  this->typeName = typeName;
  this->container = c;

}

#if 0
// this is unused for now
Paje::EventType::EventType(string & typeName, attribs_t * attribs)
{
  this->typeName = typeName;

  // walks the parent node (a Container, we hope...)
  walk_tree_up_to_root(attribs,[&](attrib_t * n, int level) {
      SemanticAttribute * attr = n->getVal();
      if (attr->id == ID_CONTAINER) {
        if (! container_type_names->count(attr->vals.name))
        {
          cerr << "Error: EventType " << typeName << "can not find it's parent container " << attr->vals.name << endl;
          exit(1);
        }

        this->container = (*container_type_names)[attr->vals.name];
        return true; // stop now, otherwise it will get to the first containers!
      }
    });
}
#endif

/*******************************************************************************
 Paje::LinkType functions
*******************************************************************************/

Paje::LinkType::LinkType(string& tn, Paje::Container * c, attribs_t * t) : 
  Paje::EventType(tn,c)
{
  this->source = NULL;
  this->dest  = NULL;
  
  attribs_t::iterator it;
  for (it = t->begin(); it != t->end(); ++it)
  {
    SemanticAttribute * attr = (*it)->getVal();
    switch(attr->id) {
      case ID_LINK_SOURCE:
        if (!container_type_names->count(attr->vals.name)) {
          cerr << "Error: Can not find source container for link type " << typeName << endl;
          exit(1);
        }
        this->source = (*container_type_names)[attr->vals.name]->getVal();
        break;
      case ID_LINK_DEST:
        if (!container_type_names->count(attr->vals.name)) {
          cerr << "Error: Can not find source container for link type " << typeName << endl;
          exit(1);
        }
        this->dest = (*container_type_names)[attr->vals.name]->getVal();
        break;
      default:
        break;
    }
  }

  if ((this->source == NULL) || (this->dest == NULL) || (this->container == NULL)) {
    cerr << "Null fields, verify!" << endl;
    exit(1);
  }
}

void Paje::LinkType::do_header(ostream &out)
{
  pajeDefineLinkType(typeName, container->typeName,
                        source->typeName,
                        dest->typeName,
                        typeName,
                        out);
}

/*******************************************************************************
 Paje::Event functions

*******************************************************************************/

Paje::Event::Event() {
  start_id = end_id = trigger_id = 0;
  this->eventType = NULL;
}

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


void Paje::Event::add_symbol_from_tree(attribs_t * attrs)
{
  char type;
  identifier_entry_t entry;

  assert(attrs->getVal()->id == ID_RASTRO_TYPE);

  // gets the type of the child variables
  type = attrs->getVal()->vals.name[0];
  switch (type) {
    CASE_TYPE1('c',c,entry.type);
    CASE_TYPE1('w',w,entry.type);
    CASE_TYPE1('i',i,entry.type);
    CASE_TYPE1('l',l,entry.type);
    CASE_TYPE1('f',f,entry.type);
    CASE_TYPE1('d',d,entry.type);
    CASE_TYPE1('s',s,entry.type);
    default:
      cerr << "Unknown type " << type << endl;
      exit(1);
  }


  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
      SemanticAttribute * attr = n->getVal();

      switch (attr->id) {
        case ID_RASTRO_TYPE:
          break; // ignore; code done above
        case ID_IDF:
          entry.field_name = attr->vals.name;
          this->identifier_names.push_back(entry);
          break;
        default:
          cerr << "Unexpected value in tree: " << attr << endl;
          exit(1);
      }
      return false;
      });

}

///
void Paje::Event::fill_from_attr(attribs_t * attrs)
{
  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
        SemanticAttribute * attr = n->getVal();
        switch (attr->id) {
          case ID_FORMAT_NAME:
            this->formatValue = attr->vals.name;
            break;
          case ID_EVENT_TYPE:
            if (eventtype_names->count(attr->vals.name)) {
              this->eventType = (*eventtype_names)[attr->vals.name];
            } else {
              cerr << "Unknown type: " << attr->vals.name << endl;
              exit(1);
            }
            break;
          case ID_RASTRO_TYPE:
            this->add_symbol_from_tree(n);
            return false; // no need to descend into this branch
            break;
          default:
            break;
        }
        return false;
      });
}


////
bool Paje::Event::load_symbols(rst_event_t *event, symbols_table_t * symbols)
{
  int count_c, count_w, count_i, count_l, count_f, count_d, count_s;
  Paje::Symbol *symbol;

      count_c= count_w= count_i= count_l= count_f= count_d= count_s = 0 ;


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

string Paje::Event::toString() {
  stringstream out;

  out << "Paje::" << this->type_identifier << "``" << this->name << "´´" << endl;
  out << "   " << "Formated to " << this->formatValue << endl;
  out << "   " << "Type: " << (this->eventType == NULL? "NULL": this->eventType->typeName) << endl;
  out << "   " << "With IDs:";
  if (start_id)   out << " Start: " << start_id;
  if (end_id)     out << " End: " << end_id;
  if (trigger_id) out << " Trigger: " << trigger_id;
  out << endl;

  out << "   " << "Fields:" << endl;
  identifier_list_t::iterator it;
  for (it = identifier_names.begin(); it != identifier_names.end(); ++it)
    out << "      " << (*it).type << " " << (*it).field_name << endl;

  return out.str();
}




/*******************************************************************************
 * Paje State functions
 ******************************************************************************/

Paje::State::State(string& name, attribs_t * attribs) {
  this->name = name;
  this->type_identifier = "State";

  this->fill_from_attr(attribs);
}

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
  Paje::Event::fill_from_attr(attrs);
  return;
}



/*******************************************************************************
 * Paje::Link class
 ******************************************************************************/

Paje::Link::Link(string &name, attribs_t * attribs) {
  this->name = name;

  this->fill_from_attr(attribs);
}


string Paje::Link::toString() {
  stringstream out;

  out << Paje::Event::toString();
  out << "   " << "Key format: " << format_key << endl;

  return out.str();
}




/////
// Trigger functions
bool Paje::Link::do_start(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  cerr << "Error: Class Link has no start action" << endl;
  return false;
}

bool Paje::Link::do_end(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  cerr << "Error: Class Link has no end  action" << endl;
  return false;
}

bool Paje::Link::do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  // Actually, Link should call PajeLink....
  cerr << "Error: Class Link has no trigger action" << endl;
  return false;
}

void Paje::Link::fill_from_attr(attribs_t * attrs) {
  Paje::Event::fill_from_attr(attrs);
  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
        SemanticAttribute * attr = n->getVal();
        switch (attr->id) {
          case ID_LINK:
            break; // name has already been set
          case ID_KEY_FORMAT:
            this->format_key = attr->vals.name;
            break;
          default:
            break;
        }


        return false;
      });
}

