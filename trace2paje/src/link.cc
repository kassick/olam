// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/link.cc"
// Created: "Ter, 04 Out 2011 14:03:18 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:06:01 -0300 (kassick)"
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


Paje::Link::Link(string &name, attribs_t * attribs) {
  this->name = name;

  this->fill_from_attr(attribs);
}


string Paje::Link::toString() {
  stringstream out;

  out << Paje::BaseEvent::toString();
  out << "   " << "Key format: " << format_key << endl;

  return out.str();
}




/////
// Trigger functions
bool Paje::Link::do_start(double timestamp,
          symbols_table_t * symbols, ostream &out) {

  LinkType * lt = (LinkType * ) this->eventType;

  string thisContainer, sourceContainer, thisValue, key;

  thisContainer = format_values(lt->container->formatName, symbols);
  sourceContainer = format_values(lt->source->formatName, symbols);
  key  = format_values(this->format_key, symbols);

  pajeStartLink(timestamp,
                   thisContainer,
                   lt->typeName,
                   sourceContainer,
                   thisValue,
                   key,
                   out);

  return false;
}

bool Paje::Link::do_end(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  LinkType * lt = (LinkType * ) this->eventType;

  string thisContainer, destContainer, thisValue, key;

  thisContainer = format_values(lt->container->formatName, symbols);
  destContainer = format_values(lt->source->formatName, symbols);
  key  = format_values(this->format_key, symbols);
  pajeEndLink(timestamp,
                 thisContainer,
                 lt->typeName,
                 destContainer,
                 thisValue,
                 key,
                 out);
  return false;
}

bool Paje::Link::do_trigger(double timestamp,
          symbols_table_t * symbols, ostream &out) {
  // Actually, Link should call PajeLink....
  cerr << "Error: Class Link has no trigger action" << endl;
  return false;
}

void Paje::Link::fill_from_attr(attribs_t * attrs) {
  Paje::BaseEvent::fill_from_attr(attrs);
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


bool Paje::Link::has_ids() const {
  return ((this->start_id != 0) && (this->end_id != 0));
}

void Paje::Link::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_START,*base_id);
  (*base_id)++;

  this->set_trigger_id(EVENT_END, *base_id);
  (*base_id)++;
}
