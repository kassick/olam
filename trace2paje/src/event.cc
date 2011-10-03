// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.cc"
// Created: "Sex, 02 Set 2011 15:23:14 -0300 (kassick)"
// Updated: "Seg, 03 Out 2011 18:45:10 -0300 (kassick)"
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
#include <stack>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "semantic_types.hh"

extern "C" {
#include <rastro.h>
}

using namespace std;
  
set<pair<string,string>> Paje::container_unique_names;


/*******************************************************************************
 Paje::EventType functions
*******************************************************************************/

void Paje::EventType::do_header(ostream &out)
{
  pajeDefineStateType(typeName, container->typeName, typeName, out);
}


Paje::EventType::EventType(const string & typeName) {
  this->typeName = typeName;
  this->container = NULL;
}

Paje::EventType::EventType(const string & typeName,Paje::Container * c) {
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
  timestamp_stack.push(DEFAULT_EVENT_PRIO);
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

  if (evt_id == start_id) {
    this->push_timestamp(timestamp);
    return do_start(timestamp,symbols,out);
  }

  if (evt_id == end_id)
  {
    bool ret = do_end(timestamp,symbols,out);
    this->pop_timestamp(timestamp);
    return ret;
  }
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
          rst_function_signature_buf << entry.type;
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
bool Paje::Event::load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols)
{
  int count_c, count_w, count_i, count_l, count_f, count_d, count_s;
  // this is a trick: state needs no symbols in the end; link needs the
  // ones needed by key. Solution: do not enforce loading ALL symbols if
  // it's an end symbol. Consistency is left to the user: all symbols in
  // the end of a link that are needed to form a key MUST BE in the
  // symbols table, but we are not gonna enforce it
  bool needs_all_symbols = true;

  Paje::Symbol *symbol;

  if (id == end_id)
    needs_all_symbols = false;
    //return true; // we don't have to load symbols for the end trigger

      count_c= count_w= count_i= count_l= count_f= count_d= count_s = 0 ;

#define CASE_TYPE(SHORT_NAME , PAJE_NAME) \
        case SHORT_NAME:                                                             \
          if (count_ ## SHORT_NAME < event->ct.n_ ## PAJE_NAME) {              \
            symbol.set_value( event->v_ ## PAJE_NAME [ count_ ## SHORT_NAME ] );\
            count_ ## SHORT_NAME ++  ;                                      \
          } else if (needs_all_symbols) {                                   \
            cerr << "Not enought values of type " << #SHORT_NAME  <<endl;   \
            return false;                                                   \
          }                                                                 \
          break; // trick?



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

#undef CASE_TYPE

  (*symbols)["EVT_NAME"].set_value(this->name.c_str());

}



/*************************************
 * Event Priorities:
 * the priorities define the order in which the triggers will be invoked
 * The natural priority order is as follows:
 *   StatePop == EventTrigger > 
 *    ContainerCreate >
 *      ContainerDestroy >
 *         StatePush
 *
 * That means that
 *    * Containers will be created before new states are pushed
 *    * Containers will be destroyed before new states are pushed (
 *      that means you can not have a state in a container if it's id
 *      is the one that triggers container create, but this makes sense
 *      cause the only state that would be valid here would be one with
 *      start_time == end_time)
 *    * States will be popped before new ones are pushed (avoiding
 *      imbrication conflicts)
 *    
 *  One important thing: the priority of StatePop is proportional to the
 *  last time that state was pushed. That means that imbrication will
 *  be mantained consistent
 *  */

void Paje::Event::push_timestamp(const double timestamp)
{
  this->timestamp_stack.push(timestamp);
}

double Paje::Event::pop_timestamp(const double timestamp)
{
  // ignore parameter, we do not need in the default case
  double ret = this->timestamp_stack.top();
  this->timestamp_stack.pop();
  return ret;
}


double Paje::Event::get_priority() const
{
  return this->timestamp_stack.top();
}


bool Paje::Event::operator<(const Paje::Event* other) const
{
  return (this->get_priority() > other->get_priority());
}


bool Paje::Event::has_ids() const {
  return (this->trigger_id != 0);
}

void Paje::Event::gen_auto_ids(long int * base_id)
{
  this->set_trigger_id(EVENT_TRIGGER,*base_id);
  (*base_id)++;
}

void Paje::Event::gen_ids_description(ostream & out)
{
  if (this->trigger_id)
    out << "ID " << this->name << " " << this->trigger_id << endl;

  if (this->start_id)
    out << "ID " << this->name << " START " << this->start_id << endl;
  
  if (this->end_id)
    out << "ID " << this->name << " END " << this->end_id << endl;
}


void Paje::Event::gen_c_header(ostream & out, 
    const string & start_suffix,
    const string & end_suffix,
    const string & evt_suffix  )
{

  string caps_name = boost::to_upper_copy(this->name);

  if (this->trigger_id)
    out << "#define " << caps_name << evt_suffix << 
            " (" << this->trigger_id << ")" << endl;

  if (this->start_id)
    out << "#define " << caps_name << start_suffix << 
            " (" << this->start_id << ")" << endl;

  if (this->end_id)
    out << "#define " << caps_name << end_suffix << 
            " (" << this->end_id << ")" << endl;
}

void Paje::Event::gen_fort_header(ostream & out,
          const string & start_suffix,
          const string & end_suffix,
          const string & evt_suffix)
{
  string caps_name = boost::to_upper_copy(this->name);
  if (this->trigger_id)
    out << "integer :: " <<  caps_name << evt_suffix << 
            " = " << this->trigger_id << endl;

  if (this->start_id)
    out << "integer :: " << caps_name << start_suffix << 
            " = " << this->start_id << endl;

  if (this->end_id)
    out << "integer :: " << caps_name << end_suffix << 
            " = " << this->end_id << endl;

}

void Paje::Event::get_rst_function_signature(set<string> & signatures)
{
  signatures.insert(this->rst_function_signature_buf.str());
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
 * Paje Dummy Event
 * This one is created for referenced events with id but that have no EVENT
 * or STATE defined. They act as nicknames for some IDs to help sometimes
 * -- like having a INIT and a FINALIZE id that does nothing but to mark
 *  where the rastro begins or ends for some given container
 ******************************************************************************/


Paje::DummyEvent::DummyEvent(const string & name, Paje::EventType * evt_type)
{
  Event();
  this->name = name;
  this->eventType = evt_type;
}

bool Paje::DummyEvent::do_start(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  return true;
}


bool Paje::DummyEvent::do_end(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  return true;
}


bool Paje::DummyEvent::do_trigger(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
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





//******************************************************************************
// Class ContainerTrigger
// called to execute PajeCreateContainer/PajeDestroyContainer when the
// correct id is seen in the rastro


Paje::ContainerCreateTrigger::ContainerCreateTrigger(Paje::Container * c,hierarchy_t * n)
{
  this->container = c;
  this->hierarchy = n;
  push_timestamp(CONTAINER_CREATE_PRIO); // can leave the EVENT_DEFAULT_PRIO on the stack, won't hurt
}


bool Paje::ContainerCreateTrigger::do_start(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName, parentName;
  
  if (this->container->parent)
    containerName =  format_values(container->parent->formatName, symbols);
  else
    containerName = PAJE_ROOT_CONTAINER;


  walk_tree_head_first(hierarchy,[&](hierarchy_t * n, int level)
      {
        Paje::Container * c = n->getVal();
        parentName = containerName;

        //cerr << "Creating " << c->typeName << endl;

        if (c->triggerParent || (c == this->container) ) {
          // do this for the current container and all it's create on parent children
          containerName = format_values(c->formatName, symbols);

          // do we already have a container with this name, of this type?
          pair<string,string> p(containerName,c->typeName);
          if (Paje::container_unique_names.count(p))
          {
              return true; // ignore
              // can happen, if a container is triggered by event1, and
              // event1 happens several times
          }

          Paje::container_unique_names.insert(p);

          pajeCreateContainer(timestamp,
                                 containerName,
                                 c->typeName,
                                 parentName,
                                 containerName,
                                 out);
          return false;
        }
        return true; // does not create on parent, this one and it's children shall not be created here, there will be an event for it

      });

  return true;
}


bool Paje::ContainerCreateTrigger::do_end(double timestamp,
    symbols_table_t * symbols, ostream &out)
{
  string containerName;
  
  // needs to destroy all children
  walk_tree_depth_first(this->hierarchy,[&](hierarchy_t * n, int level) 
      {
        Paje::Container * c = n->getVal();
        containerName = format_values(c->formatName, symbols);
        pajeDestroyContainer(timestamp,
                              c->typeName,
                              containerName,
                              out);
        return false;
      }
    );

  return true;
}


void Paje::ContainerCreateTrigger::push_timestamp(double timestamp)
{
  timestamp_stack.push(CONTAINER_DESTROY_PRIO);
}


