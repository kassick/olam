// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/baseevent.cc"
// Created: "Ter, 04 Out 2011 11:51:35 -0300 (kassick)"
// Updated: "Qui, 17 Nov 2011 18:50:45 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  baseevent.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 11:51:35 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "baseevent.hh"
#include "baseeventtype.hh"
#include "semantic_types.hh"

#include <cassert>
#include <boost/algorithm/string.hpp>
#include <algorithm>



using namespace std;




Paje::BaseEvent::BaseEvent() {
  this->formatValue = "%(EVT_NAME)";
  start_id = end_id = trigger_id = 0;
  this->eventType = NULL;
}

void Paje::BaseEvent::set_event_type(Paje::BaseEventType * evt_type)
{
  this->eventType = evt_type;

}

double Paje::BaseEvent::get_priority(event_id_t evt_id, double timestamp, symbols_table_t ** symbols)
{
  if (this->eventType) {
    if (evt_id == end_id) {
      string containerName = format_values(this->eventType->container->formatName,symbols);
      return get_timestamp(containerName,timestamp);
    }

    if (evt_id == start_id)
    {
      return DEFAULT_EVENT_PRIO;
    }
  }
 
  // defaults to return a trigger prio
  return DEFAULT_EVENT_PRIO_TRIGG;
}

/////
// Trigger functions
bool Paje::BaseEvent::do_start(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) {
  if (this->eventType)
  {
    string containerName = format_values(this->eventType->container->formatName,symbols);
    this->push_timestamp(containerName,timestamp);
  }
  
  *priority = DEFAULT_EVENT_PRIO;
  
  return false;
}



bool Paje::BaseEvent::do_end(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) {
  
  if (this->eventType)
  {
    string containerName = format_values(this->eventType->container->formatName,symbols);
  
    *priority = pop_timestamp(containerName,timestamp);
  } else {
    *priority = DEFAULT_EVENT_PRIO;
  }
  
  return false;
}

bool Paje::BaseEvent::do_trigger(double timestamp,
          symbols_table_t ** symbols,
          double * priority,
          ostream &out) {
  *priority = DEFAULT_EVENT_PRIO_TRIGG;

  return false;
}




/////
//Sets the trigger field
void Paje::BaseEvent::set_trigger_id(trigger_id_t trigger_field, event_id_t id)
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

bool Paje::BaseEvent::trigger(event_id_t evt_id, double timestamp,
    symbols_table_t ** symbols, 
    double * priority,
    ostream &out)
{
  if (evt_id == trigger_id)
    return do_trigger(timestamp,symbols,priority,out);

  if (evt_id == start_id) {
    return do_start(timestamp,symbols,priority,out);
  }

  if (evt_id == end_id)
  {
    bool ret = do_end(timestamp,symbols,priority,out);
    return ret;
  }
}


void Paje::BaseEvent::add_symbol_from_tree(attribs_t * attrs)
{
  char type;
  identifier_entry_t entry;

  assert(attrs->getVal()->id == ID_RASTRO_TYPE);

  // sadly, #SHORT_NAME[0] does not work as the compiler does not recognize
  // str_const[int_const] as a char in compilation time
#define CASE_TYPE1(CHAR_NAME,SHORT_NAME,VALUE) \
  case CHAR_NAME : \
    VALUE = SHORT_NAME;   \
    break;                \



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
          return (level > 0);
          break; // ignore; code done above
        case ID_IDF:
          entry.field_name = attr->vals.name;
          this->identifier_names.push_back(entry);
          rst_function_signature_buf << entry.type;
          break;
        default:
          cerr << "Unexpected value in tree: " << attr->toString()
              <<  " while creating " << this->name << endl;
          exit(1);
      }
      return false;
      });

}

///
void Paje::BaseEvent::fill_from_attr(attribs_t * attrs)
{
  pushlist_entry_t e;
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
          case ID_PUSH_START_PARAM:
            assert(n->childrenCount() == 2);
            for(auto an_it = n->begin(); an_it != n->end(); ++an_it)
            {
              if ((*an_it)->getVal()->id == ID_IDF) {
                e.first = (*an_it)->getVal()->vals.identifier_name;
              } else if ((*an_it)->getVal()->id == ID_AS_IDF) {
                e.second = (*an_it)->getVal()->vals.identifier_name;
              }
            }
            pushlist_start.insert(e);
            break;
          case ID_PUSH_END_PARAM:
            assert(n->childrenCount() == 2);
            for(auto an_it = n->begin(); an_it != n->end(); ++an_it)
            {
              if ((*an_it)->getVal()->id == ID_IDF) {
                e.first = (*an_it)->getVal()->vals.identifier_name;
              } else if ((*an_it)->getVal()->id == ID_AS_IDF) {
                e.second = (*an_it)->getVal()->vals.identifier_name;
              }
            }
            pushlist_end.insert(e);
            break;

          default:
            break;
        }
        return false;
      });


    fill_entry_list(ID_MAP    , map_list, attrs);
    fill_entry_list(ID_MAP_GET, get_list, attrs);
}


////
bool Paje::BaseEvent::load_symbols(event_id_t id, rst_event_t *event, symbols_table_t * symbols)
{
  int count_c, count_w, count_i, count_l, count_f, count_d, count_s;
  // this is a trick: state needs no symbols in the end; link needs the
  // ones needed by key. Solution: do not enforce loading ALL symbols if
  // it's an end symbol. Consistency is left to the user: all symbols in
  // the end of a link that are needed to form a key MUST BE in the
  // symbols table, but we are not gonna enforce it
  bool needs_all_symbols = true;
  bool ret = true;

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
            cerr << "[Warning] Not enought values of type " << #SHORT_NAME  \
                 << " for id " << id                                        \
                 << " in event " << this->name                              \
                 << endl;                                                   \
            ret = false;                                                    \
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
        ret = false;
        break;
    }

  }

#undef CASE_TYPE

  (*symbols)["EVT_NAME"].set_value(this->name.c_str());

  return ret;

}

void Paje::BaseEvent::push_symbols(Paje::event_id_t id,
    symbols_table_t * from,
                                   symbols_table_t * to)
{
  if ((id == this->trigger_id) || (id == this->start_id))
  {
      // now save the ids that must be saved in the local table
      for (auto push_it = this->pushlist_start.begin();
                push_it != this->pushlist_start.end();
                ++ push_it)
      {
        //cerr << "Push " << push_it->first << " as " << push_it->second << endl;
        (*to )[push_it->second] = (* from)[push_it->first];
      }
  } else if (id == this->end_id) {
    // now save the ids that must be saved in the local table
      for (auto push_it = this->pushlist_end.begin();
                push_it != this->pushlist_end.end();
                ++ push_it)
      {
        //cerr << "Push " << push_it->first << " as " << push_it->second << endl;
        (*to )[push_it->second] = (* from)[push_it->first];
      }
  }
}


/*************************************
 * Push Symbols to Map
 * */

void Paje::BaseEvent::map_symbols(Paje::event_id_t evt_id, 
    symbols_table_t ** symbols,
    user_defined_maps_t & usermaps)
{
  for (auto it = map_list.begin(); it != map_list.end(); ++it)
  {
    int tmp_symbol = 0;
    // either it happens on BOTH or it needs to be defined to happen on
    // START or END
    if ((it->when & USER_MAP_BOTH) || 
        ( (it->when & USER_MAP_START) && (evt_id == start_id) ) ||
        ( (it->when & USER_MAP_END  ) && (evt_id == end_id  ) )
       )
    {
      string mapname = format_values(it->map_name_format, symbols);

      // look for a symbol named idf_name on all tables
      bool mapped = false;
      while (symbols[tmp_symbol])
      {
        if (symbols[tmp_symbol]->count(it->idf_name))
        {
          string key = format_values(it->key, symbols);
          // copy the symbol from the symbol table
          (usermaps[mapname])[key] = (*(symbols[tmp_symbol]))[it->idf_name] ;
          //cerr << "Mapping " << key << " to ";
          //       (*(symbols[tmp_symbol]))[it->idf_name].format("",cerr);
          //cerr << "(" << it->idf_name << ") in map " <<  mapname << endl;


          mapped = true;
          break;
        }
        tmp_symbol++;
      }

      if (!mapped) {
        cerr << "[Warning] Could not map " << mapname << "[" << it->key << "] = " << it->idf_name << endl;
      }
    }
  }
}



void Paje::BaseEvent::get_symbols_from_map (Paje::event_id_t evt_id, 
    symbols_table_t * dest_symbols,
    symbols_table_t ** symbols,
    user_defined_maps_t & usermaps)
{
  

  for (auto it = get_list.begin(); it != get_list.end(); ++it)
  {
    // either it happens on BOTH or it needs to be defined to happen on
    // START or END

    if ((it->when & USER_MAP_BOTH) || 
        ( (it->when & USER_MAP_START) && (evt_id == start_id) ) ||
        ( (it->when & USER_MAP_END  ) && (evt_id == end_id  ) )
       )
    {
      string mapname = format_values(it->map_name_format, symbols);

      if (!(usermaps.count(mapname)))
      {
        cerr << "[Warning] No user map named " << mapname << endl;
        continue;
      }

      string key = format_values(it->key, symbols);
      if (! ((usermaps[mapname]).count(key))) {
        cerr << "[Warning] Map " << mapname << " has no key " << key << endl;
        continue;
      }

      (*dest_symbols)[it->idf_name] = ((usermaps[mapname])[key]);

    }
  }
}



/*************************************
 * BaseEvent Priorities:
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

void Paje::BaseEvent::push_timestamp(const string & containerName, const double timestamp)
{
  this->timestamp_map[containerName].push(timestamp);
}

double Paje::BaseEvent::pop_timestamp(const string & containerName, const double timestamp)
{
  if (this->timestamp_map.count(containerName))
  {
    double ret = this->timestamp_map[containerName].top();
    this->timestamp_map[containerName].pop();
    return ret;
  } else {
    return DEFAULT_EVENT_PRIO;
  }
}

double Paje::BaseEvent::get_timestamp(const string & containerName, const double timestamp)
{
  if (this->timestamp_map.count(containerName))
  {
    double ret = this->timestamp_map[containerName].top();
    return ret;
  } else {
    return DEFAULT_EVENT_PRIO;
  }
}

#if 0
double Paje::BaseEvent::get_priority(const string & containerName) const
{
  if (this->timestamp_map.count(containerName)) 
  {
    stack<double> &a = (this->timestamp_map)[containerName];
    return a.top();
  }
  else
    return DEFAULT_EVENT_PRIO;
}
#endif


#if 0
bool Paje::BaseEvent::operator<(const Paje::BaseEvent* other) const
{
  return (this->get_priority() > other->get_priority());
}
#endif


bool Paje::BaseEvent::has_ids() const {
  return (this->trigger_id != 0);
}

void Paje::BaseEvent::gen_auto_ids(long int * base_id, set<event_id_t> & unique_ids)
{
  while (unique_ids.count(*base_id))
  {
    (*base_id)++;
  }
  //this->set_trigger_id(EVENT_TRIGGER,*base_id);
}

void Paje::BaseEvent::gen_ids_description(ostream & out)
{
  if (this->trigger_id)
    out << "ID " << this->name << " " << this->trigger_id << endl;

  if (this->start_id)
    out << "ID " << this->name << " START " << this->start_id << endl;
  
  if (this->end_id)
    out << "ID " << this->name << " END " << this->end_id << endl;
}


void Paje::BaseEvent::gen_c_header(ostream & out, 
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

void Paje::BaseEvent::gen_fort_header(ostream & out,
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

void Paje::BaseEvent::get_rst_function_signature(set<string> & signatures)
{
  signatures.insert(this->rst_function_signature_buf.str());
}


string Paje::BaseEvent::toString() {
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
  for (auto it = identifier_names.begin(); it != identifier_names.end(); ++it)
    out << "      " << (*it).type << " " << (*it).field_name << endl;

  out << "   " << "Pushes:" << endl;
  for (auto it = pushlist_start.begin(); it != pushlist_start.end(); ++ it)
  {
    out << "      " << (*it).first << " => " << (*it).second << "(start)" << endl;
  }
  for (auto it = pushlist_end.begin(); it != pushlist_end.end(); ++ it)
  {
    out << "      " << (*it).first << " => " << (*it).second << "(end)" << endl;
  }

  out << "   " << "Maps:" << endl;
  for (auto it =  map_list.begin(); it != map_list.end(); ++it)
  {
    out << "      " << it->map_name_format 
        << "[" << it->key << "] = " << it->idf_name << "( ";
    if (it->when & USER_MAP_START)
      out << "start ";
    if (it->when & USER_MAP_END)
      out << "end ";
    out << ")" << endl;
  }

  out << "   " << "Gets:" << endl;
  for (auto it =  get_list.begin(); it != get_list.end(); ++it)
  {
    out << "      " << it->idf_name << " = "
      << it->map_name_format  << "[" << it->key << "]" << "( " ;

    if (it->when & USER_MAP_START)
      out << "start ";
    if (it->when & USER_MAP_END)
      out << "end ";

    out << ")" << endl;
  }

  return out.str();
}



bool Paje::BaseEvent::fits_in_event_type(
    const Paje::BaseEventType * evt_type) const
{
  // non-specialized event type should fit anywhere
  return true;
}
