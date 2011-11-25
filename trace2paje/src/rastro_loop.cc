// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_loop.cc"
// Created: "Ter, 27 Set 2011 10:23:09 -0300 (kassick)"
// Updated: "Qui, 24 Nov 2011 21:11:22 -0600 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  rastro_loop.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27-09-2011 10:23:09 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "rastro_loop.hh"
#include "semantics.hh"
#include "paje.hh"
#include "rastro_helper.hh"
#include "logging.hh"
#include <ostream>
#include <string>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>


// librastro is a mess
extern "C" {
#include <rastro.h>
}


typedef struct _serve_entry_t {
  double priority;

  //string output;

  Paje::BaseEvent * evt;

  struct _serve_entry_t & operator=(const struct _serve_entry_t &org) 
  {
    if (this == &org)
      return *this;

    this->priority = org.priority;
    this->evt = org.evt;
    //this->output = org.output;

    return *this;
  }

} serve_entry_t;

static double  _rastro_loop_events(list<string> &files_to_open, ostream &out, user_defined_maps_t & usermaps,
    bool debug, bool do_map, bool do_output)
{

  // Open files, use a tmp sync file
  char syncfile[256];
  char * fname;
  set<string> opened_files;
  double timestamp;

#define _GLOBAL_TABLE 0
#define _LOCAL_TABLE 1
#define _END_TABLE 2
  symbols_table_t *symbols[4];

  rst_event_t event;
  rst_file_t data;
  set<Paje::event_id_t> ignored_ids;

  // What on earth!?
  data.initialized = 0;


  map<pair<decltype(rst_event_t::id1),decltype(rst_event_t::id2)>, symbols_table_t *>    symbols_per_file;

  symbols[_END_TABLE] = NULL;
  // symbols[1] is local-to-file symbols,  
  symbols[_GLOBAL_TABLE] = new symbols_table_t(); // global symbols / dirty table

  
  strcpy(syncfile,"/tmp/rastroXXXXXX");
  if (mkstemp(syncfile) == -1)
  {
    ERROR(logger) << "could not make sync file " << syncfile << endl;
  }

  for_each(files_to_open.begin(), files_to_open.end(),
      [&](string & item)
      {
        if (! opened_files.count(item)) {
          INFO(logger) << "opening rastro " << item <<endl;
          fname = strdup(item.c_str());
          int ret = rst_open_file(fname, &data, syncfile, _RST_BUF_SIZE);
          free(fname);

          opened_files.insert(item);

          if (ret == -1) {
            ERROR(logger) << "trace file " << item << " could not be opened" << endl;
            exit(1);
          }
        } else {
          WARN(logger) << "double open " << item << endl;
        }
      }
    );


  // now loops on the events
  //


  vector<serve_entry_t>           evt_serve_list;
  vector<serve_entry_t>::iterator evt_serve_list_it;

  while (rst_decode_event(&data, &event)) {
    Paje::event_id_t evt_id = event.type;

    //cerr << "event type == " << event.type << " == " << evt_id << endl;
    int nevts = 0;

    pair<event_id_map_t::iterator,event_id_map_t::iterator> equal_range;
    event_id_map_t::iterator it;
    
    // load constant symbols into the tables
    (*(symbols[_GLOBAL_TABLE]))[Paje::idf1_name].set_value(event.id1);
    (*(symbols[_GLOBAL_TABLE]))[Paje::idf2_name].set_value(event.id2);

    if (symbols_per_file.count(make_pair(event.id1,event.id2))) {
      symbols[_LOCAL_TABLE] = symbols_per_file[make_pair(event.id1,event.id2)];
    } else {
      symbols[_LOCAL_TABLE] = new symbols_table_t();
      symbols_per_file[make_pair(event.id1,event.id2)] = symbols[_LOCAL_TABLE];
    }


    timestamp = (double) event.timestamp / 1000000 ;
   

    // get all events that must be served with this id
    equal_range = event_ids->equal_range(evt_id);

    // Get all events that are triggered by this id and their priorities
    for(it = equal_range.first;
        it != equal_range.second;
        ++it) 
    {
      serve_entry_t entry;

      ++nevts;
      entry.evt = it->second;

      if (debug)
        cerr << "=== Will serve id " << evt_id << " with " << entry.evt->name << endl;

      // Loads symbols from this evvent
      entry.evt->load_symbols(evt_id, &event,symbols[_GLOBAL_TABLE]);

      // loads some symbols from defined user maps
      entry.evt->get_symbols_from_map(evt_id, symbols[_GLOBAL_TABLE], symbols, usermaps);

      // get priority for this event
      entry.priority = entry.evt->get_priority(evt_id, timestamp, symbols);
      
      // push to do the output later
      evt_serve_list.push_back( entry );
    }



    // make the events list happen in the correct order
    sort(evt_serve_list.begin(), evt_serve_list.end(),
        [](serve_entry_t e1, serve_entry_t e2) {
          return (e1.priority > e2.priority);
        } );

    for(evt_serve_list_it = evt_serve_list.begin();
        evt_serve_list_it < evt_serve_list.end();
        ++evt_serve_list_it)
    {
      double dummy_prio;
      // we must repeat some work here -- different events triggered by the
      // same id MAY use conflicting symbol names...
      //
      //
      Paje::BaseEvent * evt = evt_serve_list_it->evt;
      
      if (debug)
        cerr << "=== Serving id " << evt_id << " with " << evt->name << endl;

      // Loads symbols from this evvent
      evt->load_symbols(evt_id, &event,symbols[_GLOBAL_TABLE]);

      // loads some symbols from defined user maps
      evt->get_symbols_from_map(evt_id, symbols[_GLOBAL_TABLE], symbols, usermaps);
      // saves predefined symbols in local table
      evt->push_symbols(evt_id, symbols[_GLOBAL_TABLE],symbols[_LOCAL_TABLE]);

      if (do_map) {
        // saves predefined symbols in map
        evt->map_symbols(evt_id, symbols, usermaps);
      }
     
      if (do_output) {
        // now generate paje output for this event
        evt->trigger(evt_id,
                          timestamp,
                          symbols,
                          &(dummy_prio),
                          out);
      }

    }
    

    // Just avoid warning of lost events on the mapping phase -- save the
    // space for map output
    if (do_output) {
      // warn once that there's no event for this id
      if ( (!nevts) && (!ignored_ids.count(evt_id)) ) {
        cerr << "Ignoring id " << evt_id  << endl;
        ignored_ids.insert(evt_id);
      }
    }

    evt_serve_list.clear();
  }


  // close f{{iles
  rst_close_file( &data );

  if (debug)
  {
    cerr << " === Symbols per file at the end of output === " << endl;
    for (auto it  = symbols_per_file.begin(); 
              it != symbols_per_file.end();
            ++it)
    {
      cerr << "Rastro (" << Paje::idf1_name << "=" << it->first.first 
           << " , " << Paje::idf2_name << "=" << it->first.second << ")" << endl;
        for (auto symbol_it  = it->second->begin();
                symbol_it != it->second->end();
              ++symbol_it)
      {
        cerr << "      " << symbol_it->first << " => ";
        symbol_it->second.format("",cerr);
        cerr << endl;
      }
    }



    cerr << " === Global Symbols === " << endl;
    for (auto symbol_it  = symbols[_GLOBAL_TABLE]->begin();
                symbol_it != symbols[_GLOBAL_TABLE]->end();
              ++symbol_it)
      {
        cerr << "      " << symbol_it->first << " => ";
        symbol_it->second.format("",cerr);
        cerr << endl;
      }
  }

  if (do_map || debug )
  {

    cerr << "=== Maps at the end of input ===" <<endl;
    for (auto map_it = usermaps.begin();
              map_it != usermaps.end();
              ++ map_it)
    {
      cerr << "Map " << map_it->first << endl;
      for (auto symbol_it = map_it->second.begin();
                symbol_it != map_it->second.end();
                ++symbol_it)
      {
        cerr << "   " << symbol_it->first << "  == >  " ;
        symbol_it->second.format("",cerr);
        cerr << endl;
      }
    }
  }



  return timestamp;
}




double  rastro_loop_events(list<string> &files_to_open, ostream &out, unsigned int n_maps, bool remap, bool debug)
{

  user_defined_maps_t usermaps;

  for (int i = 1; i <= n_maps; i++ )
  {
    // runs the mapping code, but generate no output
    INFO(logger) << "Mapping #" << i <<endl;
    _rastro_loop_events(files_to_open, out, usermaps, debug, true, false);
  }

  INFO(logger) << "Output" << endl;
  double timestamp = _rastro_loop_events(files_to_open, out, usermaps, debug,  (n_maps == 0) || remap, true);


  INFO(logger) << "Output ended after " << timestamp << " seconds" <<endl;


  return timestamp;
}
