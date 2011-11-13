// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_loop.cc"
// Created: "Ter, 27 Set 2011 10:23:09 -0300 (kassick)"
// Updated: "Dom, 13 Nov 2011 02:14:38 -0200 (kassick)"
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

  string output;

  Paje::BaseEvent * evt;

  struct _serve_entry_t & operator=(const struct _serve_entry_t &org) 
  {
    if (this == &org)
      return *this;

    this->priority = org.priority;
    this->evt = org.evt;
    this->output = org.output;

    return *this;
  }

} serve_entry_t;

double  rastro_loop_events(list<string> &files_to_open, ostream &out, bool debug)
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

  user_defined_maps_t usermaps;

  rst_event_t event;
  rst_file_t data;
  set<Paje::event_id_t> ignored_ids;

  map<pair<decltype(rst_event_t::id1),decltype(rst_event_t::id2)>, symbols_table_t *>    symbols_per_file;

  symbols[_END_TABLE] = NULL;
  // symbols[1] is local-to-file symbols,  
  symbols[_GLOBAL_TABLE] = new symbols_table_t(); // global symbols / dirty table

  
  strcpy(syncfile,"/tmp/rastroXXXXXX");
  if (mkstemp(syncfile) == -1)
  {
    cerr << "could not make sync file " << syncfile << endl;
  }

  for_each(files_to_open.begin(), files_to_open.end(),
      [&](string & item)
      {
        if (! opened_files.count(item)) {
          cerr << "opening rastro " << item <<endl;
          fname = strdup(item.c_str());
          int ret = rst_open_file(fname, &data, syncfile, _RST_BUF_SIZE);
          free(fname);

          opened_files.insert(item);

          if (ret == -1) {
            cerr << "Error: trace file " << item << " could not be opened" << endl;
            exit(1);
          }
        } else {
          cerr << "Warning:: double open " << item << endl;
        }
      }
    );


  // now loops on the events
  //


#if 0
  for_each(event_ids->begin(), event_ids->end(),
      [&](pair<Paje::event_id_t, Paje::BaseEvent*> p)
      {
        cerr << p.first << " => " << p.second->name << endl;
      }
      );
#endif



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

    /*
    cerr << "for id " << evt_id << " got the following limits: ";
    if (equal_range.first != event_ids->end())
       cerr << (equal_range.first)->second->name;
    else
      cerr << "(end)";
    cerr << " to ";

    if (equal_range.second != event_ids->end())
      cerr << (equal_range.second)->second->name;
    else
      cerr << "(end)";
    cerr << endl;
    */

    // Generate output strings for all events served with this id
    for(it = equal_range.first;
        it != equal_range.second;
        ++it) 
    {
      serve_entry_t entry;
      ostringstream tmp_out;

      // Serve the event
      //cerr << "serving " << it->second->name << endl;
      ++nevts;
      entry.evt = it->second;

      if (debug)
        cerr << "=== Serving id " << evt_id << " with " << entry.evt->name << endl;

      // Loads symbols from this evvent
      entry.evt->load_symbols(evt_id, &event,symbols[_GLOBAL_TABLE]);

      // loads some symbols from defined user maps
      entry.evt->get_symbols_from_map(evt_id, symbols[_GLOBAL_TABLE], symbols, usermaps);

      // saves predefined symbols in local table
      entry.evt->push_symbols(evt_id, symbols[_GLOBAL_TABLE],symbols[_LOCAL_TABLE]);

      // saves predefined symbols in map
      entry.evt->map_symbols(evt_id, symbols, usermaps);
     
      // now generate paje output for this event
      entry.evt->trigger(evt_id,
                        timestamp,
                        symbols,
                        &(entry.priority),
                        tmp_out);
      

      entry.output = tmp_out.str();


      // push to do the output later
      evt_serve_list.push_back( entry );





      //entry.output.clear();
    }


    sort(evt_serve_list.begin(), evt_serve_list.end(),
        [](serve_entry_t e1, serve_entry_t e2) {
          return (e1.priority > e2.priority);
        } );
    /*
          
          
          pair<string,Paje::BaseEvent*> p1  ,pair<string,Paje::BaseEvent*> p2) {
          Paje::BaseEvent * evt1 = p1.second;
          Paje::BaseEvent * evt2 = p2.second;
          return (evt1->get_priority(p1.first) > evt2->get_priority(p2.firat));
          }
        );*/

    for(evt_serve_list_it = evt_serve_list.begin();
        evt_serve_list_it < evt_serve_list.end();
        ++evt_serve_list_it)
    {
      out << evt_serve_list_it->output;

      /*
      pair<string,Paje::BaseEvent*> p = *it;
      
      Paje::BaseEvent * evt = p.second;
      //cerr << "11Event " << evt->name << " prio " << evt->get_priority() << endl;
      evt->load_symbols(evt_id, &event,&symbols);
      evt->trigger(evt_id, timestamp, &symbols, out); */
    }
    

    // warn once that there's no event for this id
    if ( (!nevts) && (!ignored_ids.count(evt_id)) ) {
      cerr << "Ignoring id " << evt_id  << endl;
      ignored_ids.insert(evt_id);
    }

    evt_serve_list.clear();
  }


  // close files
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
  }



  return timestamp;
}
