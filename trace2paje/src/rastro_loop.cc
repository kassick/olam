// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_loop.cc"
// Created: "Ter, 27 Set 2011 10:23:09 -0300 (kassick)"
// Updated: "Qua, 28 Set 2011 23:51:46 -0300 (kassick)"
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
#include <algorithm>
#include <map>
#include <set>
#include <iostream>


// librastro is a mess
extern "C" {
#include <rastro.h>
}


void  rastro_loop_events(list<string> &files_to_open, ostream &out)
{

  // Open files, use a tmp sync file
  char syncfile[256];
  char * fname;
  set<string> opened_files;
  
  rst_event_t event;
  rst_file_t data;
  set<Paje::event_id_t> ignored_ids;
  symbols_table_t symbols;
  
  strcpy(syncfile,"/tmp/rastroXXXXXX");
  if (mkstemp(syncfile) == -1)
  {
    cerr << "could not make sync file " << syncfile << endl;
  }

  for_each(files_to_open.begin(), files_to_open.end(),
      [&](string & item)
      {
        if (! opened_files.count(item)) {
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


  cerr << "Now loops!!" << endl;
  // now loops on the events
  //
  for_each(event_ids->begin(), event_ids->end(),
      [&](pair<Paje::event_id_t, Paje::Event*> p)
      {
        cerr << p.first << " => " << p.second->name << endl;
      }
      );



  vector<Paje::Event *> evt_serve_list;
  vector<Paje::Event *>::reverse_iterator evt_serve_list_it;
  while (rst_decode_event(&data, &event)) {
    Paje::event_id_t evt_id = event.type;
    cerr << "evemt type == " << event.type << " == " << evt_id << endl;
    int nevts = 0;

    pair<event_id_map_t::iterator,event_id_map_t::iterator> equal_range;
    event_id_map_t::iterator it;
    
    equal_range = event_ids->equal_range(evt_id);

    for(it = equal_range.first;
        it != equal_range.second;
        ++it) 
    {
      ++nevts;

      evt_serve_list.push_back(it->second);
      cerr << "Event " << it->second->name << endl;
    }


    sort(evt_serve_list.begin(), evt_serve_list.end(),
        [](Paje::Event* evt1, Paje::Event* evt2) {
          return evt1->timestamp_stack.top() < evt2->timestamp_stack.top();
        });

    for(evt_serve_list_it = evt_serve_list.rbegin();
        evt_serve_list_it < evt_serve_list.rend();
        ++evt_serve_list_it)
    {

      //pair<Paje::event_id_t, Paje::Event*> p = *it;
      
      Paje::Event * evt = *evt_serve_list_it; //p.second;
      cerr << "11Event " << evt->name << endl;
      evt->load_symbols(evt_id, &event,&symbols);
      symbols[Paje::idf1_name].set_value(event.id1);
      symbols[Paje::idf2_name].set_value(event.id2);
      double timestamp = (double) event.timestamp / 1000000 ;
      evt->trigger(evt_id, timestamp, &symbols, out);
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

  return;
}
