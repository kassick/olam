// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_loop.cc"
// Created: "Ter, 27 Set 2011 10:23:09 -0300 (kassick)"
// Updated: "Ter, 27 Set 2011 17:31:22 -0300 (kassick)"
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
  map<string,rst_file_t> file_to_data;
  char *syncfile;
  char * fname;
  
  strcpy(syncfile,"rastroXXXXXX");
  if (mkstemp(syncfile) == -1)
  {
    cerr << "could not make sync file " << syncfile << endl;
  }

  for_each(files_to_open.begin(), files_to_open.end(),
      [&](string & item)
      {
        if (! file_to_data.count(item)) {
          rst_file_t &data = file_to_data[item];
          fname = strdup(item.c_str());
          int ret = rst_open_file(fname, &data, syncfile, _RST_BUF_SIZE);
          free(fname);

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
  rst_event_t event;
  rst_file_t data;
  set<Paje::event_id_t> ignored_ids;
  symbols_table_t symbols;

  while (rst_decode_event(&data, &event)) {
    Paje::event_id_t evt_id = event.type;
    int nevts = 0;
    event_id_map_t::iterator it = event_ids->find(evt_id);

    // get every event in list
    while (it != event_ids->end()) {
      ++nevts;

      pair<Paje::event_id_t, Paje::Event*> p = *it;
      
      Paje::Event * evt = p.second;
      evt->load_symbols(&event,&symbols);
      double timestamp = (double) event.timestamp / 1000000 ;
      evt->trigger(evt_id, timestamp, &symbols, out);

      ++it;

    }

    // warn once that there's no event for this id
    if ( (!nevts) && (!ignored_ids.count(evt_id)) ) {
      cerr << "Ignoring id " << evt_id  << endl;
      ignored_ids.insert(evt_id);
    }
  }


  // close files
  for_each(file_to_data.begin(), file_to_data.end(),
      [&](pair<string, rst_file_t> p)
      {
        rst_close_file( & (p.second) );
      }
    );

  return;
}
