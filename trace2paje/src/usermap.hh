// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/usermap.hh"
// Created: "Dom, 13 Nov 2011 00:27:12 -0200 (kassick)"
// Updated: "Dom, 13 Nov 2011 00:59:49 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  usermap.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  13-11-2011 00:27:12 BRST
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#include <string>
#include <map>
#include <vector>
#include "symbols.hh"
#include "attributes.hh"

#define USER_MAP_START (2<<1) 
#define USER_MAP_END   (2<<2)
#define USER_MAP_BOTH  (USER_MAP_START | USER_MAP_END)

typedef map<string, symbols_table_t> user_defined_maps_t;

typedef struct _map_entry_t {
  string map_name_format;
  string key;
  string idf_name;

  int when;
} map_entry_t;

typedef vector<map_entry_t> map_entries_t;

void fill_entry_list(int map_id, map_entries_t &entries, attribs_t * attr);
