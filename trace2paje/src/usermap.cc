// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/usermap.cc"
// Created: "Dom, 13 Nov 2011 00:35:06 -0200 (kassick)"
// Updated: "Dom, 13 Nov 2011 01:02:59 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  usermap.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  13-11-2011 00:35:06 BRST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "usermap.hh"

void fill_entry_list(int map_id, map_entries_t &entries, attribs_t * attrs)
{
  map_entry_t entry;
  walk_tree_head_first(attrs,[&](attribs_t * n, int level) {
      SemanticAttribute * attr = n->getVal();
      entry.when = 0;

      if (map_id == attr->id) {
          walk_tree_head_first(n, [&entry](attribs_t * n1, int level1) {
            SemanticAttribute * attr1 = n1->getVal();
            switch (attr1->id) {
              case ID_MAP:
                return false;
                break;
              case ID_MAP_KEY:
                entry.key = attr1->vals.name;
                break;
              case ID_MAP_VALUE:
                entry.idf_name = attr1->vals.name;
                break;
              case ID_MAP_NAME:
                entry.map_name_format = attr1->vals.name;
                break;
              case ID_MAP_START:
                entry.when = USER_MAP_START;
                break;
              case ID_MAP_END:
                entry.when = USER_MAP_END;
                break;
              case ID_MAP_START_AND_END:
                entry.when = USER_MAP_BOTH;
                break;
            }
            return false;
          } );

          entries.push_back(entry);
          return true; // the inner walk_tree has already entered this subtree; prune
      }
    return false;
  });
}
