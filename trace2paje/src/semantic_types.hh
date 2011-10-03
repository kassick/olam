// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantic_types.hh"
// Created: "Sex, 30 Set 2011 16:13:42 -0300 (kassick)"
// Updated: "Seg, 03 Out 2011 16:33:42 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  semantic_types.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  30-09-2011 16:13:42 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */



#ifndef __SEMANTIC_TYPES_H__
#define __SEMANTIC_TYPES_H__

#include <stdlib.h>
#include <queue>
#include <iostream>
#include <map>
#include "tree.hpp"
#include "container.hh"
#include "attributes.hh"
#include "event.hh"


using namespace Paje;
using namespace std;



typedef map<string,hierarchy_t *> container_type_names_t;

typedef std::map<std::string      , Paje::Event * >      event_name_map_t;
typedef std::multimap<Paje::event_id_t , Paje::Event * > event_id_map_t;
typedef std::map<std::string      , Paje::EventType * >  event_type_name_map_t;

extern attribs_t * attributes;

extern hierarchy_t * toplevel_hierarchy;
extern container_type_names_t * container_type_names;

extern event_type_name_map_t * eventtype_names;
extern event_name_map_t      * event_names;
extern event_id_map_t        * event_ids;
extern queue<string>         files_to_parse;

extern list<pair<string,Paje::Event*>> * ordered_event_names;


extern attribs_t * late_parse_tree, *early_parse_tree;

#endif
