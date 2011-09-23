// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.hh"
// Created: "Seg, 01 Ago 2011 15:46:26 -0300 (kassick)"
// Updated: "Sex, 23 Set 2011 19:21:35 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  semantics.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  01-08-2011 15:46:26 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 *
 */


#ifndef __SEMANTICS_HH_H__
#define __SEMANTICS_HH_H__


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



typedef TreeNode<Paje::Container*> hierarchy_t;

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


extern attribs_t * late_parse_tree, *early_parse_tree;

SemanticAttribute * new_semantic_attribute();


void init_desc_parser();
void attr_to_event_types(attribs_t * attribs);
void attr_to_link_types(attribs_t * attribs);
void attr_to_links(attribs_t * attribs);
void attr_to_states(attribs_t * attribs);
hierarchy_t * attr_to_container_hierarchy(attribs_t * attr, hierarchy_t *top);
void parse_late_tree();

void hierarchy_to_paje(ostream &out);
void event_types_to_paje(ostream &out);
void map_accept_attrs(attribs_t * attribs);
void check_events_have_type();


template <typename CharT, typename Traits>
basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& out, const SemanticAttribute& r)
{
 return out<< r.toString();
}


#endif
