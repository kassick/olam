// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.hh"
// Created: "Seg, 01 Ago 2011 15:46:26 -0300 (kassick)"
// Updated: "Sex, 30 Set 2011 16:16:52 -0300 (kassick)"
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
#include "semantic_types.hh"


using namespace Paje;
using namespace std;






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
