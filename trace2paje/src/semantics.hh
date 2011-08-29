// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/semantics.hh"
// Created: "Seg, 01 Ago 2011 15:46:26 -0300 (kassick)"
// Updated: "Seg, 29 Ago 2011 19:12:54 -0300 (kassick)"
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
#include <iostream>
#include "tree.hpp"
#include "container.hh"
#include "attributes.hh"


using namespace Paje;
using namespace std;



typedef TreeNode<Paje::Container*> hierarchy_t;
typedef map<string,hierarchy_t *> container_ids_t;

extern attribs_t * attributes;

extern hierarchy_t * toplevel_hierarchy;
extern container_ids_t * container_ids;



SemanticAttribute * new_semantic_attribute();
hierarchy_t * attr_to_container_hierarchy(attribs_t * attr, hierarchy_t *top);


void init_desc_parser();
void check_unique_types();
void hierarchy_to_paje(ostream &out);


template <typename CharT, typename Traits>
basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& out, const SemanticAttribute& r)
{
 return out<< r.toString();
}


#endif
