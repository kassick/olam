// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/linktype.cc"
// Created: "Ter, 04 Out 2011 13:44:56 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 19:33:42 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  linktype.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 13:44:56 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "attributes.hh"
#include "linktype.hh"
#include "container.hh"
#include "semantic_types.hh"
#include "paje_functions.hh"

#include <string>
#include <iostream>



/*******************************************************************************
 Paje::LinkType functions
*******************************************************************************/

Paje::LinkType::LinkType(const string& tn, Paje::Container * c, attribs_t * t) : 
  Paje::BaseEventType(tn,c)
{
  this->source = NULL;
  this->dest  = NULL;
  
  attribs_t::iterator it;
  for (it = t->begin(); it != t->end(); ++it)
  {
    SemanticAttribute * attr = (*it)->getVal();
    switch(attr->id) {
      case ID_LINK_SOURCE:
        if (!container_type_names->count(attr->vals.name)) {
          cerr << "Error: Can not find source container for link type " << typeName << endl;
          exit(1);
        }
        this->source = (*container_type_names)[attr->vals.name]->getVal();
        break;
      case ID_LINK_DEST:
        if (!container_type_names->count(attr->vals.name)) {
          cerr << "Error: Can not find source container for link type " << typeName << endl;
          exit(1);
        }
        this->dest = (*container_type_names)[attr->vals.name]->getVal();
        break;
      default:
        break;
    }
  }

  if ((this->source == NULL) || (this->dest == NULL) || (this->container == NULL)) {
    cerr << "Null fields, verify!" << endl;
    exit(1);
  }
}

void Paje::LinkType::do_header(ostream &out)
{
  pajeDefineLinkType(typeName, container->typeName,
                        source->typeName,
                        dest->typeName,
                        typeName,
                        out);
}

