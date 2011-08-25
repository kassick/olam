// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/attributes.hh"
// Created: "Seg, 01 Ago 2011 16:11:04 -0300 (kassick)"
// Updated: "Qua, 03 Ago 2011 16:09:13 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  attributes.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  01-08-2011 16:11:04 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef _ATTRIBUTES_HH_H
#define _ATTRIBUTES_HH_H

#include "tree.hpp"
using namespace std;

enum _attrib_ids {
  ID_NAME,
  ID_FORMAT_NAME,
  ID_CREATE_EVENT,
  ID_CREATE_PARENT,
  ID_DESTROY_EVENT,
  ID_DESTROY_CHILDREN,
  ID_CONTAINER,
  ID_ACCEPT_LIST,
  ID_IGNORE_LIST,

  ID_NOP
};

class SemanticAttribute {
  public:
    int id;
    union _atrib_vals {
      char *name;
      char *format_name;
      char *create_event;
      char *destroy_event;
      char *container_type;
      char *identifier_name;
      bool create_parent;
      bool destroy_children;
      //Paje::Container *container;
    } vals;


    const string toString() const;
};

typedef TreeNode<SemanticAttribute*> attribs_t;


#endif
