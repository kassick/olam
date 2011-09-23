// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/attributes.hh"
// Created: "Seg, 01 Ago 2011 16:11:04 -0300 (kassick)"
// Updated: "Sex, 23 Set 2011 18:04:49 -0300 (kassick)"
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
#include "rastro_helper.hh"

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
  ID_EVENT_TYPE,
  ID_EVENT_TYPE_DEF,
  ID_STATE,
  ID_STATE_TYPE,
  ID_PAJE_TYPENAME,
  ID_RASTRO_TYPE,
  ID_RASTRO_VALUE_NAME,
  ID_IDF,
  ID_EVENT_ID = 50,
  ID_EVENT_START,
  ID_STATE_START,
  ID_STATE_END,
  ID_LINK_TYPE,
  ID_LINK_SOURCE,
  ID_LINK_DEST,
  ID_LINK,
  ID_KEY_FORMAT,

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
      unsigned long event_id;
      rastro_basic_types_t rastro_type;
      rastro_basic_val_t rastro_val;
      //Paje::Container *container;
    } vals;


    const string toString() const;
};

typedef TreeNode<SemanticAttribute*> attribs_t;


#endif
