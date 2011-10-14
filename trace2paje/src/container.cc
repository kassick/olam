// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.cc"
// Created: "Qua, 27 Jul 2011 11:07:19 -0300 (kassick)"
// Updated: "Sex, 14 Out 2011 15:44:24 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  container.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27-07-2011 11:07:19 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "tree.hpp"
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <string.h>
#include "container.hh"
#include "semantics.hh"
#include "pajeelement.hh"
#include "paje_functions.hh"



//*****
//Paje namespace classes

namespace Paje {

  void Container::do_header(ostream &out)
  {
    pajeDefineContainerType(typeName, parent->typeName, typeName,out);
  }
  
  const string Container::toString() const{
    stringstream s;

    s << "Container " << typeName << endl;
    s << "   " << "formatName: " << formatName<< endl;
    if (triggerParent)
      s << "   " << "Create on Parent"<< endl;
    else
      s << "   " << "createEvent: " << createEvent<< endl;

    if (destroyWithParent)
      s << "   " << "Destroy with Parent"<< endl;
    else
      s << "   " << "destroyEvent: " << destroyEvent<< endl;
    /*<< " accepting ";
    if (accept_all)
      s << "ALL EVENTS";
    else {
      s << "[";
      list<string>::const_iterator it = accept_list.begin();
      while (it != accept_list.end()) {
        s << *it;
        ++it;
        if (it != accept_list.end())
          s << ", ";
      }
      s << "]";

    } */
    return s.str();
  }
  
  Container::Container(string _typeName) {
    init(_typeName);
  }
  
  
  void Container::init(string _typeName)
  {
    triggerParent = destroyWithParent = true;
    this->typeName = _typeName;
    accept_all = false;


  }

  // here is to fill in  the object from the tree of attributes
  Container::Container(string _typeName,attribs_t *t) {
    //cout << "Creating from tree" <<endl;
    init(_typeName)  ;
    this->fill_from_attr(t);
  }

  void Container::fill_from_attr(attribs_t *attr_tree)
  {
    //For each node   of the subtree (non-descending into sub-containers)
    //fill in the f  ields from the attributes
    //ATTENTION: The walk function does not descend into children when cb
    //returns false -- hence here, stopping at the first non-toplevel
    //container makes sure it just reads the attributes of the current
    //container
    walk_tree_head_first(attr_tree, [&](attribs_t * n, int level)
        {
          SemanticAttribute* attr = n->getVal();
          //cerr << "Creating container " << _typeName << " " << level << endl;
          //cerr << *attr << endl;
          if ((attr-> id == ID_CONTAINER) && (level > 0))
            return true; // we're done
          this->fill_fields_from_cb(attr);
          return false; // go on
        } );


    // Fill in event types
    walk_tree(attr_tree,[&](SemanticAttribute * attr, int level) {
      if ((attr->id == ID_CONTAINER) && (level > 0))
        return true; // we're done

      if (attr->id == ID_EVENT_TYPE) {
        event_types.push_back(attr->vals.name);
      }

        /*
      if (attr->id == ID_ACCEPT_LIST) {
        //cout << "On container " << this->typeName <<endl;
        if (!strcmp("@ALL",attr->vals.identifier_name))
        {
            //cout << "Accept All\n";
          this->accept_all = true;
        } else {
          //cout << "Accept " << attr->vals.identifier_name << endl;
          accept_list.push_back(attr->vals.identifier_name);
          accept_all = false;
        }
      }
      */

      return false;

    });
  }

    void Container::fill_fields_from_cb(SemanticAttribute * attr) {
          switch (attr->id) {
            case ID_NAME:
              formatName = attr->vals.format_name;
              break;
            case ID_CREATE_EVENT:
              createEvent = attr->vals.create_event;
              triggerParent = false;
              break;
            case ID_CREATE_PARENT:
              triggerParent = true;
              break;
            case ID_DESTROY_WITH_PARENT:
              destroyWithParent = true;
              break;
            case ID_DESTROY_EVENT:
              destroyEvent = attr->vals.destroy_event;
              destroyWithParent = false;
              break;
          }
      }



}






