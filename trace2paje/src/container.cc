// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.cc"
// Created: "Qua, 27 Jul 2011 11:07:19 -0300 (kassick)"
// Updated: "Sex, 29 Jul 2011 19:04:19 -0300 (kassick)"
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
#include "container.hh"

namespace Paje {

  string Container::toPaje() {
    return "PAjeCreateContainer " + this->typeName;

  } 
  
  const string Container::toString() const{
    stringstream s;

    s << "Container with name " << typeName << " formated to " <<formatName;
    return s.str();
  }
  
  Container::Container() {
    triggerParent = destroyChildren = false;
  }
  
  
  // here is to fill in  the object from the tree of attributes


    Container::Container(attribs_t *t) {
      cout << "Creating from tree" <<endl;
      walk_tree(t, [&](SemanticAttribute* attr, int level) 
          {
            if ((attr-> id == ID_CONTAINER) && (level > 0))
              return true; // we're done
            this->fill_fields_from_cb(attr);
            return false; // go on
          } );
    }

    void Container::fill_fields_from_cb(SemanticAttribute * attr) {
          switch (attr->id) {
            case ID_NAME:
              formatName = attr->vals.format_name;
              break;
            case ID_CREATE_EVENT:
              createEvent = attr->vals.create_event;
              break;
            case ID_CREATE_PARENT:
              triggerParent = attr->vals.create_parent;
              break;
            case ID_DESTROY_CHILDREN:
              destroyChildren = attr->vals.destroy_children;
              break;
            case ID_DESTROY_EVENT:
              destroyEvent = attr->vals.destroy_event;
              break;
          }
      }



}


//Globals
attribs_t *attributes;
hierarchy_t *toplevel_hierarchy;


const string SemanticAttribute::toString() const {
  stringstream s;
  switch (id) {
    case ID_NAME:
      s << "FormatName: " <<vals.format_name;
      break;
    case ID_CREATE_EVENT:
      s << "Create on Event "<<vals.create_event;
      break;
    case ID_CREATE_PARENT:
      if (vals.create_parent)
        s << "Create on Parent";
      else
        s << "Do not create on parent -- huh!?";
      break;
    case ID_DESTROY_CHILDREN:
      if (vals.destroy_children)
        s << "Destroy on Children";
      else
        s << "Do not destroy on children -- huh!?";
      break;
    case ID_DESTROY_EVENT:
      s << "Destroy on Event " << vals.destroy_event;
      break;
  }

  return s.str();
}

//Utility functions

SemanticAttribute *new_semantic_attribute()
{
  return new SemanticAttribute();
}




/*struct add_subtree_cb {
  hierarchy_t * h top = NULL;

  void operator()(attribs_t * t, int level) const {
    
  }
};*/


hierarchy_t * attr_to_container_hierarchy(attribs_t * attr, hierarchy_t *top)
{
  hierarchy_t * top_, *h;
  attribs_t::iterator it;

  h = NULL;

  if (!attr) return NULL;

  if (attr->getVal()->id == ID_CONTAINER) {
    top_ = new hierarchy_t(attr->getVal()->vals.container);
    top->addChild(top_);
    h = top_;
  } else {
    top_ = top;
  }

  for(it = attr->begin(); it != attr->end(); ++it)
  {
    hierarchy_t * child = attr_to_container_hierarchy(*it,top_);
    /*if (child) {
      top_->addChild(child);
    }*/
  }

  return h;

}

/*
void add_tree_to_hierarchy(attribs_t *attr)
{
  hierarchy_t * top = NULL;
  walk_tree( attr, [&](SemanticAttribute * attr, int level) {
        if (attr->id == ID_CONTAINER) {
          hierarchy_t * h = new hierarchy_t(attr->vals.container);
          if (top) {
            top->addChild(h);
            top = h;
          } else {
            top = h;
          }
        }
      }
    );

  toplevel_hierarchy->addChild(top);
  cout << "Container Hierarchy for this subtree is: "<<endl;
  print_tree(top);
}
    */
void reparent_containers(Paje::Container * c, attribs_t::iterator it)
{
}

void init_desc_parser()
{
  using namespace Paje;
  Paje::Container * zero;

  zero = new Container();
  zero->typeName = "0";

  toplevel_hierarchy = new TreeNode < Paje::Container * >(zero);
  attributes = new TreeNode < SemanticAttribute *>(NULL);
}

