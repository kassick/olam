// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.hh"
// Created: "Qua, 27 Jul 2011 11:08:49 -0300 (kassick)"
// Updated: "Sex, 02 Set 2011 14:22:32 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  container.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  27-07-2011 11:08:49 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */



#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <iostream>
#include <map>
#include <string>
#include "attributes.hh"
#include "tree.hpp"

using namespace std;


namespace Paje {
  class Container;
}









namespace Paje {

  using namespace std;

  class Container {
    private:
      void fill_fields_from_cb(SemanticAttribute * attr);
      void init(string _typeName);

    public:

      Container(string);
      Container(string,attribs_t * attr_head);
      list<string> event_types;
      bool accept_all;
      string typeName;
      string formatName;
      string createEvent;
      string destroyEvent;
      bool triggerParent,destroyChildren;
      string triggerEvent;

      string toPaje();

      const string toString()const;

  }   ;



}


template <typename CharT, typename Traits>
basic_ostream<CharT, Traits>& operator<<(basic_ostream<CharT, Traits>& out, const Paje::Container& r)
{
 return out<< r.toString();
}



#endif
