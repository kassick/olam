// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:06:01 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  event.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  03-08-2011 16:14:50 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __EVENT_H__
#define __EVENT_H__

#include "attributes.hh"
#include "paje.hh"
#include "event.hh"
#include "symbols.hh"
#include "container.hh"
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <stack>
#include <list>
#include <iostream>
#include <limits>
#include "rastro_helper.hh"


#include "baseeventtype.hh"
#include "baseevent.hh"

extern "C" {
#include <rastro.h>
}


#include "baseevent.hh"
#include "linktype.hh"
#include "link.hh"
#include "state.hh"
#include "dummyevent.hh"


#define CONTAINER_CREATE_PRIO (-10)
#define CONTAINER_DESTROY_PRIO (-1)



namespace Paje {

  extern set<pair<string,string>> container_unique_names;







  class ContainerCreateTrigger: public BaseEvent {
    public:
      ContainerCreateTrigger(Paje::Container * c, hierarchy_t * n);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols, ostream &out);
      

      virtual void push_timestamp(double timestamp);

    protected:
      Container * container;
      hierarchy_t * hierarchy;

  };



} // Paje namespace





// sadly, #SHORT_NAME[0] does not work as the compiler does not recognize
// str_const[int_const] as a char in compilation time
#define CASE_TYPE1(CHAR_NAME,SHORT_NAME,VALUE) \
  case CHAR_NAME : \
    VALUE = SHORT_NAME;   \
    break;                \


#endif
