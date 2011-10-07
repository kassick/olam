// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/state.hh"
// Created: "Ter, 04 Out 2011 13:59:55 -0300 (kassick)"
// Updated: "Sex, 07 Out 2011 15:21:00 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  state.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 13:59:55 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __STATE_H__
#define __STATE_H__

#include "baseevent.hh"
#include "attributes.hh"
#include "paje_functions.hh"
#include <string>
#include <iostream>

namespace Paje {

  class State: public BaseEvent {
    public:

      State(string &_name, attribs_t * attrs) ;

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);
      
      virtual void fill_from_attr(attribs_t * attrs);

      virtual bool has_ids() const;

      virtual void gen_auto_ids(long int * base_id);
      virtual bool fits_in_event_type(const BaseEventType* evt_type) const;
  } ;


}// namespace

void check_opened_states(ostream & out);
int close_pending_states(double timestamp,string containerName, Paje::Container * c, ostream & out);

#endif

