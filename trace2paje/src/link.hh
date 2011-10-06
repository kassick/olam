// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/link.hh"
// Created: "Ter, 04 Out 2011 14:03:20 -0300 (kassick)"
// Updated: "Qui, 06 Out 2011 15:50:08 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  link.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 14:03:20 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __LINK_H__
#define __LINK_H__


#include "baseevent.hh"
#include "linktype.hh"
#include "attributes.hh"
#include <iostream>
#include <string>


namespace Paje {

  class Link: public BaseEvent {
    public:
      string format_key;
      Link(string &name, attribs_t * attribs);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);

      virtual bool do_end(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);
      
      virtual bool do_trigger(double timestamp,
          symbols_table_t * symbols, double * priority, ostream &out);
      
      virtual void fill_from_attr(attribs_t * attrs);
      virtual string toString() ;

      virtual bool has_ids() const;

      virtual void gen_auto_ids(long int * base_id);
      virtual bool fits_in_event_type(const BaseEventType* evt_type) const;
  };


}

void check_missing_links();

#endif
