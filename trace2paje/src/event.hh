// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 18:35:26 -0300 (kassick)"
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

#include "baseevent.hh"
#include "attributes.hh"
#include "semantic_types.hh"


#include <string>
#include <iostream>

namespace Paje {

  class Event: public BaseEvent 
  {
    public:
      Event(const string & name, attribs_t * attr);

      virtual bool do_start(double timestamp,
          symbols_table_t * symbols, 
          double * priotity,
          ostream &out);


      virtual bool do_end(double timestamp,
          symbols_table_t * symbols,
          double * priotity,
          ostream &out);

      virtual bool do_trigger(double timestamp,
          symbols_table_t * symbols,
          double * priority,
          ostream &out);

      virtual bool has_ids() const;

      virtual void gen_auto_ids(long int * base_id);

  };

}







#endif
