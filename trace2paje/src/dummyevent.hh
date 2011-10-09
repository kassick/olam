// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/dummyevent.hh"
// Created: "Ter, 04 Out 2011 13:48:27 -0300 (kassick)"
// Updated: "Dom, 09 Out 2011 17:57:36 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  dummyevent.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 13:48:27 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __DUMMYEVENT_H__
#define __DUMMYEVENT_H__


#include "baseevent.hh"

#include <string>
#include <iostream>


namespace Paje {

  class DummyEvent: public BaseEvent {
      public: 
        DummyEvent(const string &name, BaseEventType * evt_type);
        virtual bool do_start(double timestamp,
            symbols_table_t * symbols,
            double * priority,
            ostream &out);

        virtual bool do_end(double timestamp,
            symbols_table_t * symbols,
            double * priority,
            ostream &out);

        virtual bool do_trigger(double timestamp,
            symbols_table_t * symbols,
            double * priority,
            ostream &out);

        virtual bool has_ids() const;

        virtual void gen_auto_ids(long int * base_id, set<event_id_t> & unique_ids);
  };

} // namespace paje


#endif
