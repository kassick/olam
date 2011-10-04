// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/containertrigger.hh"
// Created: "Ter, 04 Out 2011 14:07:16 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:12:18 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  containertrigger.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 14:07:16 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __CONTAINERTRIGGER_H__
#define __CONTAINERTRIGGER_H__


#define CONTAINER_CREATE_PRIO (-10)
#define CONTAINER_DESTROY_PRIO (-1)

#include "baseevent.hh"
#include "container.hh"
#include <iostream>
#include <string>
#include <set>



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


#endif