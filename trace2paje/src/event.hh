// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/event.hh"
// Created: "Qua, 03 Ago 2011 16:14:50 -0300 (kassick)"
// Updated: "Qua, 03 Ago 2011 16:17:38 -0300 (kassick)"
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

#include "paje.hh"
#include <sstream>
#include <string>

namespace Paje {
  class EventType {
    public:
      string typeName;
      Container * cont;
  };

  class Event {
    public:
      EventType *eventType;

  };
}

