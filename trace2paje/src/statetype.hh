// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/statetype.hh"
// Created: "Qui, 06 Out 2011 15:16:13 -0300 (kassick)"
// Updated: "Qui, 06 Out 2011 16:14:46 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  statetype.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  06-10-2011 15:16:13 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __STATETYPE_H__
#define __STATETYPE_H__



#include "baseeventtype.hh"
#include "container.hh"

#include <string>
#include <iostream>

namespace Paje {

  class StateType: public BaseEventType {
    public:

      StateType(const string &typeName);
      StateType(const string &typeName,Paje::Container * c);
      //EventType(string &typeName, attribs_t * attribs);

      virtual void do_header(ostream &out);

  };



} // namespace paje


#endif
