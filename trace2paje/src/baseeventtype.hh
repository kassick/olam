// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/baseeventtype.hh"
// Created: "Ter, 04 Out 2011 12:03:12 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:22:06 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  baseeventtype.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 12:03:12 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __BASEEVENTTYPE_H__
#define __BASEEVENTTYPE_H__


#include "pajeelement.hh"
#include "container.hh"

#include <string>
#include <iostream>

namespace Paje {

  class BaseEventType: public PajeElement {
    public:
      string typeName;
      Container * container;

      BaseEventType(const string &typeName);
      BaseEventType(const string &typeName,Paje::Container * c);
      //EventType(string &typeName, attribs_t * attribs);

      virtual void do_header(ostream &out);

  };



} // namespace paje


#endif
