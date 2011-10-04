// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/linktype.hh"
// Created: "Ter, 04 Out 2011 13:43:22 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 19:33:38 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  linktype.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 13:43:22 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __LINKTYPE_H__
#define __LINKTYPE_H__

#include "baseeventtype.hh"
#include "container.hh"
#include "attributes.hh"

#include <iostream>
#include <string.h>

namespace Paje {

  class LinkType: public BaseEventType {
    public:
      Container *source, *dest;

      LinkType(const string& typeName, Paje::Container * c, attribs_t * t);

      virtual void do_header(ostream &out);
  };
}



#endif
