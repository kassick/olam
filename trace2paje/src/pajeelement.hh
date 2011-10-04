// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/pajeelement.hh"
// Created: "Ter, 04 Out 2011 14:15:57 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:36:03 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  pajeelement.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 14:15:57 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __PAJEELEMENT_H__
#define __PAJEELEMENT_H__

#include <sstream>
#include <string>
#include <iostream>

using namespace std;

namespace Paje {

  class PajeElement {
    public:
      virtual void do_header(ostream &header);
      //virtual void toPaje(stringstream &definitions, stringstream &header, stringstream &finalization);
  };
}


#endif
