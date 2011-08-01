// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/paje.hh"
// Created: "Seg, 01 Ago 2011 15:34:40 -0300 (kassick)"
// Updated: "Seg, 01 Ago 2011 16:18:22 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  paje.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  01-08-2011 15:34:40 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef _PAJE_HH_H
#define _PAJE_HH_H

#include <string>
#include <sstream>

using namespace std;

namespace Paje {

  class PajeElement {
    public:
      void toPaje(stringstream &definitions, stringstream &header, stringstream &finalization);
  };
}

#endif
