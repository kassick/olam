// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.cc"
// Created: "Qua, 27 Jul 2011 11:07:19 -0300 (kassick)"
// Updated: "Qua, 27 Jul 2011 11:32:39 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  container.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27-07-2011 11:07:19 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "tree.hh"
#include <iostream>
#include <map>
#include "container.hh"

namespace Paje {

  string Container::toPaje() {
    return "PAjeCreateContainer " + this->typeName;

  }

}



