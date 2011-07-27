// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/container.hh"
// Created: "Qua, 27 Jul 2011 11:08:49 -0300 (kassick)"
// Updated: "Qua, 27 Jul 2011 11:32:31 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  container.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  27-07-2011 11:08:49 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#include "tree.hh"
#include <iostream>
#include <map>
#include <string>


namespace Paje {

  using namespace std;

  class Container {

    public:
      string typeName;
      string formatName;
      bool triggerParent;
      string triggerEvent;
      void * eventTypes;

      string toPaje();

}; 




}
